"""Validate local input files and explicit statistical designs before queueing."""
from pathlib import Path
from urllib.parse import unquote, urlparse
import csv, re
import numpy as np
import pandas as pd

ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}')

def local_path(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('文件路径不能为空')
    value = value.strip()
    if value.startswith('file:'):
        value = unquote(urlparse(value).path).lstrip('/') if re.match(r'file:/+[A-Za-z]:', value) else unquote(urlparse(value).path)
    if re.match(r'^[A-Za-z]:[\\/]', value):
        value = '/mnt/' + value[0].lower() + '/' + value[3:].replace('\\', '/')
    return Path(value).expanduser().resolve()

def local_file(value):
    path = local_path(value)
    if not path.is_file():
        raise ValueError(f'文件不存在：{path}')
    return path

def read_table(path, **kwargs):
    # Many sequencing reports name plain tab-separated text "xls".
    return pd.read_csv(path, sep=None, engine='python', encoding='utf-8-sig', **kwargs)

def bounded(value, low, high, name, cast=float):
    try:
        value = cast(value)
    except (TypeError, ValueError):
        raise ValueError(f'{name} 必须为有效数字')
    if not np.isfinite(value) or not low <= value <= high:
        raise ValueError(f'{name} 超出范围 {low}–{high}')
    return value

def validate(cfg):
    if not isinstance(cfg, dict):
        raise ValueError('请求必须为 JSON 对象')
    cfg = dict(cfg)
    mode = cfg.get('mode', 'counts')
    if mode not in ('counts', 'fastq', 'salmon'):
        raise ValueError('未知输入类型')
    m = read_table(local_file(cfg.get('samples')), dtype=str).fillna('')
    if not {'sample', 'group'}.issubset(m.columns):
        raise ValueError('样本表需要 sample、group 列')
    if m['sample'].duplicated().any() or not all(ID.fullmatch(s) for s in m['sample']):
        raise ValueError('样本名须唯一且只含字母、数字、下划线、点或短横线')
    excluded = cfg.get('excluded_samples', [])
    if not isinstance(excluded, list) or not set(excluded).issubset(set(m['sample'])):
        raise ValueError('排除样本必须来自样本表')
    m = m[~m['sample'].isin(excluded)].copy()
    if len(m) < 4 or m.group.nunique() < 2 or min(m.groupby('group').size()) < 2:
        raise ValueError('至少两个分组，每组至少2个生物学重复，建议≥3')
    if (m.group == '').any():
        raise ValueError('分组不可为空')
    comparisons = cfg.get('comparisons') or [dict(name='comparison', control=cfg.get('control'), treatment=cfg.get('treatment'))]
    if not isinstance(comparisons, list) or not 1 <= len(comparisons) <= 30:
        raise ValueError('需要1–30个比较组合')
    names = set()
    for c in comparisons:
        if not isinstance(c, dict):
            raise ValueError('比较组合格式错误')
        c['name'] = c.get('name') or 'comparison'
        if not ID.fullmatch(c['name']) or c['name'] in names:
            raise ValueError('比较名须唯一且只含字母、数字、下划线、点或短横线')
        if not c.get('control') or not c.get('treatment') or c['control'] == c['treatment']:
            raise ValueError('请选择不同的对照组和处理组')
        if not {c['control'], c['treatment']}.issubset(set(m.group)):
            raise ValueError('比较中的分组不在样本表中')
        names.add(c['name'])
    covariates = cfg.get('covariates')
    if covariates is None:
        covariates = [dict(name='batch', type='categorical')] if 'batch' in m and m.batch.nunique() > 1 else []
    if not isinstance(covariates, list):
        raise ValueError('协变量必须为列表')
    matrix = [np.ones((len(m), 1))]
    cov_names = set()
    for c in covariates:
        if not isinstance(c, dict) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', str(c.get('name', ''))):
            raise ValueError('协变量列名格式错误')
        name = c['name']
        c.setdefault('type', 'categorical')
        if name in ('sample', 'group') or name in cov_names or name not in m or (m[name] == '').any():
            raise ValueError('协变量缺失、重复或列名无效：' + name)
        cov_names.add(name)
        if c.get('type', 'categorical') == 'numeric':
            values = pd.to_numeric(m[name], errors='raise').to_numpy(float)
            if not np.isfinite(values).all() or np.std(values) == 0:
                raise ValueError('数值协变量无变异或含非有限值：' + name)
            matrix.append(values[:, None])
        elif c.get('type', 'categorical') == 'categorical':
            if m[name].nunique() < 2:
                raise ValueError('分类协变量只有一个水平，请移除此协变量：' + name)
            matrix.append(pd.get_dummies(m[name], drop_first=True, dtype=float).to_numpy())
        else:
            raise ValueError('协变量类型只能为 categorical 或 numeric')
    matrix.append(pd.get_dummies(m.group, drop_first=True, dtype=float).to_numpy())
    design = np.concatenate(matrix, axis=1)
    if np.linalg.matrix_rank(design) < design.shape[1] or len(m) <= design.shape[1]:
        raise ValueError('批次/协变量与分组混杂，或没有残差自由度')
    cfg.update(mode=mode, comparisons=comparisons, covariates=covariates,
               control=comparisons[0]['control'], treatment=comparisons[0]['treatment'],
               threads=bounded(cfg.get('threads', 4), 1, 8, '线程', int),
               memory=bounded(cfg.get('memory', 24), 4, 28, '内存 GB', int),
               alpha=bounded(cfg.get('alpha', .05), 1e-6, .999, 'P阈值'),
               lfc=bounded(cfg.get('lfc', 1), 0, 10, 'LFC阈值'),
               permutations=bounded(cfg.get('permutations', 1000), 100, 10000, 'GSEA预估排列', int),
               min_count=bounded(cfg.get('min_count', 10), 1, 1000, '低表达计数阈值', int),
               min_samples=bounded(cfg.get('min_samples', 2), 1, len(m), '最低表达样本数', int),
               min_set=bounded(cfg.get('min_set', 15), 3, 100, '基因集最小大小', int),
               max_set=bounded(cfg.get('max_set', 500), 20, 5000, '基因集最大大小', int),
               gmt_id=cfg.get('gmt_id', 'ENSEMBL'), significance=cfg.get('significance', 'padj'),
               shrinkage=cfg.get('shrinkage', 'none'), pathway_activity=bool(cfg.get('pathway_activity', True)),
               wgcna=bool(cfg.get('wgcna', False)), traits=cfg.get('traits', []))
    if cfg['min_set'] >= cfg['max_set'] or cfg['significance'] not in ('padj', 'pvalue') or cfg['shrinkage'] not in ('none', 'normal'):
        raise ValueError('统计/基因集参数无效')
    if not isinstance(cfg['traits'], list):
        raise ValueError('表型列必须为列表')
    for trait in cfg['traits']:
        if trait not in m or trait in ('sample', 'subject', 'fastq1', 'fastq2', 'quant') or (m[trait] == '').any():
            raise ValueError('表型列缺失或无效：' + str(trait))
    if mode == 'counts':
        table = read_table(local_file(cfg.get('counts')), index_col=0, converters={0:str})
        if table.index.has_duplicates or table.index.isna().any() or (table.index.astype(str).str.strip()=='').any() or not set(m['sample']).issubset(table.columns):
            raise ValueError('计数矩阵基因须唯一，列须包含所有样本')
        values = table[m['sample']].to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError('DESeq2 需要非负整数原始计数；不能使用 TPM/FPKM')
        if (values.sum(axis=0) == 0).any():
            raise ValueError('样本总计数不能为0')
        if (values > np.iinfo(np.int32).max).any():
            raise ValueError('单基因原始计数超过 DESeq2 整数范围，请核对是否误用其他数据')
    else:
        mapping = read_table(local_file(cfg.get('tx2gene')), dtype=str)
        if mapping.shape[1] != 2 or mapping.isna().any().any() or mapping.iloc[:, 0].str.replace(r'\|.*$', '', regex=True).str.replace(r'\..*$', '', regex=True).duplicated().any():
            raise ValueError('tx2gene 需要带表头的两列，转录本（去版本号后）须唯一')
        if mode == 'fastq':
            local_file(cfg.get('transcripts'))
            if cfg.get('library', 'A') not in ['A', 'IU', 'ISF', 'ISR', 'U', 'SF', 'SR']:
                raise ValueError('不支持的文库链特异性类型')
            if 'fastq1' not in m:
                raise ValueError('FASTQ 模式需要 fastq1 列，双端添加 fastq2 列')
            for col in ('fastq1', 'fastq2'):
                if col in m:
                    if (m[col] == '').any():
                        raise ValueError('同一任务不能混合单端和双端')
                    m[col] = [str(local_file(f)) for f in m[col]]
            files=[path for col in ('fastq1','fastq2') if col in m for path in m[col]]
            if len(files)!=len(set(files)):
                raise ValueError('同一个 FASTQ 被重复使用；R1/R2 和不同样本必须使用各自文件，不能复制为生物学重复')
            cfg['preprocess'] = cfg.get('preprocess', 'fastp')
            if cfg['preprocess'] not in ('fastp', 'none'):
                raise ValueError('未知预处理方式')
            cfg['quality'] = bounded(cfg.get('quality', 20), 5, 35, 'fastp Q阈值', int)
            cfg['length'] = bounded(cfg.get('length', 30), 20, 250, 'fastp 最短长度', int)
            cfg['fragment_mean'] = bounded(cfg.get('fragment_mean',200),30,2000,'单端片段长度均值')
            cfg['fragment_sd'] = bounded(cfg.get('fragment_sd',80),1,1000,'单端片段长度标准差')
        else:
            if 'quant' not in m:
                raise ValueError('Salmon 导入需要 quant 列（quant.sf 路径）')
            m['quant'] = [str(local_file(f)) for f in m['quant']]
            if m['quant'].duplicated().any():
                raise ValueError('同一个 quant.sf 被重复当作多个样本，请核对生物学重复')
            reference_names=None
            mapping_ids=set(mapping.iloc[:,0].str.replace(r'\|.*$','',regex=True).str.replace(r'\..*$','',regex=True))
            for file in m['quant']:
                q = read_table(file,dtype={'Name':str})
                if not {'Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads'}.issubset(q.columns):
                    raise ValueError('quant.sf 缺少 Salmon 标准字段')
                numeric=q[['Length','EffectiveLength','TPM','NumReads']].apply(pd.to_numeric,errors='coerce')
                if q.Name.isna().any() or q.Name.duplicated().any() or not np.isfinite(numeric.to_numpy()).all() or (numeric[['TPM','NumReads']]<0).any().any() or (numeric[['Length','EffectiveLength']]<=0).any().any() or numeric.NumReads.sum()==0:
                    raise ValueError('Salmon 定量无效')
                if reference_names is not None and q.Name.tolist()!=reference_names:
                    raise ValueError('各 quant.sf 的转录本列表或顺序不同，请使用同一参考索引的完整量化文件')
                reference_names=q.Name.tolist()
                normalized=q.Name.str.replace(r'\|.*$','',regex=True).str.replace(r'\..*$','',regex=True)
                if normalized.duplicated().any() or not (set(normalized)&mapping_ids):
                    raise ValueError('Salmon 转录本与 tx2gene 无法唯一匹配，请核对版本与管道符前的 ID')
    if cfg.get('annotation'):
        ann = read_table(local_file(cfg['annotation']), dtype=str).fillna('')
        if not {'gene', 'symbol', 'entrez', 'go', 'ontology'}.issubset(ann.columns):
            raise ValueError('注释表需要 gene symbol entrez go ontology 五列')
        if not set(ann.ontology).issubset({'BP', 'CC', 'MF', ''}):
            raise ValueError('GO ontology 只能为 BP/CC/MF')
        known = ann[ann.go != '']
        if not known.go.str.fullmatch(r'GO:\d{7}').all():
            raise ValueError('GO编号须为 GO:七位数字')
        if not known.ontology.isin(['BP','CC','MF']).all() or (ann.gene=='').any():
            raise ValueError('注释基因不能为空，每个 GO 必须标明 BP/CC/MF')
        if (known.groupby('go').ontology.nunique()>1).any():
            raise ValueError('同一个 GO 编号不能同时属于不同 ontology')
    if cfg.get('gmt'):
        if cfg['gmt_id'] not in ('ENSEMBL', 'ENTREZID', 'SYMBOL'):
            raise ValueError('GMT ID 类型无效')
        terms = set()
        for line in local_file(cfg['gmt']).read_text().splitlines():
            parts = line.split('\t')
            if len(parts) < 3 or parts[0] in terms:
                raise ValueError('GMT 每行至少3列，通路名不能重复')
            terms.add(parts[0])
        if not terms:
            raise ValueError('GMT 文件为空')
    if cfg.get('tf_network'):
        net=read_table(local_file(cfg['tf_network']))
        if not {'source','target','weight'}.issubset(net):raise ValueError('调控网络需要 source target weight 三列')
        if net[['source','target']].isna().any().any() or net[['source','target']].duplicated().any():raise ValueError('调控网络的 source/target 缺失或重复')
        weight=pd.to_numeric(net.weight,errors='raise')
        if not np.isfinite(weight).all() or (weight==0).any():raise ValueError('调控网络权重必须为有限非零数')
        if cfg.get('tf_id','SYMBOL') not in ('ENSEMBL','ENTREZID','SYMBOL'):raise ValueError('调控网络ID类型无效')
        if cfg.get('tf_species','Mus musculus')!='Mus musculus':raise ValueError('本工作站仅接受已转换至小鼠的调控网络')
    cfg['tf_id']=cfg.get('tf_id','SYMBOL')
    for key in ('tf_network','counts', 'transcripts', 'tx2gene', 'gmt', 'annotation', 'samples'):
        if cfg.get(key):
            cfg[key] = str(local_file(cfg[key]))
    cfg['deg_profile']=cfg.get('deg_profile','strict')
    if cfg['deg_profile'] not in ('nominal','strict'):raise ValueError('差异方案只能为 nominal 或 strict')
    cfg.update(alpha=.05,lfc=0 if cfg['deg_profile']=='nominal' else 1,significance='pvalue' if cfg['deg_profile']=='nominal' else 'padj')
    cfg['wgcna_genes'] = bounded(cfg.get('wgcna_genes', 3000), 100, 6000, 'WGCNA基因数', int)
    cfg['wgcna_power'] = bounded(cfg.get('wgcna_power', 6), 1, 30, 'WGCNA软阈值', int)
    cfg['wgcna_power_mode'] = cfg.get('wgcna_power_mode', 'manual')
    if cfg['wgcna_power_mode'] not in ('auto', 'manual'):
        raise ValueError('WGCNA 软阈值模式无效')
    if mode == 'fastq' and cfg.get('genome'):
        cfg['genome'] = str(local_file(cfg['genome']))
    if mode == 'fastq':
        paired = 'fastq2' in m
        if (paired and cfg.get('library', 'A') in ('U', 'SF', 'SR')) or (not paired and cfg.get('library', 'A') in ('IU', 'ISF', 'ISR')):
            raise ValueError('文库类型与单端/双端 FASTQ 不一致，无法判断时使用自动 A')
    return cfg, m
