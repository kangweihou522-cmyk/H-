"""Validate standalone downstream analyses without pretending they used FASTQ."""
import re
import numpy as np
import pandas as pd
from validation import ID,read_table,local_file,bounded,validate_resources

STAGES=('enrichment','expression','activity','network')


def validate_module(request):
    cfg=dict(request);stage=cfg['analysis_stage'];m=pd.DataFrame(columns=['sample','group'])
    cfg.update(mode='module',analysis_stage=stage)
    if stage!='enrichment' or cfg.get('expression'):
        m=read_table(local_file(cfg.get('samples')),dtype=str).fillna('')
        if not {'sample','group'}.issubset(m):raise ValueError('样本表需要 sample、group')
        if m['sample'].duplicated().any() or not all(ID.fullmatch(x) for x in m['sample']):raise ValueError('样本名须唯一且格式有效')
        excluded=cfg.get('excluded_samples',[])
        if not isinstance(excluded,list) or not set(excluded).issubset(set(m['sample'])):raise ValueError('排除样本不在样本表中')
        m=m[~m['sample'].isin(excluded)].copy()
        if len(m)<4 or (m.group=='').any():raise ValueError('表达矩阵单步分析至少4个样本，group不能为空')
        matrix=read_table(local_file(cfg.get('expression')),index_col=0,converters={0:str})
        if matrix.index.has_duplicates or matrix.index.isna().any() or (matrix.index.astype(str).str.strip()=='').any() or not set(m['sample']).issubset(matrix.columns):raise ValueError('表达矩阵基因须唯一，列须包含全部样本')
        values=matrix[m['sample']].apply(pd.to_numeric,errors='coerce').to_numpy(float)
        if len(matrix)<10 or not np.isfinite(values).all() or np.sum(np.var(values,axis=1)>1e-12)<10:raise ValueError('需要至少10个有变异的基因，表达矩阵须为有限数值')
        if cfg.get('expression_unit','vst') not in ('vst','log2_normalized','log2_tpm'):raise ValueError('仅接受已转换的VST/log2表达，不接受原始TPM或计数直接当作VST')
        cfg['expression_unit']=cfg.get('expression_unit','vst')
        cfg['samples']=str(local_file(cfg['samples']));cfg['expression']=str(local_file(cfg['expression']))
    else:
        cfg['samples']='';cfg['expression']=''
    if stage=='enrichment':
        result=read_table(local_file(cfg.get('differential_input')),dtype={'gene':str})
        if not {'gene','log2FoldChange','pvalue'}.issubset(result):raise ValueError('富集入口需要完整差异结果：gene、log2FoldChange、pvalue；可提供padj，GSEA另需stat')
        if len(result)<10 or result.gene.isna().any() or result.gene.duplicated().any():raise ValueError('请导入包含未显著基因的完整差异结果，gene须唯一，至少10行')
        for key in ('log2FoldChange','pvalue','padj','stat','baseMean'):
            if key not in result:continue
            x=pd.to_numeric(result[key],errors='coerce')
            if (result[key].notna()&x.isna()).any() or np.isinf(x).any():raise ValueError('差异结果数值列无效：'+key)
            if key in ('pvalue','padj') and ((x<0)|(x>1)).any():raise ValueError('P/padj须在0–1之间或为空')
        if result.pvalue.notna().sum()<2:raise ValueError('完整差异结果没有足够的可用P值')
        cfg['differential_input']=str(local_file(cfg['differential_input']))
        cfg['run_gsea']=bool(cfg.get('run_gsea',True) and 'stat' in result and pd.to_numeric(result['stat'],errors='coerce').notna().sum()>=10)
    covariates=cfg.get('covariates')
    if covariates is None:covariates=[dict(name='batch',type='categorical')] if 'batch' in m and m.batch.nunique()>1 else []
    if not isinstance(covariates,list):raise ValueError('协变量须为列表')
    cfg['covariates']=covariates if len(m) else []
    columns=[np.ones((len(m),1))];names=set()
    for item in cfg['covariates']:
        name=item.get('name');kind=item.get('type','categorical');item['type']=kind
        if not isinstance(name,str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name) or name in names or name in ('sample','group') or name not in m or (m[name]=='').any():raise ValueError('协变量缺失、重复或无效')
        names.add(name)
        if kind=='numeric':
            x=pd.to_numeric(m[name],errors='coerce').to_numpy(float)
            if not np.isfinite(x).all() or np.std(x)==0:raise ValueError('数值协变量缺失或无变异')
            columns.append(x[:,None])
        elif kind=='categorical':
            if m[name].nunique()<2:raise ValueError('分类协变量只有一个水平')
            columns.append(pd.get_dummies(m[name],drop_first=True,dtype=float).to_numpy())
        else:raise ValueError('协变量类型无效')
    if len(m):columns.append(pd.get_dummies(m.group,drop_first=True,dtype=float).to_numpy())
    comparisons=cfg.get('comparisons') or [dict(name='comparison',control=cfg.get('control'),treatment=cfg.get('treatment'))]
    if stage in ('expression','network') or (stage=='activity' and m.group.nunique()<2):comparisons=[]
    if not isinstance(comparisons,list) or len(comparisons)>30:raise ValueError('比较组合无效')
    seen=set()
    for cmp in comparisons:
        name=cmp.get('name','comparison');cmp['name']=name
        if not isinstance(name,str) or not ID.fullmatch(name) or name in seen:raise ValueError('比较名须有效且唯一')
        seen.add(name)
        if not cmp.get('control') or not cmp.get('treatment') or cmp['control']==cmp['treatment']:raise ValueError('请声明原始比较方向：不同的对照组/处理组')
        if len(m) and not {cmp['control'],cmp['treatment']}.issubset(set(m.group)):raise ValueError('比较中的分组不在样本表中')
    if stage=='activity' and comparisons:
        design=np.concatenate(columns,axis=1)
        if min(m.groupby('group').size())<2 or np.linalg.matrix_rank(design)<design.shape[1] or len(m)<=design.shape[1]:raise ValueError('通路活性比较需要组内重复、满秩设计与残差自由度')
    if stage=='enrichment' and len(comparisons)!=1:raise ValueError('独立富集每次导入一个比较的完整差异结果')
    cfg['comparisons']=comparisons
    cfg['control']=comparisons[0]['control'] if comparisons else m.group.iloc[0] if len(m) else ''
    cfg['treatment']=comparisons[0]['treatment'] if comparisons else ''
    cfg['traits']=cfg.get('traits',[])
    if not isinstance(cfg['traits'],list):raise ValueError('表型列须为列表')
    for trait in cfg['traits']:
        if trait in ('sample','subject','fastq1','fastq2','quant') or trait not in m or (m[trait]=='').any():raise ValueError('表型列缺失或无效')
    for key,default,lo,hi,integer in [('threads',4,1,8,True),('memory',24,4,28,True),('min_set',15,3,100,True),('max_set',500,20,5000,True),('permutations',1000,100,10000,True),('wgcna_genes',3000,100,6000,True),('wgcna_power',6,1,30,True)]:
        cfg[key]=bounded(cfg.get(key,default),lo,hi,key,int if integer else float)
    if cfg['min_set']>=cfg['max_set']:raise ValueError('基因集大小范围无效')
    cfg['deg_profile']=cfg.get('deg_profile','strict')
    if cfg['deg_profile'] not in ('strict','nominal'):raise ValueError('差异筛选方案无效')
    cfg.update(alpha=.05,lfc=0 if cfg['deg_profile']=='nominal' else 1,significance='pvalue' if cfg['deg_profile']=='nominal' else 'padj',gmt_id=cfg.get('gmt_id','ENSEMBL'),tf_id=cfg.get('tf_id','SYMBOL'),preprocess='none',shrinkage='none',min_count=10,min_samples=2)
    cfg['pathway_activity']=stage=='activity' and bool(cfg.get('pathway_activity',True))
    cfg['wgcna']=stage=='network'
    cfg['wgcna_power_mode']=cfg.get('wgcna_power_mode','auto')
    if cfg['wgcna_power_mode'] not in ('auto','manual'):raise ValueError('软阈值模式无效')
    if stage!='activity':cfg['tf_network']=''
    if stage=='activity' and not cfg.get('gmt') and not cfg.get('tf_network'):raise ValueError('通路活性需GMT，或提供小鼠调控网络计算ULM')
    validate_resources(cfg)
    for key in ('annotation','gmt','tf_network'):
        if cfg.get(key):cfg[key]=str(local_file(cfg[key]))
    return cfg,m
