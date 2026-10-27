"""One validated plot specification for offline exports and interactive exploration."""
from pathlib import Path
import html
import hashlib
import json
import os
import re
import textwrap

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, is_color_like
from matplotlib.patches import Ellipse
from adjustText import adjust_text
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import chi2
import networkx as nx
import numpy as np
import pandas as pd
from dataio import read_frame

PALETTE = ['#126657','#b64c4c','#317ba0','#bd9a4b','#8668ab','#c96d91','#73913b','#77858c']
LABELS = dict(base_content='A/T/G/C 含量曲线', quality='碱基质量曲线', fastp_summary='原始/过滤 reads', qc_rates='Q20 / Q30 / GC 与保留率', mapping='Salmon 映射率', library_size='文库大小', distribution='表达分布', pca='PCA 样本诊断', pca_scree='PCA 方差解释', correlation='样本相关性热图', sample_distance='样本距离热图', dispersion='离散度诊断', volcano='火山图', ma='MA 图', p_hist='P 值分布诊断', heatmap='差异基因聚类热图', gene_expression='基因表达', go='GO 富集', gsea='GO GSEA', pathway='KEGG / GMT 富集', pathway_gsea='KEGG / GMT GSEA', gsea_curve='GSEA 富集曲线', leading_edge='GSEA 核心基因表达热图', enrichment_network='通路基因重叠网络', tf_activity='调控因子活性热图', tf_de='差异调控活性', activity='通路活性热图', activity_de='差异通路活性', pathway_traits='通路–表型关联热图', gene_trait='基因–表型散点图', wgcna_soft='WGCNA 软阈值诊断', wgcna_traits='模块–表型关联', wgcna_modules='共表达模块', wgcna_hubs='模块 hub 基因', wgcna_network='共表达子网络', deg_summary='各比较双方案差异统计', deg_overlap='跨比较 DEG 重叠',wgcna_enrichment='共表达模块功能富集')
TITLES = dict(zip(LABELS, ['Base content','Read quality','Sequencing QC','Sequencing QC rates','Salmon mapping rate','Library size','Expression distribution','Sample PCA','PCA explained variance','Sample correlation','Sample distance','DESeq2 dispersion','Differential expression','MA plot','P-value distribution','Differential gene heatmap','Gene expression','GO over-representation','GO GSEA','GMT over-representation','GMT GSEA','GSEA enrichment curve','Leading-edge expression','Gene-set overlap network','Regulatory activity (ULM)','Differential regulatory activity','Pathway activity (ssGSEA)','Differential pathway activity','Pathway / trait association','Gene / phenotype exploration','WGCNA soft-threshold diagnostics','Module / trait association','WGCNA modules','Module hub genes','Coexpression subnetwork','Differential gene overview','DEG overlap across contrasts','Module over-representation']))
NOTES = {
    'volcano':'点表示基因。悬停显示 ID、Symbol、LFC、P、padj；点击查看表达，输入基因可定位，标注框可拖动。改变图中阈值仅影响显示。',
    'heatmap':'默认逐行 Z-score；颜色表示该基因在不同样本中的相对高低，不能比较不同基因的绝对表达。分组色带跟随聚类后的样本顺序。',
    'pca':'最高变异 500 个 VST 基因，不自动去批次。椭圆仅为样本协方差的描述性 95% 区域，≥3 个非共线样本才显示，不能作为组间显著性检验。',
    'p_hist':'展示所有有可用 P 值的保留基因。独立过滤可能使 padj 缺失；低表达、离群点和模型失配均需结合其他诊断判断。',
    'leading_edge':'使用所选 GSEA 基因集的 leading-edge 成员，映射回保留的基因并显示 VST。可关闭逐行标准化。',
    'enrichment_network':'边表示基因集合的 Jaccard 重叠，不是通路调控方向或蛋白互作。默认使用显著 GSEA leading-edge；可调整阈值与成员定义。',
    'gene_trait':'原始散点用于检查样本和离群点；如有调整后结果，同时列出控制分组及全部协变量后的偏 Pearson r、P 与全局 BH padj。两者可能不同，均不是因果证据。',
    'wgcna_network':'从已保存的共表达边中显示有限子网络；边权为 signed adjacency，不是蛋白互作，不代表完整网络。grey 为未分配基因，不列为 hub。',
    'wgcna_traits':'原始 Pearson 关联；单元格可同时显示相关系数与表内 BH padj。它未调整协变量，与调整后基因–表型结果不同。',
    'wgcna_enrichment':'以进入网络的全部基因为背景（再按数据库/ontology覆盖限定），分别检验非grey模块。BH覆盖同数据库/ontology的所有模块×术语，不是DEG富集。',
    'activity':'ssGSEA 分数用于同一通路的跨样本比较；默认逐行标准化，不能把不同通路的分数当作统一活性单位。',
    'tf_activity':'ULM 根据导入的小鼠调控网络和靶基因表达估计活性，不是转录因子自身的 RNA 丰度。',
    'deg_overlap':'集合来自每个比较当前选择的 DEG 方案；先确认各比较方向及筛选方案。相交不代表独立重复验证。',
}


def module_colors(modules, neutral):
    # Unknown R/X11 names need a stable fallback across plots and repeated rows.
    return [neutral if n=='grey' else n if is_color_like(n) else PALETTE[int(hashlib.sha256(str(n).encode()).hexdigest()[:8],16)%len(PALETTE)] for n in modules]


def data_root(p):
    return p.parent.parent if p.parent.name=='contrasts' else p


def get_settings(p, kind, opts=None):
    if kind not in LABELS:
        raise ValueError('未知图表')
    root = data_root(p)
    analysis = json.loads((root/'config.json').read_text()) if (root/'config.json').exists() else {}
    f = p/(kind+'.settings.json')
    saved = json.loads(f.read_text()) if f.exists() else {}
    # Migrate the former universal 9x6 default; preserve other custom dimensions.
    if 'style_version' not in saved and saved.get('width')==9 and saved.get('height')==6 and kind in ('volcano','pca'):
        saved.pop('width'); saved.pop('height')
    width,height = (7,6.5) if kind in ('volcano','pca','enrichment_network','wgcna_network') else (9,6.5)
    cfg = dict(style_version=3,title='',width=width,height=height,font=11,dpi=300,up='#c56361',down='#397d9e',neutral='#c3cac6',alpha=.05,lfc=analysis.get('lfc',1),top=20,xlabel='',ylabel='',direction='all',significant_only=True,cluster=True,cluster_rows=True,cluster_cols=True,labels=5,genes='',gene='',term='',database='go',pcx=1,pcy=2,sample='',expression_scale='vst',gene_scale='normalized',gene_labels='symbol',annotated_genes='',label_positions={},point_size=28,point_alpha=.8,label_font=10,grid=True,legend=True,x_min='',x_max='',y_min='',y_max='',p_column=analysis.get('significance','padj'),group_colors=analysis.get('group_colors',{}),sample_labels=True,ellipse=False,symmetric_x=True,threshold_lines=True,show_points=True,connect_pairs=False,plot_style='box',heatmap_scale='row_z',color_map='blue_red',color_min='',color_max='',show_values=False,sample_order='original',ontology='all',enrichment_style='dot',enrichment_x='gene_ratio',sort_by='padj',wrap=36,read_end='R1',read_phase='both',module='',trait='',network_threshold=.2,network_nodes=40,network_members='leading_edge',correlation_method='pearson',profile_display='both')
    if (p/'active-profile.json').exists():
        active = json.loads((p/'active-profile.json').read_text())
        cfg.update(p_column=active['column'], alpha=active['alpha'], lfc=active['lfc'])
    if not (root/'normalized.csv').exists() and (root/'vst.csv').exists():cfg['gene_scale']='vst'
    cfg.update(saved); cfg.update(opts or {})
    cfg['style_version']=3
    ranges = [('width',3,20),('height',3,20),('font',6,30),('dpi',72,600),('top',3,100),('alpha',1e-6,.999),('lfc',0,10),('labels',0,30),('pcx',1,5),('pcy',1,5),('point_size',5,150),('point_alpha',.05,1),('label_font',6,22),('wrap',12,80),('network_threshold',.01,1),('network_nodes',5,120)]
    for key,lo,hi in ranges:
        try: cfg[key] = float(cfg[key])
        except (ValueError,TypeError): raise ValueError('图表参数不是数字：'+key)
        if not np.isfinite(cfg[key]) or not lo<=cfg[key]<=hi:
            raise ValueError('图表参数超出范围：'+key)
    for key in ('up','down','neutral'):
        if not re.fullmatch('#[0-9a-fA-F]{6}',str(cfg[key])):
            raise ValueError('颜色需要 #RRGGBB')
    for key in ('title','xlabel','ylabel','gene','term','sample','genes','annotated_genes','x_min','x_max','y_min','y_max','color_min','color_max','module','trait'):
        if not isinstance(cfg[key],str) or len(cfg[key])>(5000 if key in ('genes','annotated_genes') else 180):
            raise ValueError('图表文本无效/过长：'+key)
    for lo,hi in [('x_min','x_max'),('y_min','y_max'),('color_min','color_max')]:
        for k in (lo,hi):
            if cfg[k]!='' and not np.isfinite(float(cfg[k])): raise ValueError('范围必须为有限数')
        if cfg[lo]!='' and cfg[hi]!='' and float(cfg[lo])>=float(cfg[hi]): raise ValueError('范围最小值必须小于最大值')
    choices = dict(expression_scale=('vst','normalized','tpm'),gene_scale=('normalized','vst','tpm'),gene_labels=('symbol','id'),p_column=('padj','pvalue'),direction=('all','up','down'),database=('go','pathway'),plot_style=('box','violin'),heatmap_scale=('row_z','none'),color_map=('blue_red','purple_green','viridis','custom'),sample_order=('original','group'),ontology=('all','BP','CC','MF'),enrichment_style=('dot','bar','lollipop'),enrichment_x=('gene_ratio','rich_factor','overlap','neglog10'),sort_by=('padj','NES','overlap','gene_ratio'),read_end=('R1','R2'),read_phase=('before','after','both'),network_members=('leading_edge','set'),correlation_method=('pearson','spearman'),profile_display=('both','active'))
    for key, options in choices.items():
        if cfg[key] not in options: raise ValueError('图表选择无效：'+key)
    for key in ('cluster','cluster_rows','cluster_cols','significant_only','grid','legend','sample_labels','ellipse','symmetric_x','threshold_lines','show_points','connect_pairs','show_values'):
        if not isinstance(cfg[key], bool): raise ValueError('图表开关必须为布尔值：'+key)
    if not isinstance(cfg['group_colors'],dict) or any(not re.fullmatch('#[0-9a-fA-F]{6}',str(x)) for x in cfg['group_colors'].values()): raise ValueError('分组颜色格式无效')
    if not isinstance(cfg['label_positions'],dict) or len(cfg['label_positions'])>200: raise ValueError('标签位置无效')
    for name,pos in cfg['label_positions'].items():
        if not isinstance(pos,dict) or any(not isinstance(pos.get(k,0),(float,int)) or not np.isfinite(pos.get(k,0)) or abs(pos.get(k,0))>600 for k in ('ax','ay')): raise ValueError('标签位置超出范围')
    return cfg


def chart_catalog(p):
    charts = []
    def add(kind,stage,contrast=None): charts.append(dict(kind=kind,title=LABELS[kind],stage=stage,contrast=contrast))
    sources = dict(fastp_summary='fastp-qc.csv',qc_rates='fastp-qc.csv',quality='fastp-qc.csv',base_content='fastp-qc.csv',mapping='salmon-qc.csv',library_size='sample-qc.csv',distribution='vst.csv',pca='vst.csv',pca_scree='pca-variance.csv',correlation='sample-correlation.csv',sample_distance='sample-distance.csv',dispersion='dispersion.csv',deg_summary='analysis-summary.json',activity='pathway-activity.csv',tf_activity='tf-activity.csv',pathway_traits='pathway-trait.csv',wgcna_soft='wgcna-soft-threshold.csv',wgcna_traits='wgcna-module-trait.csv',wgcna_modules='wgcna-modules.csv',wgcna_hubs='wgcna-modules.csv',wgcna_network='wgcna-network.csv',wgcna_enrichment='wgcna-module-enrichment.csv')
    for kind,file in sources.items():
        if kind=='deg_summary' and (p/file).exists():
            contrasts=json.loads((p/file).read_text()).get('contrasts') or {}
            if not isinstance(contrasts,dict) or not any('up' in c for c in contrasts.values()):continue
        if (p/file).exists(): add(kind,'qc' if kind in ('fastp_summary','qc_rates','quality','base_content','mapping') else 'mining' if kind in ('activity','tf_activity','pathway_traits') or kind.startswith('wgcna') else 'expression')
    targets = [(q.name,q) for q in sorted((p/'contrasts').iterdir()) if q.is_dir()] if (p/'contrasts').is_dir() else [(None,p)]
    if (p/'vst.csv').exists() and not (p/'differential.csv').exists() and not (p/'contrasts').exists():add('gene_expression','expression')
    if len([q for _,q in targets if (q/'deg-all.csv').exists()])>1: add('deg_overlap','differential')
    for name,q in targets:
        for kind,file in dict(volcano='differential.csv',ma='differential.csv',p_hist='differential.csv',heatmap='differential.csv',gene_expression='differential.csv',go='go.csv',gsea='gsea.csv',pathway='pathway.csv',pathway_gsea='pathway_gsea.csv',activity_de='pathway-activity-differential.csv',tf_de='tf-activity-differential.csv').items():
            if not (q/file).exists(): continue
            if kind in ('heatmap','gene_expression') and not (p/'vst.csv').exists():continue
            if kind=='ma' and ('baseMean' not in read_frame(q/file) or not (read_frame(q/file).baseMean>0).any()):continue
            if kind in ('go','gsea','pathway','pathway_gsea','activity_de','tf_de') and read_frame(q/file).empty: continue
            add(kind,'differential' if kind in ('volcano','ma','p_hist','heatmap','gene_expression') else 'mining' if kind in ('activity_de','tf_de') else 'enrichment',name)
        valid_gsea = any((q/f).exists() and not read_frame(q/f).dropna(subset=['padj','NES']).empty for f in ('gsea.csv','pathway_gsea.csv'))
        if (p/'gene-sets.json').exists() and valid_gsea:
            for kind in ('gsea_curve','leading_edge','enrichment_network'): add(kind,'enrichment',name)
        if (p/'gene-trait-adjusted.csv').exists(): add('gene_trait','mining',name)
    return charts


def colors_for(m,root,cfg):
    custom = cfg['group_colors']
    return {g:custom.get(g,PALETTE[i%len(PALETTE)]) for i,g in enumerate(m.group.unique())}


def gene_labels(root,genes,cfg):
    if cfg['gene_labels']=='symbol' and (root/'gene-annotation.csv').exists():
        ann = read_frame(root/'gene-annotation.csv').drop_duplicates('gene').set_index('gene')
        return [str(ann.loc[g,'symbol']) if g in ann.index and pd.notna(ann.loc[g,'symbol']) else str(g) for g in genes]
    return list(genes)


def ordered_matrix(x,m,cfg,scale=True):
    x = x.copy()
    if scale and cfg['heatmap_scale']=='row_z': x = x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1).replace(0,1).fillna(1),axis=0)
    if cfg['sample_order']=='group' and set(x.columns).issubset(set(m['sample'])):
        order = m.sort_values('group',kind='stable')['sample'];x = x[[s for s in order if s in x.columns]]
    if cfg['cluster'] and not x.empty:
        if cfg['cluster_rows'] and len(x)>1: x = x.iloc[leaves_list(linkage(x.to_numpy(),method='average'))]
        if cfg['cluster_cols'] and x.shape[1]>1: x = x.iloc[:,leaves_list(linkage(x.T.to_numpy(),method='average'))]
    return x


def resolve_genes(root,r,queries):
    from results import gene_data
    genes = []
    for name in re.split(r'[,;\n]+', queries):
        if not name.strip(): continue
        try: gene = gene_data(root,name.strip(),r)['gene']
        except ValueError: continue
        if gene not in genes: genes.append(gene)
    if not genes: raise ValueError('指定的基因未在本次保留基因中找到')
    return genes


def gsea_database(p,cfg):
    database = cfg['database']
    filename = 'gsea.csv' if database=='go' else 'pathway_gsea.csv'
    if not cfg['term'] and (not (p/filename).exists() or read_frame(p/filename).empty): database = 'pathway' if database=='go' else 'go'
    return database


def load_plot(p,kind,cfg):
    from results import gene_data, curve_data
    root = data_root(p);m = read_frame(root/'samples.csv');v = read_frame(root/'vst.csv',index_col=0) if (root/'vst.csv').exists() else pd.DataFrame()
    r = read_frame(p/'differential.csv') if (p/'differential.csv').exists() else None
    contrast = p.name if p!=root else None
    if kind in ('volcano','ma'):
        x = r.dropna(subset=['log2FoldChange',cfg['p_column']]).copy()
        x = x[np.isfinite(x.log2FoldChange)&np.isfinite(x[cfg['p_column']])]
        sign = x[cfg['p_column']]<cfg['alpha']
        x['direction_plot'] = np.where(sign&(x.log2FoldChange>cfg['lfc']),'Up',np.where(sign&(x.log2FoldChange<-cfg['lfc']),'Down','Other'))
        x['color'] = x.direction_plot.map({'Up':cfg['up'],'Down':cfg['down'],'Other':cfg['neutral']})
        x['plot_x'] = x.log2FoldChange if kind=='volcano' else x.baseMean
        if kind=='ma': x = x[x.plot_x>0]
        manual = {g.strip().casefold() for g in re.split(r'[,;\n]+',cfg['annotated_genes']) if g.strip()}
        x['manual_label'] = x.gene.astype(str).str.casefold().isin(manual) | (x.symbol.fillna('').astype(str).str.casefold().isin(manual) if 'symbol' in x else False)
        x['plot_y'] = -np.log10(x[cfg['p_column']].clip(lower=1e-300)) if kind=='volcano' else x.log2FoldChange
        return 'scatter',x,('log2 fold change','−log10 BH adjusted P' if cfg['p_column']=='padj' else '−log10 unadjusted P') if kind=='volcano' else ('Mean normalized expression','log2 fold change')
    if kind=='pca':
        if (root/'pca-scores.csv').exists():
            x = read_frame(root/'pca-scores.csv');variance = read_frame(root/'pca-variance.csv').variance.to_numpy()
        else:
            # Same top-variable-gene convention as DESeq2 plotPCA for legacy jobs.
            genes = v.var(axis=1).nlargest(min(500,len(v))).index
            a = v.loc[genes,m['sample']].T.to_numpy(dtype=float,copy=True);a -= a.mean(axis=0)
            u,s,_ = np.linalg.svd(a,full_matrices=False);variance = s*s/(s*s).sum()
            x = pd.DataFrame({'sample':m['sample'],'group':m.group,**{f'PC{i+1}':u[:,i]*s[i] for i in range(min(5,u.shape[1]))}})
        a,b = int(cfg['pcx']),int(cfg['pcy'])
        if f'PC{a}' not in x or f'PC{b}' not in x or a==b: raise ValueError('请选择存在且不同的主成分')
        x['plot_x']=x[f'PC{a}'];x['plot_y']=x[f'PC{b}']
        return 'pca',x,(f'PC{a} ({variance[a-1]:.1%})',f'PC{b} ({variance[b-1]:.1%})')
    if kind in ('heatmap','leading_edge'):
        if cfg['genes']: genes = resolve_genes(root,contrast,cfg['genes'])
        elif kind=='leading_edge':
            database = gsea_database(p,cfg);curve = curve_data(root,database,cfg['term'],contrast)
            ann = read_frame(root/'gene-annotation.csv',dtype=str).fillna('')
            analysis = json.loads((root/'config.json').read_text())
            key = 'id' if database=='go' else {'ENSEMBL':'id','SYMBOL':'symbol','ENTREZID':'entrez'}[analysis['gmt_id']]
            genes = list(ann.loc[ann[key].isin(curve['leading']),'gene'])[:int(cfg['top'])]
        else:
            significant = r.significant if 'significant' in r else (r[cfg['p_column']]<cfg['alpha'])&(r.log2FoldChange.abs()>cfg['lfc'])
            genes = list(r[significant].nsmallest(int(cfg['top']),cfg['p_column']).gene)
        x = v.loc[[g for g in dict.fromkeys(genes) if g in v.index],m['sample']]
        return 'heatmap',ordered_matrix(x,m,cfg),('Samples','Row Z-score (VST)' if cfg['heatmap_scale']=='row_z' else 'VST expression')
    if kind in ('correlation','sample_distance'):
        x = v[m['sample']].corr(method=cfg['correlation_method']) if kind=='correlation' else read_frame(root/'sample-distance.csv',index_col=0).loc[m['sample'],m['sample']]
        if cfg['cluster'] and len(x)>1:
            distance = (1-x).clip(lower=0) if kind=='correlation' else x
            condensed_source=distance.to_numpy(copy=True)
            np.fill_diagonal(condensed_source,0)
            order = leaves_list(linkage(squareform(condensed_source,checks=False),method='average'));x = x.iloc[order,order]
        return 'matrix',x,('Samples','Samples')
    if kind in ('activity','tf_activity'):
        x = read_frame(root/('tf-activity.csv' if kind=='tf_activity' else 'pathway-activity.csv'),index_col=0)
        x = x.loc[x.var(axis=1).nlargest(int(cfg['top'])).index,m['sample']]
        return 'heatmap',ordered_matrix(x,m,cfg),('Samples','Row Z-score' if cfg['heatmap_scale']=='row_z' else 'ULM score' if kind=='tf_activity' else 'ssGSEA score')
    if kind in ('wgcna_traits','pathway_traits'):
        x = read_frame(root/('wgcna-module-trait.csv' if kind=='wgcna_traits' else 'pathway-trait.csv'))
        key = 'module' if kind=='wgcna_traits' else 'term'
        if kind=='pathway_traits':
            terms = x.groupby(key).padj.min().nsmallest(int(cfg['top'])).index;x = x[x[key].isin(terms)]
        return 'association',x,('Traits','Modules' if key=='module' else 'Pathways')
    if kind in ('go','pathway','gsea','pathway_gsea','activity_de','tf_de','wgcna_enrichment'):
        filename = {'activity_de':'pathway-activity-differential.csv','tf_de':'tf-activity-differential.csv'}.get(kind,kind+'.csv')
        x = read_frame(root/'wgcna-module-enrichment.csv') if kind=='wgcna_enrichment' else read_frame(p/filename)
        if kind=='wgcna_enrichment':
            x=x[x.database==cfg['database']]
            if cfg['module']: x=x[x.module==cfg['module']]
            if cfg['ontology']!='all': x=x[x.ontology==cfg['ontology']]
            x=x.copy();x['description']=x.module+' · '+x.description
        if kind in ('go','pathway') and 'direction' in x:
            x = x[x.direction==cfg['direction']]
            if cfg['ontology']!='all': x=x[x.ontology==cfg['ontology']]
        if kind in ('activity_de','tf_de'): x = x.rename(columns={'adj.P.Val':'padj','logFC':'NES'})
        x = x.dropna(subset=['padj'])
        if cfg['significant_only']: x = x[x.padj<cfg['alpha']]
        sort = cfg['sort_by'] if cfg['sort_by'] in x else 'padj'
        x = x.sort_values(sort,ascending=sort=='padj',kind='stable').head(int(cfg['top'])).copy()
        if 'description' not in x: x['description']=x.term
        if kind in ('go','pathway','wgcna_enrichment'):
            x['plot_x'] = -np.log10(x.padj.clip(lower=1e-300)) if cfg['enrichment_x']=='neglog10' else x[cfg['enrichment_x']]
            axis = {'gene_ratio':'Gene ratio (k / n)','rich_factor':'Rich factor (k / K)','overlap':'Overlap genes','neglog10':'−log10 BH adjusted P'}[cfg['enrichment_x']]
            return 'dot',x,(axis,'Gene sets')
        return 'barh',x,('Activity score difference' if kind in ('activity_de','tf_de') else 'Normalized enrichment score (NES)','Gene sets')
    if kind=='gsea_curve':
        return 'curve',curve_data(root,gsea_database(p,cfg),cfg['term'],contrast),('Gene rank','Running enrichment score')
    if kind in ('gene_expression','gene_trait'):
        candidates = r.dropna(subset=['padj']) if r is not None and 'padj' in r else pd.DataFrame()
        if candidates.empty and r is not None:candidates=r.sort_values('baseMean',ascending=False)
        gene = cfg['gene'] or (str(candidates.iloc[0].gene) if len(candidates) else str(v.mean(axis=1).idxmax()))
        data = gene_data(root,gene,contrast)
        scale = {'vst':'vst','normalized':'normalized','tpm':'TPM'}[cfg['gene_scale']]
        if scale not in data['samples'][0]: raise ValueError('当前任务没有所选表达标度：'+cfg['gene_scale'])
        data['scale'] = scale
        if kind=='gene_trait':
            traits = read_frame(root/'gene-trait-adjusted.csv')
            trait = cfg['trait'] or str(traits.sort_values('padj').iloc[0].trait)
            if trait not in m or pd.to_numeric(m[trait],errors='coerce').isna().any(): raise ValueError('请选择存在的连续数值表型')
            data['trait']=trait
            data['adjusted']=json.loads(traits[(traits.gene==data['gene'])&(traits.trait==trait)].to_json(orient='records'))
            return 'trait',data,(trait,{'vst':'VST','normalized':'Normalized count','TPM':'TPM'}[scale])
        return 'gene',data,('Group',{'vst':'VST expression','normalized':'Normalized count','TPM':'TPM'}[scale])
    if kind in ('enrichment_network','wgcna_network'):
        return 'network',network_data(p,kind,cfg),('','')
    if kind=='distribution':
        if cfg['expression_scale']=='vst': x=v;label='VST expression'
        else:
            path = root/('tpm.csv' if cfg['expression_scale']=='tpm' else 'normalized.csv')
            if not path.exists(): raise ValueError('当前任务没有所选表达标度；可选择VST/已导入转换表达。原始计数不能直接推算TPM。')
            x=np.log2(read_frame(path,index_col=0)+1);label='log2(TPM + 1)' if cfg['expression_scale']=='tpm' else 'log2(normalized count + 1)'
        return 'box',x[m['sample']],('Samples',label)
    if kind in ('quality','base_content'):
        sample = cfg['sample'] or m['sample'].iloc[0];f=root/'qc'/sample/'fastp.json'
        if not f.exists(): raise ValueError('该样本无 fastp 数据')
        return kind,json.loads(f.read_text()),('Position in read','Mean Phred quality' if kind=='quality' else 'Base fraction')
    if kind=='deg_summary':
        contrasts=json.loads((root/'analysis-summary.json').read_text())['contrasts'];rows=[]
        for name,c in contrasts.items():
            profiles = c.get('profiles',{}) if cfg['profile_display']=='both' else {c.get('active_profile','active'):c}
            for profile,spec in profiles.items(): rows.append(dict(name=name+' / '+profile,up=spec['up'],down=spec['down']))
        return 'deg',pd.DataFrame(rows),('Comparison / profile','Significant genes')
    if kind=='deg_overlap':
        sets={q.name:set(read_frame(q/'deg-all.csv').gene) for q in sorted((root/'contrasts').iterdir()) if (q/'deg-all.csv').exists()}
        x=pd.DataFrame(index=sets,columns=sets,dtype=float)
        for a in sets:
            for b in sets: x.loc[a,b]=len(sets[a]&sets[b])/len(sets[a]|sets[b]) if sets[a]|sets[b] else 0
        return 'matrix',x,('Comparisons','DEG Jaccard overlap')
    if kind=='wgcna_modules':
        return 'modules',read_frame(root/'wgcna-modules.csv').groupby('module').size().rename('genes').reset_index(),('Module','Genes')
    if kind=='wgcna_hubs':
        x=read_frame(root/'wgcna-modules.csv');x=x[(x.module!='grey')&x.kME.notna()]
        if cfg['module']: x=x[x.module==cfg['module']]
        x=x.reindex(x.kME.abs().sort_values(ascending=False).index).head(int(cfg['top']))
        return 'hubs',x,('Signed module membership (kME)','Genes')
    source=dict(library_size='sample-qc.csv',mapping='salmon-qc.csv',fastp_summary='fastp-qc.csv',qc_rates='fastp-qc.csv',dispersion='dispersion.csv',pca_scree='pca-variance.csv',wgcna_soft='wgcna-soft-threshold.csv',p_hist=None)
    if kind in source:
        x=r.copy() if kind=='p_hist' else read_frame(root/source[kind])
        axes=dict(library_size=('Sample','Estimated fragments / counts'),mapping=('Sample','Mapped fragments (%)'),fastp_summary=('Sample','Reads (R1 + R2)'),qc_rates=('Sample','Percent (%)'),dispersion=('Mean normalized expression','Dispersion'),pca_scree=('Principal component','Explained variance (%)'),wgcna_soft=('Power','Signed scale-free fit R²'),p_hist=('P value','Number of genes'))[kind]
        return kind,x,axes
    raise ValueError('未知图表')


def network_data(p,kind,cfg):
    root=data_root(p);graph=nx.Graph()
    if kind=='wgcna_network':
        edges=read_frame(root/'wgcna-network.csv');modules=read_frame(root/'wgcna-modules.csv').set_index('gene')
        eligible=modules.index[(modules.module!='grey')&modules.kME.notna()]
        if cfg['module']: eligible=modules.index[(modules.module==cfg['module'])&modules.kME.notna()]
        eligible=modules.loc[eligible].kME.abs().nlargest(int(cfg['network_nodes'])).index
        edges=edges[edges.source.isin(eligible)&edges.target.isin(eligible)&(edges.weight>=cfg['network_threshold'])]
        labels=dict(zip(modules.index,gene_labels(root,modules.index,cfg)))
        for row in edges.itertuples():
            graph.add_edge(row.source,row.target,weight=float(row.weight))
        for name in graph:
            graph.nodes[name].update(label=labels[name],value=float(modules.loc[name,'kME']),module=str(modules.loc[name,'module']))
    else:
        database=gsea_database(p,cfg);filename='gsea.csv' if database=='go' else 'pathway_gsea.csv'
        frame=read_frame(p/filename).dropna(subset=['padj','NES'])
        if cfg['significant_only']: frame=frame[frame.padj<cfg['alpha']]
        if cfg['ontology']!='all': frame=frame[frame.ontology==cfg['ontology']]
        frame=frame.sort_values('padj',kind='stable').head(min(int(cfg['top']),int(cfg['network_nodes'])))
        sets=json.loads((root/'gene-sets.json').read_text()).get(database,{})
        # Restrict full sets to genes actually present in the tested ranked universe.
        ranks=read_frame(p/('gsea-ranks.csv' if database=='go' else 'pathway_gsea-ranks.csv'),dtype={'gene':str})
        universe=set(ranks.gene)
        members={}
        for row in frame.itertuples():
            raw=str(row.leadingEdge).split(';') if cfg['network_members']=='leading_edge' and pd.notna(row.leadingEdge) else sets.get(row.term,[])
            if isinstance(raw,str): raw=[raw]
            members[row.term]=set(raw)&universe
            graph.add_node(row.term,label=str(row.description),value=float(row.NES),padj=float(row.padj))
        for i,a in enumerate(members):
            for b in list(members)[i+1:]:
                union=members[a]|members[b];overlap=len(members[a]&members[b]);weight=overlap/len(union) if union else 0
                if weight>=cfg['network_threshold']: graph.add_edge(a,b,weight=weight,overlap=overlap)
    pos=nx.spring_layout(graph,seed=42,weight='weight') if len(graph) else {}
    nodes=[dict(id=name,x=float(pos[name][0]),y=float(pos[name][1]),**attrs) for name,attrs in graph.nodes(data=True)]
    edges=[dict(source=a,target=b,**attrs) for a,b,attrs in graph.edges(data=True)]
    return dict(nodes=nodes,edges=edges,definition='Signed adjacency; saved limited subnetwork' if kind=='wgcna_network' else 'Jaccard '+cfg['network_members'])


def curve_series(data,cfg,kind):
    traces=[];ends=['read1' if cfg['read_end']=='R1' else 'read2']
    phases=['before_filtering','after_filtering'] if cfg['read_phase']=='both' else [cfg['read_phase']+'_filtering']
    for end in ends:
        for phase in phases:
            record=data.get(end+'_'+phase,{})
            if kind=='quality':
                values=record.get('quality_curves',{}).get('mean',[])
                if values: traces.append((cfg['read_end']+' '+phase.split('_')[0],values,cfg['neutral'] if phase.startswith('before') else cfg['down'],'solid'))
            else:
                for base,color in [('A','#126657'),('T','#b64c4c'),('G','#bd9a4b'),('C','#317ba0')]:
                    values=record.get('content_curves',{}).get(base,[])
                    if values: traces.append((base+' '+phase.split('_')[0],values,color,'dash' if phase.startswith('before') else 'solid'))
    if not traces: raise ValueError('所选 R1/R2 或过滤阶段没有对应 fastp 曲线')
    return traces


def heatmap_view(data,kind,root,cfg):
    annotations=None
    if kind in ('wgcna_traits','pathway_traits'):
        key='module' if kind=='wgcna_traits' else 'term'
        matrix=data.pivot(index=key,columns='trait',values='correlation')
        padj=data.pivot(index=key,columns='trait',values='padj').reindex_like(matrix)
        annotations=[[f'{matrix.iloc[i,j]:.2f}\nq={padj.iloc[i,j]:.2g}' for j in range(len(matrix.columns))] for i in range(len(matrix))]
    else: matrix=data
    labels=gene_labels(root,matrix.index,cfg) if kind in ('heatmap','leading_edge') else list(matrix.index)
    if annotations is None: annotations=[[f'{val:.2f}' for val in row] for row in matrix.to_numpy()]
    lo,hi=(None,None)
    if kind in ('heatmap','leading_edge','activity','tf_activity') and cfg['heatmap_scale']=='row_z': lo,hi=-2,2
    if kind in ('correlation','deg_overlap'): lo,hi=(-1,1) if kind=='correlation' else (0,1)
    if kind in ('wgcna_traits','pathway_traits'): lo,hi=-1,1
    if cfg['color_min']!='': lo=float(cfg['color_min'])
    if cfg['color_max']!='': hi=float(cfg['color_max'])
    return matrix,labels,annotations,lo,hi


def palette(cfg):
    if cfg['color_map']=='viridis': return plt.get_cmap('viridis'), 'Viridis'
    colors=['#6671a5','#ffffff','#58a58b'] if cfg['color_map']=='purple_green' else [cfg['down'],'#ffffff',cfg['up']]
    return LinearSegmentedColormap.from_list('workstation',colors), [[0,colors[0]],[.5,colors[1]],[1,colors[2]]]


def ellipse_points(frame):
    if len(frame)<3: return None
    values=frame[['plot_x','plot_y']].to_numpy(float);cov=np.cov(values,rowvar=False)
    if np.linalg.matrix_rank(cov)<2: return None
    w,u=np.linalg.eigh(cov);theta=np.linspace(0,2*np.pi,100)
    return np.array([np.cos(theta),np.sin(theta)]).T@np.diag(np.sqrt(w*chi2.ppf(.95,2)))@u.T+values.mean(axis=0)


def wrapped(values,cfg,web=False):
    return [('<br>' if web else '\n').join(html.escape(s) if web else s for s in textwrap.wrap(str(v),width=int(cfg['wrap']))) for v in values]


def selected_labels(data,cfg):
    significant=data[data.direction_plot!='Other']
    return pd.concat([data[data.manual_label],significant.nsmallest(int(cfg['labels']),cfg['p_column'])]).drop_duplicates('gene')


def render(p,kind,opts):
    cfg=get_settings(p,kind,opts);root=data_root(p);m=read_frame(root/'samples.csv');groups=colors_for(m,root,cfg)
    if kind in ('gsea_curve','leading_edge','enrichment_network'):cfg['database']=gsea_database(p,cfg)
    category,data,axes=load_plot(p,kind,cfg)
    axes=tuple(label.replace('VST',expression_label(root)) for label in axes)
    available={f.name for f in font_manager.fontManager.ttflist}
    family=next((x for x in ['Noto Sans CJK SC','Noto Sans CJK JP','Microsoft YaHei','WenQuanYi Zen Hei'] if x in available),'DejaVu Sans')
    with plt.rc_context({'font.family':family,'font.size':cfg['font'],'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False}):
        fig,ax=plt.subplots(figsize=(cfg['width'],cfg['height']),layout='constrained')
        try:
            if category=='scatter':
                for direction,color in [('Other',cfg['neutral']),('Down',cfg['down']),('Up',cfg['up'])]:
                    x=data[data.direction_plot==direction]
                    ax.scatter(x.plot_x,x.plot_y,c=color,s=cfg['point_size'],alpha=cfg['point_alpha'],edgecolors='#64716a',linewidths=.45,label=f'{direction} ({len(x)})')
                if kind=='volcano':
                    if cfg['threshold_lines']:
                        ax.axhline(-np.log10(cfg['alpha']),ls='--',color='#909a94',lw=.8)
                        for val in sorted({-cfg['lfc'],cfg['lfc']}): ax.axvline(val,ls='--',color='#909a94',lw=.8)
                    if cfg['symmetric_x'] and len(data):
                        extent=max(1,float(data.plot_x.abs().max())*1.12);ax.set_xlim(-extent,extent)
                    ax.set_ylim(bottom=0)
                else:
                    ax.set_xscale('log');ax.axhline(0,color='#909a94',ls='--',lw=.8)
                texts=[]
                for row in selected_labels(data,cfg).itertuples():
                    label=str(row.symbol) if hasattr(row,'symbol') and pd.notna(row.symbol) else row.gene
                    pos=cfg['label_positions'].get(row.gene)
                    box=dict(boxstyle='round,pad=.2',facecolor='white',edgecolor='#8b938e',alpha=.95)
                    if pos: ax.annotate(label,(row.plot_x,row.plot_y),xytext=(pos.get('ax',25)*.75,-pos.get('ay',-25)*.75),textcoords='offset points',fontsize=cfg['label_font'],fontstyle='italic',bbox=box,arrowprops=dict(arrowstyle='-',color='#69756d',lw=.7))
                    else: texts.append(ax.text(row.plot_x,row.plot_y,label,fontsize=cfg['label_font'],fontstyle='italic',bbox=box))
                if texts:
                    np.random.seed(42);adjust_text(texts,x=data.plot_x.to_numpy(),y=data.plot_y.to_numpy(),ax=ax,time_lim=.5,arrowprops=dict(arrowstyle='-',color='#69756d',lw=.7))
            elif category=='pca':
                for g in data.group.unique():
                    x=data[data.group==g];ax.scatter(x.plot_x,x.plot_y,color=groups[g],label=g,s=cfg['point_size'],alpha=cfg['point_alpha'])
                    ellipse=ellipse_points(x) if cfg['ellipse'] else None
                    if ellipse is not None: ax.plot(ellipse[:,0],ellipse[:,1],color=groups[g],lw=1,alpha=.65)
                if cfg['sample_labels']:
                    for row in data.itertuples(): ax.annotate(row.sample,(row.plot_x,row.plot_y),fontsize=cfg['label_font'],xytext=(4,4),textcoords='offset points')
            elif category in ('heatmap','matrix','association'):
                matrix,labels,annotations,lo,hi=heatmap_view(data,kind,root,cfg)
                if matrix.empty: empty_plot(ax,'No eligible genes / sets at this threshold')
                else:
                    cmap,_=palette(cfg);im=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=lo,vmax=hi);fig.colorbar(im,ax=ax,shrink=.75)
                    ax.set_xticks(range(len(matrix.columns)),matrix.columns,rotation=50,ha='right')
                    ax.set_yticks(range(len(matrix)),labels,fontsize=cfg['label_font'])
                    if cfg['show_values'] and matrix.size<=600:
                        for i in range(len(matrix)):
                            for j in range(len(matrix.columns)): ax.text(j,i,annotations[i][j],ha='center',va='center',fontsize=min(cfg['label_font'],9))
                    sample_group=dict(zip(m['sample'],m.group))
                    for j,name in enumerate(matrix.columns):
                        if name in sample_group: ax.plot(j,-.65,marker='s',color=groups[sample_group[name]],markersize=7,clip_on=False)
            elif category in ('dot','barh'):
                if data.empty: empty_plot(ax,'No gene sets at the selected threshold')
                else:
                    y=np.arange(len(data));labels=wrapped(data.description,cfg)
                    if category=='dot':
                        color=-np.log10(data.padj.clip(lower=1e-300))
                        if cfg['enrichment_style']=='bar': ax.barh(y,data.plot_x,color=cfg['down'])
                        else:
                            if cfg['enrichment_style']=='lollipop': ax.hlines(y,0,data.plot_x,color=cfg['neutral'],lw=2)
                            im=ax.scatter(data.plot_x,y,s=cfg['point_size']*(1+np.sqrt(data.overlap)),c=color,cmap='YlOrRd',alpha=cfg['point_alpha']);fig.colorbar(im,ax=ax,label='−log10 BH adjusted P',shrink=.75)
                    else:
                        color=[cfg['up'] if a>=0 else cfg['down'] for a in data.NES]
                        if cfg['enrichment_style']=='lollipop': ax.hlines(y,0,data.NES,color=color,lw=2);ax.scatter(data.NES,y,c=color,s=cfg['point_size'])
                        else: ax.barh(y,data.NES,color=color)
                        ax.axvline(0,color='#999999',lw=.7)
                    ax.set_yticks(y,labels,fontsize=cfg['label_font']);ax.invert_yaxis()
            elif category=='curve':
                fig.clear();ax,hit_ax,rank_ax=fig.subplots(3,1,sharex=True,gridspec_kw={'height_ratios':[3,.55,1]})
                ax.plot(data['rank'],data['running'],color=cfg['up'],lw=1.8);ax.axhline(0,color='#999999',ls='--',lw=.7)
                hit_ax.vlines(np.array(data['hits'])+1,0,1,color=cfg['down'],lw=.6);hit_ax.axis('off')
                rank_ax.fill_between(data['rank'],data['stat'],color=cfg['neutral']);rank_ax.set_xlabel('Gene rank');rank_ax.set_ylabel(ranking_label(root))
                ax.text(.02,.03,f"NES={data['result']['NES']:.3g}; padj={data['result']['padj']:.3g}",transform=ax.transAxes,fontsize=9)
            elif category in ('gene','trait'):
                x=pd.DataFrame(data['samples']);scale=data['scale'];rng=np.random.default_rng(42)
                if category=='gene':
                    if cfg['connect_pairs'] and 'subject' in x:
                        order={g:i for i,g in enumerate(x.group.unique())}
                        for _,subject in x.groupby('subject'):
                            points=subject.assign(plot_x=subject.group.map(order)).sort_values('plot_x')
                            ax.plot(points.plot_x,points[scale],color=cfg['neutral'],lw=.8,zorder=1)
                    for i,g in enumerate(x.group.unique()):
                        values=x.loc[x.group==g,scale].to_numpy()
                        if cfg['plot_style']=='violin' and len(values)>1 and np.std(values)>0:
                            body=ax.violinplot(values,positions=[i],widths=.7,showmedians=True)
                            for b in body['bodies']: b.set_facecolor(groups[g]);b.set_alpha(.25)
                        else: ax.boxplot([values],positions=[i],widths=.45,showfliers=False)
                        if cfg['show_points']: ax.scatter(i+rng.uniform(-.13,.13,len(values)),values,c=groups[g],s=cfg['point_size'],alpha=cfg['point_alpha'],zorder=3)
                    ax.set_xticks(range(x.group.nunique()),list(x.group.unique()))
                else:
                    for g in x.group.unique():
                        part=x[x.group==g];ax.scatter(part[data['trait']],part[scale],color=groups[g],label=g,s=cfg['point_size'],alpha=cfg['point_alpha'])
                    if cfg['sample_labels']:
                        for row in x.to_dict('records'): ax.annotate(row['sample'],(row[data['trait']],row[scale]),fontsize=8,xytext=(3,3),textcoords='offset points')
                    if data['adjusted']:
                        a=data['adjusted'][0];ax.text(.02,.98,f"Adjusted r={a['correlation']:.3g}; global q={a['padj']:.3g}",ha='left',va='top',transform=ax.transAxes,fontsize=9)
            elif category=='network':
                nodes={n['id']:n for n in data['nodes']}
                for e in data['edges']:
                    a,b=nodes[e['source']],nodes[e['target']];ax.plot([a['x'],b['x']],[a['y'],b['y']],c=cfg['neutral'],lw=.5+2*e['weight'],zorder=1)
                for node in nodes.values():
                    ax.scatter(node['x'],node['y'],s=cfg['point_size']*3,color=cfg['up'] if node['value']>=0 else cfg['down'],alpha=cfg['point_alpha'],zorder=2)
                    if cfg['sample_labels']: ax.annotate('\n'.join(textwrap.wrap(node['label'],22)),(node['x'],node['y']),fontsize=cfg['label_font'],xytext=(4,4),textcoords='offset points')
                if not nodes: empty_plot(ax,'No network edges / sets pass these filters')
                ax.margins(.2);ax.set_axis_off()
            elif category=='box':
                for i,s in enumerate(data.columns):
                    values=data[s].to_numpy();color=groups[m.set_index('sample').loc[s,'group']]
                    if cfg['plot_style']=='violin' and np.std(values)>0:
                        parts=ax.violinplot(values,positions=[i],showmedians=True)
                        for body in parts['bodies']: body.set_facecolor(color)
                    else: ax.boxplot([values],positions=[i],widths=.5,showfliers=False)
                ax.set_xticks(range(len(data.columns)),data.columns,rotation=45,ha='right')
            elif category in ('quality','base_content'):
                for label,values,color,dash in curve_series(data,cfg,kind): ax.plot(range(1,len(values)+1),values,label=label,color=color,ls='--' if dash=='dash' else '-')
                if category=='base_content': ax.set_ylim(0,1)
            elif category=='p_hist':
                for col,color in [('pvalue',cfg['down']),('padj',cfg['up'])]:
                    if col in data: ax.hist(data[col].dropna(),bins=np.linspace(0,1,21),histtype='step',lw=1.8,color=color,label=col)
            elif category=='pca_scree':
                x=data.head(int(cfg['top']));ax.bar(x.PC,x.variance*100,color=cfg['down']);ax.plot(x.PC,np.cumsum(x.variance)*100,'o-',color=cfg['up'],label='Cumulative');ax.set_ylim(0,101)
            elif category=='library_size': ax.bar(data['sample'],data.library_size,color=[groups[g] for g in data.group]);ax.tick_params(axis='x',labelrotation=45)
            elif category=='mapping':
                if 'percent_mapped' in data: ax.bar(data['sample'],data.percent_mapped,color=cfg['down']);ax.set_ylim(0,100)
                else: empty_plot(ax,'Mapping metadata unavailable')
            elif category in ('fastp_summary','deg','qc_rates'):
                series=[('raw_reads','Raw reads',cfg['neutral']),('clean_reads','Clean reads',cfg['down'])] if category=='fastp_summary' else [('up','Up',cfg['up']),('down','Down',cfg['down'])] if category=='deg' else [('q20_rate','Q20',cfg['down']),('q30_rate','Q30',cfg['up']),('gc_content','GC',PALETTE[0]),('retained_pct','Retained',PALETTE[3])]
                labels=data['name'] if category=='deg' else data['sample'];xx=np.arange(len(data));width=.8/len(series)
                for i,(col,label,color) in enumerate(series):
                    ax.bar(xx-.4+width*(i+.5),data[col]*(100 if category=='qc_rates' and col!='retained_pct' else 1),width=width,label=label,color=color)
                ax.set_xticks(xx,labels,rotation=40,ha='right')
                if category=='qc_rates': ax.set_ylim(0,100)
            elif category=='wgcna_soft':
                ax.plot(data.Power,-np.sign(data.slope)*data['SFT.R.sq'],'o-',color=cfg['up'],label='Signed fit');ax.axhline(.85,ls='--',color='#89958d');ax.set_ylim(-1,1)
                twin=ax.twinx();twin.plot(data.Power,data['mean.k.'],'s--',color=cfg['down'],label='Mean connectivity');twin.set_ylabel('Mean connectivity')
            elif category=='modules': ax.bar(data.module,data.genes,color=module_colors(data.module,cfg['neutral']));ax.tick_params(axis='x',labelrotation=45)
            elif category=='hubs':
                ax.barh(range(len(data)),data.kME,color=module_colors(data.module,cfg['neutral']));ax.set_yticks(range(len(data)),gene_labels(root,data.gene,cfg));ax.invert_yaxis()
                if data.empty: empty_plot(ax,'No assigned hub genes')
            elif category=='dispersion':
                x=data[(data.baseMean>0)&(data.dispersion>0)]
                if 'dispersion_gene' in x: ax.scatter(x.baseMean,x.dispersion_gene,s=cfg['point_size']/3,c=cfg['neutral'],alpha=.5,label='Gene-wise')
                ax.scatter(x.baseMean,x.dispersion,s=cfg['point_size']/3,c=cfg['down'],alpha=cfg['point_alpha'],label='Final')
                if 'dispersion_fit' in x:
                    z=x.sort_values('baseMean');ax.plot(z.baseMean,z.dispersion_fit,color=cfg['up'],lw=1.3,label='Fitted trend')
                ax.set_xscale('log');ax.set_yscale('log')
            apply_axes(ax,cfg,kind,category,axes)
            if not cfg['title'] and category in ('gene','trait'):
                ax.set_title(data['gene']+(' · '+data['trait'] if category=='trait' else ''),pad=14)
            if not cfg['title'] and category=='curve':
                ax.set_title(data['term'],pad=14)
            handles,labels=ax.get_legend_handles_labels()
            if handles and cfg['legend']: ax.legend(frameon=False,fontsize=min(10,cfg['font']))
            for fmt in ('svg','pdf','png'): fig.savefig(p/(kind+'.'+fmt),dpi=int(cfg['dpi']))
            if isinstance(data,pd.DataFrame): data.to_csv(p/(kind+'.plot-data.csv'))
            else: (p/(kind+'.plot-data.json')).write_text(json.dumps(data,ensure_ascii=False,allow_nan=False))
            (p/(kind+'.settings.json')).write_text(json.dumps(cfg,indent=2,ensure_ascii=False))
        finally: plt.close(fig)


def empty_plot(ax,message):
    ax.text(.5,.5,message,transform=ax.transAxes,ha='center',va='center',wrap=True);ax.set_axis_off()


def apply_axes(ax,cfg,kind,category,axes):
    if cfg['grid'] and category not in ('heatmap','matrix','association','network'):
        ax.grid(True,color='#e8ece9',lw=.6);ax.set_axisbelow(True)
    for axis in ('x','y'):
        low,high=cfg[axis+'_min'],cfg[axis+'_max']
        if low!='' or high!='':
            limits=getattr(ax,'get_'+axis+'lim')();lo=float(low) if low!='' else limits[0];hi=float(high) if high!='' else limits[1]
            if getattr(ax,'get_'+axis+'scale')()=='log' and min(lo,hi)<=0: raise ValueError('对数坐标范围必须大于 0')
            getattr(ax,'set_'+axis+'lim')(lo,hi)
    ax.set_xlabel(cfg['xlabel'] or axes[0]);ax.set_ylabel(cfg['ylabel'] or axes[1]);ax.set_title(cfg['title'] or TITLES[kind],pad=14)


def interactive_data(p,kind,args,options=None):
    from results import contrast_dir, records
    target=contrast_dir(p,args.get('contrast'));overrides={k:args[k] for k in ('term','gene','database','sample','annotated_genes','trait','module') if args.get(k)}
    overrides.update(options or {})
    cfg=get_settings(target,kind,overrides)
    if kind in ('gsea_curve','leading_edge','enrichment_network'):cfg['database']=gsea_database(target,cfg)
    category,data,axes=load_plot(target,kind,cfg)
    axes=tuple(label.replace('VST',expression_label(data_root(target))) for label in axes)
    root=data_root(target);m=read_frame(root/'samples.csv');groups=colors_for(m,root,cfg)
    layout=dict(title=dict(text=cfg['title'] or TITLES[kind],x=.5,font=dict(size=cfg['font']+4)),xaxis=dict(title=dict(text=axes[0]),showgrid=cfg['grid'],gridcolor='#e8ece9',zeroline=False,automargin=True),yaxis=dict(title=dict(text=axes[1]),showgrid=cfg['grid'],gridcolor='#e8ece9',zeroline=False,automargin=True),margin=dict(l=75,r=45,t=85,b=85),font=dict(size=cfg['font']),paper_bgcolor='#ffffff',plot_bgcolor='#ffffff',showlegend=cfg['legend'],legend=dict(orientation='h',y=-.18,x=.5,xanchor='center'),hovermode='closest')
    traces=[];marker=dict(size=np.sqrt(cfg['point_size'])*1.5,opacity=cfg['point_alpha'])
    if category=='scatter':
        # Keep gene ID as customdata for point click compatibility; extra statistics
        # live in hovertext so both label and point remain independently editable.
        names=(data.symbol.fillna(data.gene) if 'symbol' in data else data.gene)
        hover=[]
        for (_,row),name in zip(data.iterrows(),names):
            hover.append(html.escape(str(name))+'<br>'+html.escape(str(row.gene))+f'<br>log2FC={row.log2FoldChange:.4g}'+''.join(f'<br>{c}={row[c]:.4g}' if c in row and pd.notna(row[c]) else '<br>'+c+'=NA' for c in ('pvalue','padj')))
        traces=[dict(type='scattergl',mode='markers',x=data.plot_x.tolist(),y=data.plot_y.tolist(),text=names.tolist(),hovertext=hover,customdata=data.gene.tolist(),showlegend=False,marker=dict(color=data.color.tolist(),size=np.where(data.manual_label,12,marker['size']).tolist(),opacity=cfg['point_alpha'],line=dict(width=.6,color='#637169')),hovertemplate='%{hovertext}<extra></extra>')]
        for name,color in [('Up',cfg['up']),('Down',cfg['down']),('Other',cfg['neutral'])]:
            traces.append(dict(type='scatter',mode='markers',x=[None],y=[None],name=f'{name} ({int((data.direction_plot==name).sum())})',marker=dict(color=color,size=8),hoverinfo='skip'))
        layout['annotations']=[]
        for i,row in enumerate(selected_labels(data,cfg).itertuples()):
            pos=cfg['label_positions'].get(row.gene,dict(ax=30 if i%2==0 else -30,ay=-30-22*(i%3)))
            label=str(row.symbol) if hasattr(row,'symbol') and pd.notna(row.symbol) else row.gene
            # Plotly annotation coordinates on log axes use log10(data coordinates).
            xx=float(np.log10(row.plot_x)) if kind=='ma' else float(row.plot_x)
            layout['annotations'].append(dict(name=row.gene,x=xx,y=float(row.plot_y),text='<i>'+html.escape(label)+'</i>',showarrow=True,arrowhead=0,arrowcolor='#68756d',ax=pos.get('ax',25),ay=pos.get('ay',-25),bgcolor='white',bordercolor='#919b94',borderwidth=1,borderpad=3,font=dict(size=cfg['label_font'],color='#26382e')))
        if kind=='volcano':
            if cfg['threshold_lines']:
                layout['shapes']=[dict(type='line',x0=val,x1=val,y0=0,y1=1,yref='paper',line=dict(dash='dash',color='#77857c',width=1)) for val in sorted({-cfg['lfc'],cfg['lfc']})]
                layout['shapes'].append(dict(type='line',x0=0,x1=1,xref='paper',y0=-np.log10(cfg['alpha']),y1=-np.log10(cfg['alpha']),line=dict(dash='dash',color='#77857c',width=1)))
            if cfg['symmetric_x'] and len(data):
                extent=max(1,float(data.plot_x.abs().max())*1.12);layout['xaxis']['range']=[-extent,extent]
            layout['yaxis']['rangemode']='tozero'
        else: layout['xaxis']['type']='log'
    elif category=='pca':
        for g in data.group.unique():
            x=data[data.group==g];traces.append(dict(type='scatter',mode='markers+text' if cfg['sample_labels'] else 'markers',name=g,x=x.plot_x.tolist(),y=x.plot_y.tolist(),text=x['sample'].tolist(),textposition='top center',textfont=dict(size=cfg['label_font']),marker=dict(color=groups[g],**marker),hovertemplate='%{text}<br>%{x:.3g}, %{y:.3g}<extra>%{fullData.name}</extra>'))
            ellipse=ellipse_points(x) if cfg['ellipse'] else None
            if ellipse is not None: traces.append(dict(type='scatter',mode='lines',x=ellipse[:,0].tolist(),y=ellipse[:,1].tolist(),line=dict(color=groups[g],width=1),showlegend=False,hoverinfo='skip'))
    elif category in ('heatmap','matrix','association'):
        matrix,labels,annotations,lo,hi=heatmap_view(data,kind,root,cfg);_,colorscale=palette(cfg)
        trace=dict(type='heatmap',z=matrix.to_numpy().tolist(),x=list(matrix.columns),y=labels,colorscale=colorscale,hovertemplate='%{y}<br>%{x}: %{z:.4g}<extra></extra>',colorbar=dict(thickness=12,len=.7),showscale=True)
        if lo is not None: trace['zmin']=lo
        if hi is not None: trace['zmax']=hi
        if cfg['show_values'] and matrix.size<=600: trace.update(text=[[v.replace('\n','<br>') for v in row] for row in annotations],texttemplate='%{text}',textfont=dict(size=min(cfg['label_font'],11)))
        if kind in ('heatmap','leading_edge','activity','tf_activity'):
            trace['customdata']=[[g]*len(matrix.columns) for g in matrix.index]
            trace['hovertemplate']='%{y}<br>%{customdata}<br>%{x}: %{z:.4g}<extra></extra>'
        if category=='association':
            trace['text']=[[v.replace('\n','<br>') for v in row] for row in annotations];trace['hovertemplate']='%{y}<br>%{x}<br>%{text}<extra></extra>'
        traces=[trace];layout['yaxis']['autorange']='reversed';layout['yaxis']['tickfont']=dict(size=cfg['label_font']);layout['margin']['l']=130
        if matrix.empty: layout['annotations']=[dict(x=.5,y=.5,xref='paper',yref='paper',text='没有符合当前条件的基因或基因集',showarrow=False)]
        # A separate categorical annotation strip follows the actual plotted order.
        sample_group=dict(zip(m['sample'],m.group))
        if len(matrix.columns) and all(s in sample_group for s in matrix.columns):
            traces.append(dict(type='scatter',mode='markers',x=list(matrix.columns),y=[.5]*len(matrix.columns),yaxis='y2',marker=dict(symbol='square',size=10,color=[groups[sample_group[s]] for s in matrix.columns]),text=[sample_group[s] for s in matrix.columns],showlegend=False,hovertemplate='%{x}: %{text}<extra></extra>'))
            layout['yaxis']['domain']=[0,.94]
            layout['yaxis2']=dict(domain=[.965,1],range=[0,1],visible=False,fixedrange=True,anchor='x')
    elif category in ('dot','barh'):
        labels=wrapped(data.description,cfg,True);ids=data.term.tolist()
        if category=='dot':
            hover=[f'{html.escape(str(row.term))}<br>{html.escape(str(row.description))}<br>k={row.overlap}, K={row.size}, n={row.selected}, N={row.background}<br>padj={row.padj:.4g}' for row in data.itertuples()]
            if cfg['enrichment_style']=='bar': traces=[dict(type='bar',orientation='h',x=data.plot_x.tolist(),y=labels,customdata=ids,marker=dict(color=cfg['down']),hovertext=hover,hovertemplate='%{hovertext}<extra></extra>')]
            else:
                if cfg['enrichment_style']=='lollipop':
                    xx=[];yy=[]
                    for val,label in zip(data.plot_x,labels): xx.extend([0,val,None]);yy.extend([label,label,None])
                    traces.append(dict(type='scatter',mode='lines',x=xx,y=yy,line=dict(color=cfg['neutral'],width=2),hoverinfo='skip',showlegend=False))
                traces.append(dict(type='scatter',mode='markers',x=data.plot_x.tolist(),y=labels,customdata=ids,marker=dict(size=(np.sqrt(cfg['point_size'])*(1+np.sqrt(data.overlap)/2)).tolist(),color=(-np.log10(data.padj.clip(lower=1e-300))).tolist(),colorscale='YlOrRd',showscale=True,opacity=cfg['point_alpha'],colorbar=dict(title=dict(text='−log10 q'),thickness=12,len=.7)),hovertext=hover,hovertemplate='%{hovertext}<extra></extra>',showlegend=False))
        else:
            color=[cfg['up'] if val>=0 else cfg['down'] for val in data.NES]
            traces=[dict(type='bar',orientation='h',x=data.NES.tolist(),y=labels,customdata=ids,marker=dict(color=color),text=[f'q={v:.3g}' for v in data.padj],hovertemplate='%{y}<br>effect=%{x:.4g}<br>%{text}<extra></extra>',showlegend=False)]
            if cfg['enrichment_style']=='lollipop':
                traces[0].update(type='scatter',mode='markers',marker=dict(color=color,size=marker['size']));traces[0].pop('orientation',None)
                for val,label,c in zip(data.NES,labels,color): traces.append(dict(type='scatter',mode='lines',x=[0,val],y=[label,label],line=dict(color=c,width=2),showlegend=False,hoverinfo='skip'))
        layout['yaxis']['autorange']='reversed';layout['yaxis']['tickfont']=dict(size=cfg['label_font']);layout['margin']['l']=min(270,int(cfg['wrap'])*6);layout['showlegend']=False
    elif category=='curve':
        traces=[dict(type='scatter',mode='lines',name=data['term'],x=data['rank'],y=data['running'],line=dict(color=cfg['up'],width=2)),dict(type='scatter',mode='markers',name='Gene-set hits',x=(np.array(data['hits'])+1).tolist(),y=[.5]*len(data['hits']),text=[data['gene'][i] for i in data['hits']],marker=dict(symbol='line-ns-open',size=15,color=cfg['down']),yaxis='y2',hovertemplate='%{text}<extra></extra>'),dict(type='scatter',mode='lines',name=ranking_label(root),x=data['rank'],y=data['stat'],fill='tozeroy',line=dict(color=cfg['neutral'],width=1),yaxis='y3')]
        layout['yaxis'].update(domain=[.42,1]);layout['yaxis2']=dict(domain=[.3,.37],visible=False,range=[0,1]);layout['yaxis3']=dict(domain=[0,.25],title=dict(text=ranking_label(root)),zeroline=True);layout['xaxis']['anchor']='y3'
        layout['title']['text']=cfg['title'] or html.escape(data['term'])+f"<br><sup>NES {data['result']['NES']:.3g} · BH padj {data['result']['padj']:.3g}</sup>";layout['showlegend']=False
    elif category in ('gene','trait'):
        frame=pd.DataFrame(data['samples']);scale=data['scale']
        if category=='gene':
            for g in frame.group.unique():
                x=frame[frame.group==g];trace=dict(type=cfg['plot_style'],name=g,y=x[scale].tolist(),text=x['sample'].tolist(),marker=dict(color=groups[g],**marker),hovertemplate='%{text}<br>%{y:.3g}<extra></extra>')
                if cfg['plot_style']=='box': trace.update(boxpoints='all' if cfg['show_points'] else False,jitter=.3)
                else: trace.update(points='all' if cfg['show_points'] else False,box=dict(visible=True),meanline=dict(visible=True),jitter=.2)
                traces.append(trace)
            if cfg['connect_pairs'] and 'subject' in frame:
                for subject,x in frame.groupby('subject'):
                    if len(x)>1: traces.append(dict(type='scatter',mode='lines',x=x.group.tolist(),y=x[scale].tolist(),line=dict(color=cfg['neutral'],width=1),name=str(subject),showlegend=False,hoverinfo='skip'))
        else:
            for g in frame.group.unique():
                x=frame[frame.group==g];traces.append(dict(type='scatter',mode='markers+text' if cfg['sample_labels'] else 'markers',name=g,x=x[data['trait']].tolist(),y=x[scale].tolist(),text=x['sample'].tolist(),textposition='top center',marker=dict(color=groups[g],**marker)))
            if data['adjusted']:
                a=data['adjusted'][0];layout['annotations']=[dict(x=.5,y=1.04,xref='paper',yref='paper',text=f"Adjusted {expression_label(root)} r={a['correlation']:.3g}; global q={a['padj']:.3g}",showarrow=False)]
        layout['title']['text']=cfg['title'] or html.escape(data['gene'])
    elif category=='network':
        nodes={node['id']:node for node in data['nodes']}
        xx=[];yy=[]
        for edge in data['edges']:
            a,b=nodes[edge['source']],nodes[edge['target']];xx.extend([a['x'],b['x'],None]);yy.extend([a['y'],b['y'],None])
        traces=[dict(type='scatter',mode='lines',x=xx,y=yy,line=dict(color=cfg['neutral'],width=1),hoverinfo='skip',showlegend=False)]
        traces.append(dict(type='scatter',mode='markers+text' if cfg['sample_labels'] else 'markers',x=[n['x'] for n in nodes.values()],y=[n['y'] for n in nodes.values()],customdata=list(nodes),text=[html.escape(n['label']) for n in nodes.values()],textfont=dict(size=cfg['label_font']),textposition='top center',hovertext=[html.escape(n['id'])+f"<br>value={n['value']:.3g}" for n in nodes.values()],hovertemplate='%{hovertext}<extra></extra>',marker=dict(size=marker['size']*2,color=[cfg['up'] if n['value']>=0 else cfg['down'] for n in nodes.values()],opacity=cfg['point_alpha']),showlegend=False))
        layout['xaxis'].update(visible=False);layout['yaxis'].update(visible=False);layout['margin']=dict(l=35,r=35,t=70,b=35)
        if not nodes: layout['annotations']=[dict(x=.5,y=.5,xref='paper',yref='paper',text='没有符合条件的节点/边；可调整筛选与阈值',showarrow=False)]
    elif category=='box':
        for s in data.columns:
            trace=dict(type=cfg['plot_style'],name=s,y=data[s].tolist(),marker=dict(color=groups[m.set_index('sample').loc[s,'group']]))
            trace.update(boxpoints=False) if cfg['plot_style']=='box' else trace.update(points=False,box=dict(visible=True))
            traces.append(trace)
        layout['showlegend']=False
    elif category in ('quality','base_content'):
        for label,values,color,dash in curve_series(data,cfg,kind): traces.append(dict(type='scatter',mode='lines',name=label,x=list(range(1,len(values)+1)),y=values,line=dict(color=color,dash=dash)))
        if category=='base_content': layout['yaxis']['range']=[0,1]
    elif category=='p_hist':
        for col,color in [('pvalue',cfg['down']),('padj',cfg['up'])]:
            if col in data: traces.append(dict(type='histogram',name=col,x=data[col].dropna().tolist(),xbins=dict(start=0,end=1,size=.05),opacity=.55,marker=dict(color=color)))
        layout['barmode']='overlay'
    elif category=='pca_scree':
        x=data.head(int(cfg['top']));traces=[dict(type='bar',name='Per PC',x=x.PC.tolist(),y=(100*x.variance).tolist(),marker=dict(color=cfg['down'])),dict(type='scatter',mode='lines+markers',name='Cumulative',x=x.PC.tolist(),y=(100*np.cumsum(x.variance)).tolist(),line=dict(color=cfg['up']))];layout['yaxis']['range']=[0,101]
    elif category=='library_size': traces=[dict(type='bar',x=data['sample'].tolist(),y=data.library_size.tolist(),marker=dict(color=[groups[g] for g in data.group]),showlegend=False)]
    elif category=='mapping':
        traces=[dict(type='bar',x=data['sample'].tolist(),y=data.percent_mapped.tolist() if 'percent_mapped' in data else [],marker=dict(color=cfg['down']),showlegend=False)];layout['yaxis']['range']=[0,100]
        if 'percent_mapped' not in data: layout['annotations']=[dict(x=.5,y=.5,xref='paper',yref='paper',text='导入文件未提供映射率元数据，不能推算',showarrow=False)]
    elif category in ('fastp_summary','deg','qc_rates'):
        series=[('raw_reads','Raw reads',cfg['neutral']),('clean_reads','Clean reads',cfg['down'])] if category=='fastp_summary' else [('up','Up',cfg['up']),('down','Down',cfg['down'])] if category=='deg' else [('q20_rate','Q20',cfg['down']),('q30_rate','Q30',cfg['up']),('gc_content','GC',PALETTE[0]),('retained_pct','Retained',PALETTE[3])]
        for col,label,color in series: traces.append(dict(type='bar',name=label,x=(data['name'] if category=='deg' else data['sample']).tolist(),y=(data[col]*(100 if category=='qc_rates' and col!='retained_pct' else 1)).tolist(),marker=dict(color=color)))
        layout['barmode']='group'
        if category=='qc_rates': layout['yaxis']['range']=[0,100]
    elif category=='wgcna_soft':
        traces=[dict(type='scatter',mode='lines+markers',x=data.Power.tolist(),y=(-np.sign(data.slope)*data['SFT.R.sq']).tolist(),name='Signed scale-free fit',line=dict(color=cfg['up'])),dict(type='scatter',mode='lines+markers',x=data.Power.tolist(),y=data['mean.k.'].tolist(),name='Mean connectivity',line=dict(color=cfg['down'],dash='dash'),yaxis='y2')]
        layout['yaxis']['range']=[-1,1];layout['yaxis2']=dict(overlaying='y',side='right',title=dict(text='Mean connectivity'));layout['margin']['r']=85
        layout['shapes']=[dict(type='line',xref='paper',x0=0,x1=1,y0=.85,y1=.85,line=dict(dash='dash',color='#89958d'))]
    elif category=='modules': traces=[dict(type='bar',x=data.module.tolist(),y=data.genes.tolist(),marker=dict(color=module_colors(data.module,cfg['neutral'])),showlegend=False)]
    elif category=='hubs':
        traces=[dict(type='bar',orientation='h',x=data.kME.tolist(),y=gene_labels(root,data.gene,cfg),customdata=data.gene.tolist(),marker=dict(color=module_colors(data.module,cfg['neutral'])),showlegend=False)];layout['yaxis']['autorange']='reversed'
    elif category=='dispersion':
        x=data[(data.baseMean>0)&(data.dispersion>0)]
        for col,label,color in [('dispersion_gene','Gene-wise',cfg['neutral']),('dispersion','Final',cfg['down'])]:
            if col in x: traces.append(dict(type='scattergl',mode='markers',x=x.baseMean.tolist(),y=x[col].tolist(),name=label,marker=dict(color=color,size=marker['size']/2,opacity=cfg['point_alpha'])))
        if 'dispersion_fit' in x:
            x=x.sort_values('baseMean');traces.append(dict(type='scatter',mode='lines',x=x.baseMean.tolist(),y=x.dispersion_fit.tolist(),name='Fitted trend',line=dict(color=cfg['up'])))
        layout['xaxis']['type']='log';layout['yaxis']['type']='log'
    for axis in ('x','y'):
        if cfg[axis+'label']: layout[axis+'axis']['title']=dict(text=cfg[axis+'label'])
        low,high=cfg[axis+'_min'],cfg[axis+'_max']
        if low!='' or high!='':
            limits=[float(low) if low!='' else None,float(high) if high!='' else None]
            if layout[axis+'axis'].get('type')=='log':
                if any(val is not None and val<=0 for val in limits): raise ValueError('对数坐标范围必须大于 0')
                limits=[float(np.log10(val)) if val is not None else None for val in limits]
            layout[axis+'axis']['range']=limits
    config=dict(responsive=True,displaylogo=False,editable=True,edits=dict(annotationPosition=True,annotationText=False,axisTitleText=False,titleText=False,legendPosition=True),toImageButtonOptions=dict(format='svg',filename=kind,width=int(cfg['width']*96),height=int(cfg['height']*96)))
    return json_safe(dict(data=traces,layout=layout,config=config,settings=cfg,category=category,metadata=data if isinstance(data,dict) else None,note=NOTES.get(kind,'悬停查看数据，拖动缩放。保存后同步更新 SVG、PDF、PNG、绘图数据和报告。').replace('VST',expression_label(root))))


def json_safe(value):
    if isinstance(value,dict): return {key:json_safe(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)): return [json_safe(item) for item in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value): return None
    return value


def expression_label(root):
    cfg=json.loads((root/'config.json').read_text()) if (root/'config.json').exists() else {}
    return {'vst':'VST','log2_normalized':'log2(normalized count + 1)','log2_tpm':'log2(TPM + 1)'}.get(cfg.get('expression_unit','vst'),'VST')


def ranking_label(root):
    cfg=json.loads((root/'config.json').read_text()) if (root/'config.json').exists() else {}
    return 'Imported stat' if cfg.get('analysis_stage')=='enrichment' else 'Wald stat'
