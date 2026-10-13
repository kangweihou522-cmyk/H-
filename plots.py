"""Scientific exports and offline Plotly data share the same settings and input tables."""
from pathlib import Path
import html
from adjustText import adjust_text
import json, re, os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, is_color_like
from scipy.cluster.hierarchy import linkage, leaves_list
import numpy as np
import pandas as pd

TITLES={'base_content':'Base content','quality':'Read quality','fastp_summary':'Sequencing QC','mapping':'Salmon mapping rate','library_size':'Library size','distribution':'Expression distribution (VST)','pca':'Sample PCA','correlation':'Sample Pearson correlation (VST)','sample_distance':'Sample distance (VST)','dispersion':'DESeq2 dispersion','volcano':'Differential expression','ma':'MA plot','heatmap':'Differential gene heatmap','gene_expression':'Gene expression','go':'GO enrichment','gsea':'GO GSEA','pathway':'GMT / KEGG enrichment','pathway_gsea':'GMT / KEGG GSEA','gsea_curve':'GSEA enrichment curve','tf_activity':'Regulatory activity (decoupleR ULM)','tf_de':'Differential regulatory activity','activity':'Pathway activity (ssGSEA)','activity_de':'Differential pathway activity','wgcna_soft':'WGCNA soft-threshold diagnostics','wgcna_traits':'Module / trait association','wgcna_modules':'WGCNA modules','deg_summary':'Differential gene overview'}
LABELS={'base_content':'A/T/G/C 含量曲线','quality':'碱基质量曲线','fastp_summary':'原始/过滤 reads','mapping':'Salmon 映射率','library_size':'文库大小','distribution':'表达分布','pca':'PCA 样本诊断','correlation':'样本相关性热图','sample_distance':'样本距离热图','dispersion':'离散度诊断','volcano':'火山图','ma':'MA 图','heatmap':'差异基因聚类热图','gene_expression':'基因表达','go':'GO 富集点图','gsea':'GO GSEA','pathway':'KEGG / GMT 富集点图','pathway_gsea':'KEGG / GMT GSEA','gsea_curve':'GSEA 富集曲线','tf_activity':'调控因子活性热图','tf_de':'差异调控活性','activity':'通路活性热图','activity_de':'差异通路活性','wgcna_soft':'WGCNA 软阈值诊断','wgcna_traits':'模块–表型关联','wgcna_modules':'共表达模块','deg_summary':'各比较差异基因统计'}
PALETTE=['#126657','#b64c4c','#317ba0','#bd9a4b','#8668ab','#c96d91','#73913b','#77858c']

def module_colors(modules,neutral):
    # WGCNA uses R names such as lightcyan1 which are not CSS/Matplotlib colors.
    return [neutral if name=='grey' else name if is_color_like(name) else PALETTE[i%len(PALETTE)] for i,name in enumerate(modules)]

def data_root(p):
    return p.parent.parent if p.parent.name=='contrasts' else p

def get_settings(p,kind,opts=None):
    root=data_root(p);analysis=json.loads((root/'config.json').read_text()) if (root/'config.json').exists() else {}
    f=p/(kind+'.settings.json');saved=json.loads(f.read_text()) if f.exists() else {}
    cfg=dict(title='',width=9,height=6,font=11,dpi=300,up='#b64c4c',down='#317ba0',neutral='#bac5c0',alpha=analysis.get('alpha',.05),lfc=analysis.get('lfc',1),top=20,xlabel='',ylabel='',direction='all',significant_only=True,cluster=True,labels=5,genes='',gene='',term='',database='go',pcx=1,pcy=2,sample='',expression_scale='vst',gene_labels='symbol',annotated_genes='',label_positions={},point_size=36,point_alpha=.85,label_font=10,grid=True,x_min='',x_max='',y_min='',y_max='',p_column=analysis.get('significance','padj'),group_colors=analysis.get('group_colors',{}))
    if (p/'active-profile.json').exists():
        active=json.loads((p/'active-profile.json').read_text());cfg.update(p_column=active['column'],alpha=active['alpha'],lfc=active['lfc'])
    cfg.update(saved);cfg.update(opts or {})
    for key,lo,hi in [('width',3,20),('height',3,20),('font',6,30),('dpi',72,600),('top',3,100),('alpha',1e-6,.999),('lfc',0,10),('labels',0,30),('pcx',1,5),('pcy',1,5),('point_size',5,150),('point_alpha',.15,1),('label_font',6,22)]:
        cfg[key]=float(cfg[key])
        if not lo<=cfg[key]<=hi: raise ValueError('图表参数超出范围：'+key)
    for key in ('up','down','neutral'):
        if not re.fullmatch('#[0-9a-fA-F]{6}',str(cfg[key])):raise ValueError('颜色需要 #RRGGBB')
    for key in ('title','xlabel','ylabel','gene','term','sample','genes','annotated_genes','x_min','x_max','y_min','y_max'):
        if not isinstance(cfg[key],str) or len(cfg[key])>(5000 if key in ('genes','annotated_genes') else 180):raise ValueError('图表文本无效/过长')
    if not isinstance(cfg['label_positions'],dict):raise ValueError('标签位置无效')
    for gene,position in cfg['label_positions'].items():
        if not isinstance(position,dict) or any(not isinstance(position.get(k,0),(int,float)) or abs(position.get(k,0))>600 for k in ('ax','ay')):raise ValueError('标签位置超出范围')
    for axis in ('x','y'):
        low,high=cfg[axis+'_min'],cfg[axis+'_max']
        if low!='' and not np.isfinite(float(low)):raise ValueError('坐标范围无效')
        if high!='' and not np.isfinite(float(high)):raise ValueError('坐标范围无效')
        if low!='' and high!='' and float(low)>=float(high):raise ValueError('坐标最小值必须小于最大值')
    if not isinstance(cfg['group_colors'],dict) or any(not re.fullmatch('#[0-9a-fA-F]{6}',str(x)) for x in cfg['group_colors'].values()):raise ValueError('分组颜色格式无效')
    if cfg['expression_scale'] not in ('vst','normalized','tpm') or cfg['gene_labels'] not in ('symbol','id'):raise ValueError('表达/基因标签选择无效')
    if cfg['p_column'] not in ('padj','pvalue'):raise ValueError('显著性列无效')
    if cfg['direction'] not in ('all','up','down') or cfg['database'] not in ('go','pathway'):raise ValueError('图表筛选无效')
    return cfg

def chart_catalog(p):
    charts=[]
    def add(kind,stage,contrast=None):
        charts.append(dict(kind=kind,title=LABELS[kind],stage=stage,contrast=contrast))
    for kind,file in [('fastp_summary','fastp-qc.csv'),('quality','fastp-qc.csv'),('base_content','fastp-qc.csv'),('mapping','salmon-qc.csv'),('library_size','sample-qc.csv'),('distribution','vst.csv'),('pca','vst.csv'),('correlation','sample-correlation.csv'),('sample_distance','sample-distance.csv'),('dispersion','dispersion.csv'),('deg_summary','analysis-summary.json'),('activity','pathway-activity.csv'),('tf_activity','tf-activity.csv'),('wgcna_soft','wgcna-soft-threshold.csv'),('wgcna_traits','wgcna-module-trait.csv'),('wgcna_modules','wgcna-modules.csv')]:
        if (p/file).exists():add(kind,'qc' if kind in ('fastp_summary','quality','base_content','mapping') else 'mining' if kind in ('activity','tf_activity','wgcna_soft','wgcna_traits','wgcna_modules') else 'expression')
    targets=[(q.name,q) for q in sorted((p/'contrasts').iterdir()) if q.is_dir()] if (p/'contrasts').is_dir() else [(None,p)]
    for name,q in targets:
        for kind,file in [('volcano','differential.csv'),('ma','differential.csv'),('heatmap','differential.csv'),('gene_expression','differential.csv'),('go','go.csv'),('gsea','gsea.csv'),('pathway','pathway.csv'),('pathway_gsea','pathway_gsea.csv'),('activity_de','pathway-activity-differential.csv'),('tf_de','tf-activity-differential.csv')]:
            if not (q/file).exists():continue
            if kind in ('go','gsea','pathway','pathway_gsea','activity_de','tf_de') and pd.read_csv(q/file).empty:continue
            add(kind,'differential' if kind in ('volcano','ma','heatmap','gene_expression') else 'mining' if kind in ('activity_de','tf_de') else 'enrichment',name)
        if (p/'gene-sets.json').exists() and (q/'gsea.csv').exists() and not pd.read_csv(q/'gsea.csv').dropna(subset=['padj','NES']).empty:add('gsea_curve','enrichment',name)
        elif (p/'gene-sets.json').exists() and (q/'pathway_gsea.csv').exists() and not pd.read_csv(q/'pathway_gsea.csv').dropna(subset=['padj','NES']).empty:add('gsea_curve','enrichment',name)
    return charts

def colors_for(m,root,cfg):
    analysis=json.loads((root/'config.json').read_text()) if (root/'config.json').exists() else {}
    custom=cfg.get('group_colors') or analysis.get('group_colors',{})
    groups=list(m.group.unique());return {g:custom.get(g,PALETTE[i%len(PALETTE)]) for i,g in enumerate(groups)}

def load_plot(p,kind,cfg):
    root=data_root(p);m=pd.read_csv(root/'samples.csv');v=pd.read_csv(root/'vst.csv',index_col=0)
    r=pd.read_csv(p/'differential.csv') if (p/'differential.csv').exists() else None
    if kind in ('volcano','ma'):
        x=r.dropna(subset=['log2FoldChange',cfg['p_column']]).copy();sign=x[cfg['p_column']]<cfg['alpha']
        x['color']=np.where(sign&(x.log2FoldChange>cfg['lfc']),cfg['up'],np.where(sign&(x.log2FoldChange<-cfg['lfc']),cfg['down'],cfg['neutral']))
        x['plot_x']=x.log2FoldChange if kind=='volcano' else x.baseMean
        manual={g.strip().casefold() for g in re.split(r'[,;\n]+',cfg['annotated_genes']) if g.strip()}
        x['manual_label']=x.gene.astype(str).str.casefold().isin(manual) | (x.symbol.fillna('').astype(str).str.casefold().isin(manual) if 'symbol' in x else False)
        x['plot_y']=-np.log10(x[cfg['p_column']].clip(lower=1e-300)) if kind=='volcano' else x.log2FoldChange
        return 'scatter',x,('log2 fold change','−log10 BH adjusted P' if cfg['p_column']=='padj' else '−log10 unadjusted P') if kind=='volcano' else ('Mean normalized expression','log2 fold change')
    if kind=='pca':
        if (root/'pca-scores.csv').exists():x=pd.read_csv(root/'pca-scores.csv');variance=pd.read_csv(root/'pca-variance.csv').variance.to_numpy()
        else:
            a=v[m['sample']].T.to_numpy(dtype=float,copy=True);a-=a.mean(axis=0);u,s,_=np.linalg.svd(a,full_matrices=False);variance=s*s/(s*s).sum();x=pd.DataFrame({'sample':m['sample'],'group':m.group,**{f'PC{i+1}':u[:,i]*s[i] for i in range(min(5,u.shape[1]))}})
        a,b=int(cfg['pcx']),int(cfg['pcy'])
        if f'PC{a}' not in x or f'PC{b}' not in x or a==b:raise ValueError('请选择存在且不同的主成分')
        x['plot_x']=x[f'PC{a}'];x['plot_y']=x[f'PC{b}']
        return 'pca',x,(f'PC{a} ({variance[a-1]:.1%})',f'PC{b} ({variance[b-1]:.1%})')
    if kind=='heatmap':
        manual=[x.strip() for x in re.split(r'[,;\n]+',cfg['genes']) if x.strip()]
        if manual:
            aliases=dict(zip(r.symbol.dropna(),r[r.symbol.notna()].gene)) if 'symbol' in r else {}
            genes=[g if g in v.index else aliases.get(g) for g in manual];genes=[g for g in genes if g in v.index]
            if not genes:raise ValueError('指定的基因未在保留基因中找到')
        else:genes=list(r[r.significant if 'significant' in r else ((r[cfg['p_column']]<cfg['alpha'])&(r.log2FoldChange.abs()>cfg['lfc']))].dropna(subset=[cfg['p_column']]).nsmallest(int(cfg['top']),cfg['p_column']).gene)
        x=v.loc[v.index.intersection(genes),m['sample']].copy()
        x=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1).replace(0,1),axis=0)
        if cfg['cluster'] and x.shape[0]>1 and x.shape[1]>1:
            x=x.iloc[leaves_list(linkage(x.to_numpy(),method='average')),leaves_list(linkage(x.T.to_numpy(),method='average'))]
        return 'heatmap',x,('Samples','Row Z-score (VST)')
    if kind in ('correlation','sample_distance'):
        x=pd.read_csv(root/('sample-correlation.csv' if kind=='correlation' else 'sample-distance.csv'),index_col=0)
        if cfg['cluster'] and len(x)>1:
            order=leaves_list(linkage(v[m['sample']].T.to_numpy(),method='average'));x=x.iloc[order,order]
        return 'matrix',x,('Samples','Samples')
    if kind in ('activity','tf_activity'):
        x=pd.read_csv(root/('tf-activity.csv' if kind=='tf_activity' else 'pathway-activity.csv'),index_col=0);x=x.loc[x.var(axis=1).nlargest(int(cfg['top'])).index,m['sample']]
        return 'heatmap',x,('Samples','ULM score' if kind=='tf_activity' else 'ssGSEA score')
    if kind=='wgcna_soft':return 'soft',pd.read_csv(root/'wgcna-soft-threshold.csv'),('Power','Signed scale-free fit R²')
    if kind=='wgcna_traits':
        x=pd.read_csv(root/'wgcna-module-trait.csv').pivot(index='module',columns='trait',values='correlation')
        return 'matrix',x,('Traits','Modules')
    if kind in ('go','pathway','gsea','pathway_gsea','activity_de','tf_de'):
        x=pd.read_csv(p/(('tf-activity-differential.csv' if kind=='tf_de' else 'pathway-activity-differential.csv') if kind in ('activity_de','tf_de') else kind+'.csv'))
        if kind in ('go','pathway') and 'direction' in x:x=x[x.direction==cfg['direction']]
        if kind in ('activity_de','tf_de'):x=x.rename(columns={'adj.P.Val':'padj','logFC':'NES'})
        if cfg['significant_only']:x=x[x.padj<=cfg['alpha']]
        x=x.sort_values('padj').head(int(cfg['top']))
        if 'description' not in x:x['description']=x.term
        return ('dot' if kind in ('go','pathway') else 'barh'),x,('Gene ratio' if kind in ('go','pathway') else 'NES' if kind not in ('activity_de','tf_de') else 'Activity score difference','Gene sets')
    if kind=='gsea_curve':
        from results import curve_data
        database=cfg['database']
        if not cfg['term'] and database=='go' and (p/'gsea.csv').exists() and pd.read_csv(p/'gsea.csv').empty:database='pathway'
        x=curve_data(root,database,cfg['term'],p.name if p!=root else None)
        return 'curve',x,('Gene rank','Running enrichment score')
    if kind=='gene_expression':
        from results import gene_data
        candidates=r.dropna(subset=['padj'])
        if candidates.empty:candidates=r.sort_values('baseMean',ascending=False)
        gene=cfg['gene'] or str(candidates.iloc[0].gene)
        x=gene_data(root,gene,p.name if p!=root else None)
        return 'gene',x,('Group','Normalized count')
    if kind in ('distribution','library_size','mapping','fastp_summary','wgcna_modules','deg_summary','quality','base_content','dispersion'):
        if kind=='distribution':
            if cfg['expression_scale']=='vst':x=v;label='VST expression'
            else:
                path=root/('tpm.csv' if cfg['expression_scale']=='tpm' else 'normalized.csv')
                if not path.exists():raise ValueError('本项目没有TPM；原始计数不含基因长度，不能直接推算TPM')
                x=np.log2(pd.read_csv(path,index_col=0)+1);label='log2(TPM + 1)' if cfg['expression_scale']=='tpm' else 'log2(normalized count + 1)'
            return 'box',x,('Samples',label)
        if kind in ('quality','base_content'):
            sample=cfg['sample'] or m['sample'].iloc[0];f=root/'qc'/sample/'fastp.json'
            if not f.exists():raise ValueError('该样本无 fastp 数据')
            x=json.loads(f.read_text());return ('quality' if kind=='quality' else 'base_content'),x,('Position in read','Mean Phred quality' if kind=='quality' else 'Base fraction')
        if kind=='deg_summary':
            x=pd.DataFrame(list(json.loads((root/'analysis-summary.json').read_text())['contrasts'].values()));return 'deg',x,('Comparison','Significant genes')
        if kind=='wgcna_modules':x=pd.read_csv(root/'wgcna-modules.csv').groupby('module').size().rename('genes').reset_index();return 'modules',x,('Module','Genes')
        file={'library_size':'sample-qc.csv','mapping':'salmon-qc.csv','fastp_summary':'fastp-qc.csv','dispersion':'dispersion.csv'}[kind]
        return kind,pd.read_csv(root/file),('Sample','Reads / fragments' if kind!='mapping' else 'Mapped fragments (%)')
    raise ValueError('未知图表')

def render(p,kind,opts):
    if kind not in TITLES:raise ValueError('未知图表')
    cfg=get_settings(p,kind,opts);root=data_root(p);m=pd.read_csv(root/'samples.csv');groups=colors_for(m,root,cfg)
    category,data,axes=load_plot(p,kind,cfg)
    available={f.name for f in font_manager.fontManager.ttflist};family=next((x for x in ['Noto Sans CJK SC','Noto Sans CJK JP','Microsoft YaHei','WenQuanYi Zen Hei'] if x in available),'DejaVu Sans')
    with plt.rc_context({'font.family':family,'font.size':cfg['font'],'svg.fonttype':'none','pdf.fonttype':42}):
        fig,ax=plt.subplots(figsize=(cfg['width'],cfg['height']),layout='constrained')
        try:
            if category=='scatter':
                ax.scatter(data.plot_x,data.plot_y,c=data.color,s=cfg['point_size'],alpha=cfg['point_alpha'],edgecolors='#5c6962',linewidths=.65,rasterized=False)
                if kind=='volcano':
                    ax.axhline(-np.log10(cfg['alpha']),ls='--',color='#888888');ax.axvline(cfg['lfc'],ls='--',color='#888888');ax.axvline(-cfg['lfc'],ls='--',color='#888888')
                else:ax.set_xscale('log');ax.axhline(0,color='#999999',ls='--')
                chosen=pd.concat([data[data.manual_label],data.nsmallest(int(cfg['labels']),cfg['p_column'])]).drop_duplicates('gene')
                texts=[];fig.canvas.draw()
                for i,(_,row) in enumerate(chosen.iterrows()):
                    label=str(row.symbol) if 'symbol' in row and pd.notna(row.symbol) else row.gene
                    position=cfg['label_positions'].get(row.gene)
                    if position:
                        text=ax.annotate(label,(row.plot_x,row.plot_y),xytext=(position.get('ax',25)*.75,-position.get('ay',-25)*.75),textcoords='offset points',fontsize=cfg['label_font'],fontstyle='italic',bbox=dict(boxstyle='round,pad=.2',facecolor='white',edgecolor='#8b938e',alpha=.95),arrowprops=dict(arrowstyle='-',color='#69756d',lw=.7))
                    else:
                        text=ax.text(row.plot_x,row.plot_y,label,fontsize=cfg['label_font'],fontstyle='italic',bbox=dict(boxstyle='round,pad=.2',facecolor='white',edgecolor='#8b938e',alpha=.95));texts.append(text)
                if texts:np.random.seed(42);adjust_text(texts,x=data.plot_x.to_numpy(),y=data.plot_y.to_numpy(),ax=ax,time_lim=.7,arrowprops=dict(arrowstyle='-',color='#69756d',lw=.7))
                if cfg['grid']:ax.grid(True,color='#e8ece9',lw=.6);ax.set_axisbelow(True)
            elif category=='pca':
                for g in data.group.unique():
                    x=data[data.group==g];ax.scatter(x.plot_x,x.plot_y,color=groups[g],label=g,s=65)
                for _,row in data.iterrows():ax.annotate(row['sample'],(row.plot_x,row.plot_y),fontsize=8,xytext=(4,4),textcoords='offset points')
                ax.legend(frameon=False)
            elif category in ('heatmap','matrix'):
                if data.empty:ax.text(.5,.5,'No eligible genes / sets',ha='center');ax.axis('off')
                else:
                    cmap=LinearSegmentedColormap.from_list('groups',[cfg['down'],'#ffffff',cfg['up']])
                    kwargs={'vmin':-2,'vmax':2} if kind=='heatmap' else {'vmin':-1,'vmax':1} if kind=='wgcna_traits' else {'vmin':0,'vmax':1} if kind=='correlation' else {}
                    image=ax.imshow(data,aspect='auto',cmap=cmap,**kwargs);fig.colorbar(image,ax=ax)
                    ylabels=list(data.index)
                    if kind=='heatmap' and cfg['gene_labels']=='symbol' and (root/'gene-annotation.csv').exists():
                        ann=pd.read_csv(root/'gene-annotation.csv').set_index('gene');ylabels=[str(ann.loc[g,'symbol']) if g in ann.index and pd.notna(ann.loc[g,'symbol']) else g for g in data.index]
                    ax.set_xticks(range(len(data.columns)),data.columns,rotation=50,ha='right');ax.set_yticks(range(len(data)),ylabels,fontsize=8)
                    if kind in ('correlation','wgcna_traits') and data.size<300:
                        for i in range(len(data)):
                            for j in range(len(data.columns)):ax.text(j,i,f'{data.iloc[i,j]:.2f}',ha='center',va='center',fontsize=7)
                    # Group annotation strips follow the clustered sample order.
                    sample_group=dict(zip(m['sample'],m.group))
                    for j,name in enumerate(data.columns):
                        if name in sample_group:ax.plot(j,-.65,marker='s',color=groups[sample_group[name]],markersize=7,clip_on=False)
            elif category in ('dot','barh'):
                if data.empty:ax.text(.5,.5,'No significant sets at the selected threshold',ha='center',va='center');ax.axis('off')
                else:
                    y=np.arange(len(data));labels=[str(x)[:70] for x in data.description]
                    if category=='dot':
                        im=ax.scatter(data.gene_ratio,y,s=30+data.overlap*6,c=-np.log10(data.padj.clip(lower=1e-300)),cmap='YlOrRd');fig.colorbar(im,ax=ax,label='−log10 BH adjusted P')
                    else:ax.barh(y,data.NES,color=[cfg['up'] if a>=0 else cfg['down'] for a in data.NES]);ax.axvline(0,color='#999999')
                    ax.set_yticks(y,labels,fontsize=8);ax.invert_yaxis()
            elif category=='curve':
                fig.clear();ax,hit_ax,rank_ax=fig.subplots(3,1,sharex=True,gridspec_kw={'height_ratios':[3,.5,1]})
                ax.plot(data['rank'],data['running'],color=cfg['up']);ax.axhline(0,color='#999999',ls='--');hit_ax.vlines(np.array(data['hits'])+1,0,1,color=cfg['down'],lw=.6);hit_ax.axis('off');rank_ax.fill_between(data['rank'],data['stat'],color=cfg['neutral']);rank_ax.set_xlabel('Gene rank');rank_ax.set_ylabel('Wald stat')
                ax.text(.02,.02,f"NES={data['result']['NES']:.3g}; padj={data['result']['padj']:.3g}",transform=ax.transAxes,fontsize=9)
            elif category=='gene':
                x=pd.DataFrame(data['samples']);rng=np.random.default_rng(42)
                for i,g in enumerate(x.group.unique()):
                    values=x[x.group==g].normalized;ax.boxplot([values],positions=[i],widths=.45,showfliers=False);ax.scatter(i+rng.uniform(-.13,.13,len(values)),values,c=groups[g],s=45,zorder=3)
                ax.set_xticks(range(x.group.nunique()),list(x.group.unique()));ax.set_title(data['gene'])
            elif category=='box':
                ax.boxplot([data[s].to_numpy() for s in m['sample']],tick_labels=m['sample'],showfliers=False);ax.tick_params(axis='x',labelrotation=45)
            elif category=='library_size':ax.bar(data['sample'],data.library_size,color=[groups[g] for g in data.group]);ax.tick_params(axis='x',labelrotation=45)
            elif category=='mapping':
                if 'percent_mapped' in data:ax.bar(data['sample'],data.percent_mapped,color=cfg['down']);ax.set_ylim(0,100)
                else:ax.text(.5,.5,'Mapping metadata unavailable',ha='center');ax.axis('off')
            elif category=='fastp_summary':
                x=np.arange(len(data));ax.bar(x-.2,data.raw_reads,width=.4,label='Raw reads',color=cfg['neutral']);ax.bar(x+.2,data.clean_reads,width=.4,label='Clean reads',color=cfg['down']);ax.set_xticks(x,data['sample'],rotation=45);ax.legend(frameon=False)
            elif category=='quality':
                for key,color,label in [('read1_before_filtering',cfg['neutral'],'Raw R1'),('read1_after_filtering',cfg['down'],'Clean R1')]:
                    values=data.get(key,{}).get('quality_curves',{}).get('mean',[])
                    if values:ax.plot(range(1,len(values)+1),values,label=label,color=color)
                ax.legend(frameon=False)
            elif category=='base_content':
                for base,color in [('A','#126657'),('T','#b64c4c'),('G','#bd9a4b'),('C','#317ba0')]:
                    values=data.get('read1_before_filtering',{}).get('content_curves',{}).get(base,[])
                    if values:ax.plot(range(1,len(values)+1),values,color=color,label=base)
                ax.legend(frameon=False);ax.set_ylim(0,1)
            elif category=='deg':
                x=np.arange(len(data));ax.bar(x-.2,data.up,width=.4,color=cfg['up'],label='Up');ax.bar(x+.2,data.down,width=.4,color=cfg['down'],label='Down');ax.set_xticks(x,data.name,rotation=30);ax.legend(frameon=False)
            elif category=='soft':
                y=-np.sign(data.slope)*data['SFT.R.sq'];ax.plot(data.Power,y,'o-',color=cfg['up']);ax.axhline(.85,ls='--',color='#89958d');ax.set_ylim(-1,1)
                twin=ax.twinx();twin.plot(data.Power,data['mean.k.'],'s--',color=cfg['down']);twin.set_ylabel('Mean connectivity')
            elif category=='modules':ax.bar(data.module,data.genes,color=module_colors(data.module,cfg['neutral']));ax.tick_params(axis='x',labelrotation=45)
            elif category=='dispersion':
                data=data[(data.baseMean>0)&(data.dispersion>0)];ax.scatter(data.baseMean,data.dispersion,s=8,c=cfg['down'],alpha=.5);ax.set_xscale('log');ax.set_yscale('log');axes=('Mean normalized expression','Dispersion')
            for axis in ('x','y'):
                if cfg[axis+'_min']!='' or cfg[axis+'_max']!='':
                    limits=getattr(ax,'get_'+axis+'lim')()
                    getattr(ax,'set_'+axis+'lim')(float(cfg[axis+'_min']) if cfg[axis+'_min']!='' else limits[0],float(cfg[axis+'_max']) if cfg[axis+'_max']!='' else limits[1])
            ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
            ax.set_xlabel(cfg['xlabel'] or axes[0]);ax.set_ylabel(cfg['ylabel'] or axes[1]);ax.set_title(cfg['title'] or TITLES[kind])
            for fmt in ('svg','pdf','png'):fig.savefig(p/(kind+'.'+fmt),dpi=int(cfg['dpi']))
            if isinstance(data,pd.DataFrame):data.to_csv(p/(kind+'.plot-data.csv'))
            else:(p/(kind+'.plot-data.json')).write_text(json.dumps(data,ensure_ascii=False))
            (p/(kind+'.settings.json')).write_text(json.dumps(cfg,indent=2,ensure_ascii=False))
        finally:plt.close(fig)

def interactive_data(p,kind,args):
    from results import contrast_dir, records
    target=contrast_dir(p,args.get('contrast'));cfg=get_settings(target,kind)
    for key in ('term','gene','database','sample','annotated_genes'):
        if args.get(key):cfg[key]=args[key]
    category,data,axes=load_plot(target,kind,cfg);root=data_root(target);m=pd.read_csv(root/'samples.csv');groups=colors_for(m,root,cfg)
    layout=dict(title=cfg['title'] or TITLES[kind],xaxis=dict(title=axes[0]),yaxis=dict(title=axes[1]),margin=dict(l=90,r=40,t=65,b=100),font=dict(size=cfg['font']),paper_bgcolor='#ffffff',plot_bgcolor='#ffffff')
    traces=[]
    if category in ('scatter','pca'):
        if category=='scatter':
            traces=[dict(type='scattergl',mode='markers',x=data.plot_x.tolist(),y=data.plot_y.tolist(),text=(data['symbol'].fillna(data.gene) if 'symbol' in data else data.gene).tolist(),customdata=data.gene.tolist(),marker=dict(color=data.color.tolist(),size=np.where(data.manual_label,13,max(5,np.sqrt(cfg['point_size']))).tolist(),opacity=cfg['point_alpha'],line=dict(width=.8,color='#526158')),hovertemplate='%{text}<br>x=%{x:.3g}<br>y=%{y:.3g}<extra></extra>')]
            chosen=pd.concat([data[data.manual_label],data.nsmallest(int(cfg['labels']),cfg['p_column'])]).drop_duplicates('gene')
            layout['annotations']=[]
            for i,(_,row) in enumerate(chosen.iterrows()):
                pos=cfg['label_positions'].get(row.gene,dict(ax=25 if i%2==0 else -25,ay=-25-18*(i%3)))
                label=str(row.symbol) if 'symbol' in row and pd.notna(row.symbol) else row.gene
                layout['annotations'].append(dict(name=row.gene,x=float(row.plot_x),y=float(row.plot_y),text='<i>'+html.escape(label)+'</i>',showarrow=True,arrowhead=0,arrowcolor='#68756d',ax=pos.get('ax',25),ay=pos.get('ay',-25),bgcolor='white',bordercolor='#919b94',borderwidth=1,borderpad=3,font=dict(size=cfg['label_font'],color='#26382e')))
            layout['xaxis']['showgrid']=cfg['grid'];layout['yaxis']['showgrid']=cfg['grid']
            if kind=='volcano':
                layout['shapes']=[dict(type='line',x0=cfg['lfc'],x1=cfg['lfc'],y0=0,y1=1,yref='paper',line=dict(dash='dash',color='#77857c',width=1)),dict(type='line',x0=-cfg['lfc'],x1=-cfg['lfc'],y0=0,y1=1,yref='paper',line=dict(dash='dash',color='#77857c',width=1)),dict(type='line',x0=0,x1=1,xref='paper',y0=-np.log10(cfg['alpha']),y1=-np.log10(cfg['alpha']),line=dict(dash='dash',color='#77857c',width=1))]
            if kind=='ma':layout['xaxis']['type']='log'
        else:
            for g in data.group.unique():
                x=data[data.group==g];traces.append(dict(type='scatter',mode='markers+text',name=g,x=x.plot_x.tolist(),y=x.plot_y.tolist(),text=x['sample'].tolist(),textposition='top center',marker=dict(color=groups[g],size=11)))
    elif category in ('heatmap','matrix'):
        traces=[dict(type='heatmap',z=data.to_numpy().tolist(),x=list(data.columns),y=list(data.index),colorscale=[[0,cfg['down']],[.5,'#ffffff'],[1,cfg['up']]],hovertemplate='%{y}<br>%{x}: %{z:.4g}<extra></extra>')]
        if kind=='heatmap':
            if data.empty:layout['annotations']=[dict(x=.5,y=.5,xref='paper',yref='paper',text='所选方案没有可展示的差异基因；可切换方案或指定基因',showarrow=False)]
            traces[0].update(zmin=-2,zmax=2,customdata=[[g]*len(data.columns) for g in data.index])
            if cfg['gene_labels']=='symbol' and (root/'gene-annotation.csv').exists():
                ann=pd.read_csv(root/'gene-annotation.csv').set_index('gene');traces[0]['y']=[str(ann.loc[g,'symbol']) if g in ann.index and pd.notna(ann.loc[g,'symbol']) else g for g in data.index]
            traces[0]['hovertemplate']='%{y}<br>%{customdata}<br>%{x}: %{z:.4g}<extra></extra>'
        if kind=='correlation':traces[0].update(zmin=0,zmax=1)
        if kind=='activity':traces[0]['customdata']=[[g]*len(data.columns) for g in data.index]
        if kind=='wgcna_traits':traces[0].update(zmin=-1,zmax=1)
        layout['yaxis']['autorange']='reversed'
    elif category in ('dot','barh'):
        if category=='dot':traces=[dict(type='scatter',mode='markers',x=data.gene_ratio.tolist(),y=data.description.tolist(),customdata=data.term.tolist(),marker=dict(size=(10+np.sqrt(data.overlap)*5).tolist(),color=(-np.log10(data.padj.clip(lower=1e-300))).tolist(),colorscale='YlOrRd',showscale=True))]
        else:traces=[dict(type='bar',orientation='h',x=data.NES.tolist(),y=data.description.tolist(),customdata=data.term.tolist(),marker=dict(color=[cfg['up'] if x>=0 else cfg['down'] for x in data.NES]))]
        layout['yaxis']['autorange']='reversed';layout['margin']['l']=240
    elif category=='curve':
        traces=[dict(type='scatter',mode='lines',name=data['term'],x=data['rank'],y=data['running'],line=dict(color=cfg['up'])),dict(type='scatter',mode='markers',name='Gene-set hits',x=(np.array(data['hits'])+1).tolist(),y=[0]*len(data['hits']),text=[data['gene'][i] for i in data['hits']],marker=dict(symbol='line-ns-open',size=14,color=cfg['down']))]
        layout['title']=data['term']+f" · NES {data['result']['NES']:.3g}, padj {data['result']['padj']:.3g}"
    elif category=='gene':
        frame=pd.DataFrame(data['samples'])
        for g in frame.group.unique():
            x=frame[frame.group==g];traces.append(dict(type='box',name=g,y=x.normalized.tolist(),text=x['sample'].tolist(),boxpoints='all',jitter=.3,marker=dict(color=groups[g]),hovertemplate='%{text}<br>%{y:.3g}<extra></extra>'))
        layout['title']=data['gene']
    elif category=='box':
        for s in m['sample']:traces.append(dict(type='box',name=s,y=data[s].tolist(),boxpoints=False,marker=dict(color=groups[m.set_index('sample').loc[s,'group']])))
    elif category=='quality':
        for key,color,label in [('read1_before_filtering',cfg['neutral'],'Raw R1'),('read1_after_filtering',cfg['down'],'Clean R1')]:
            values=data.get(key,{}).get('quality_curves',{}).get('mean',[]);traces.append(dict(type='scatter',mode='lines',name=label,x=list(range(1,len(values)+1)),y=values,line=dict(color=color)))
    elif category=='base_content':
        for base,color in [('A','#126657'),('T','#b64c4c'),('G','#bd9a4b'),('C','#317ba0')]:
            values=data.get('read1_before_filtering',{}).get('content_curves',{}).get(base,[])
            traces.append(dict(type='scatter',mode='lines',name=base,x=list(range(1,len(values)+1)),y=values,line=dict(color=color)))
    elif category=='library_size':traces=[dict(type='bar',x=data['sample'].tolist(),y=data.library_size.tolist(),marker=dict(color=[groups[g] for g in data.group]))]
    elif category=='mapping':traces=[dict(type='bar',x=data['sample'].tolist(),y=data.percent_mapped.tolist() if 'percent_mapped' in data else [],marker=dict(color=cfg['down']))]
    elif category=='fastp_summary':traces=[dict(type='bar',name=label,x=data['sample'].tolist(),y=data[col].tolist(),marker=dict(color=color)) for col,label,color in [('raw_reads','Raw reads',cfg['neutral']),('clean_reads','Clean reads',cfg['down'])]]
    elif category=='deg':traces=[dict(type='bar',name=label,x=data.name.tolist(),y=data[col].tolist(),marker=dict(color=color)) for col,label,color in [('up','Up',cfg['up']),('down','Down',cfg['down'])]]
    elif category=='soft':traces=[dict(type='scatter',mode='lines+markers',x=data.Power.tolist(),y=(-np.sign(data.slope)*data['SFT.R.sq']).tolist(),name='Signed scale-free fit',line=dict(color=cfg['up']))];layout['yaxis']['range']=[-1,1]
    elif category=='modules':traces=[dict(type='bar',x=data.module.tolist(),y=data.genes.tolist(),marker=dict(color=module_colors(data.module,cfg['neutral']))) ]
    elif category=='dispersion':
        x=data[(data.baseMean>0)&(data.dispersion>0)];traces=[dict(type='scattergl',mode='markers',x=x.baseMean.tolist(),y=x.dispersion.tolist(),marker=dict(color=cfg['down'],size=4))];layout['xaxis']['type']='log';layout['yaxis']['type']='log'
    if cfg['xlabel']:layout['xaxis']['title']=cfg['xlabel']
    if cfg['ylabel']:layout['yaxis']['title']=cfg['ylabel']
    if cfg['title']:layout['title']=cfg['title']
    for axis in ('x','y'):
        if cfg[axis+'_min']!='' and cfg[axis+'_max']!='':
            limits=[float(cfg[axis+'_min']),float(cfg[axis+'_max'])]
            if layout[axis+'axis'].get('type')=='log':
                if min(limits)<=0:raise ValueError('对数轴坐标范围必须大于0')
                limits=np.log10(limits).tolist()
            layout[axis+'axis']['range']=limits
    return dict(data=traces,layout=layout,config=dict(responsive=True,displaylogo=False,editable=True,edits=dict(annotationPosition=True,annotationText=False,axisTitleText=False,titleText=False,legendPosition=True),toImageButtonOptions=dict(format='svg',filename=kind)),settings=cfg,category=category,metadata=data if category in ('gene','curve') else None)
