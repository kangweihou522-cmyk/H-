from pathlib import Path
import json, re, os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).parent/'.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

def render(p,kind,opts):
 if kind not in ('pca','volcano','heatmap','go','gsea','pathway','pathway_gsea'):raise ValueError('未知图表')
 saved=p/(kind+'.settings.json');old=json.loads(saved.read_text()) if saved.exists() else {}
 analysis=json.loads((p/'config.json').read_text()) if (p/'config.json').exists() else {}
 cfg={**dict(title='',width=8,height=6,font=11,dpi=300,up='#b64c4c',down='#317ba0',neutral='#b8c2c0',alpha=analysis.get('alpha',.05),lfc=analysis.get('lfc',1),top=20,xlabel='',ylabel=''),**old,**opts}
 for key,lo,hi in [('width',3,20),('height',3,20),('font',6,30),('dpi',72,600),('top',5,60),('alpha',.000001,.999),('lfc',0,10)]:
  cfg[key]=float(cfg[key])
  if not lo<=cfg[key]<=hi:raise ValueError('图表参数超出范围：'+key)
 for key in ('up','down','neutral'):
  if not re.fullmatch('#[0-9a-fA-F]{6}',cfg[key]):raise ValueError('颜色需要 #RRGGBB')
 for key in ('title','xlabel','ylabel'):
  if not isinstance(cfg[key],str) or len(cfg[key])>180:raise ValueError('标题过长')
 available={f.name for f in font_manager.fontManager.ttflist}
 family=next((n for n in ['Noto Sans CJK SC','Noto Sans CJK JP','Microsoft YaHei','WenQuanYi Zen Hei'] if n in available),'DejaVu Sans')
 with plt.rc_context({'font.family':family,'font.size':cfg['font'],'svg.fonttype':'none','pdf.fonttype':42}):
  fig,ax=plt.subplots(figsize=(cfg['width'],cfg['height']),layout='constrained')
  try:
   r=pd.read_csv(p/'differential.csv');m=pd.read_csv(p/'samples.csv');v=pd.read_csv(p/'vst.csv',index_col=0)
   if kind=='volcano':
    good=r.dropna(subset=['padj','log2FoldChange']).copy();x=good.log2FoldChange;y=-np.log10(good.padj.clip(lower=1e-300))
    colors=np.where((good.padj<cfg['alpha'])&(x>=cfg['lfc']),cfg['up'],np.where((good.padj<cfg['alpha'])&(x<=-cfg['lfc']),cfg['down'],cfg['neutral']))
    ax.scatter(x,y,c=colors,s=16,alpha=.8);ax.axhline(-np.log10(cfg['alpha']),ls='--',c='#777777');ax.axvline(cfg['lfc'],ls='--',c='#777777');ax.axvline(-cfg['lfc'],ls='--',c='#777777')
    for i in good.nsmallest(min(10,int(cfg['top'])),'padj').index:
     ax.annotate(str(r.loc[i,'symbol']) if pd.notna(r.loc[i,'symbol']) else str(r.loc[i,'gene']),(x.loc[i],y.loc[i]),fontsize=8)
    ax.set(xlabel='log2 fold change',ylabel='−log10 adjusted P');data=good
   elif kind=='pca':
    a=v[m['sample']].T.to_numpy(dtype=float,copy=True);a-=a.mean(axis=0);u,s,_=np.linalg.svd(a,full_matrices=False);scores=u[:,:2]*s[:2];variance=s*s/(s*s).sum()
    for i,g in enumerate(m.group.unique()):
     mask=m.group==g;ax.scatter(scores[mask,0],scores[mask,1],label=g,c=cfg['up'] if i==0 else cfg['down'],s=65)
    for i,name in enumerate(m['sample']):ax.annotate(name,scores[i],xytext=(5,5),textcoords='offset points',fontsize=8)
    ax.legend(frameon=False);ax.set(xlabel=f'PC1 ({variance[0]:.1%})',ylabel=f'PC2 ({variance[1]:.1%})');data=pd.DataFrame({'sample':m['sample'],'group':m.group,'PC1':scores[:,0],'PC2':scores[:,1]})
   elif kind=='heatmap':
    genes=r.dropna(subset=['padj']).nsmallest(int(cfg['top']),'padj').gene;data=v.loc[v.index.intersection(genes),m['sample']]
    if data.empty:ax.text(.5,.5,'No genes with adjusted P',ha='center');ax.axis('off')
    else:
     data=data.sub(data.mean(axis=1),axis=0).div(data.std(axis=1).replace(0,1),axis=0)
     im=ax.imshow(data,aspect='auto',cmap=LinearSegmentedColormap.from_list('groups',[cfg['down'],'#ffffff',cfg['up']]),vmin=-2,vmax=2)
     ax.set_xticks(range(len(data.columns)),data.columns,rotation=45,ha='right');ax.set_yticks(range(len(data)),data.index,fontsize=8);fig.colorbar(im,ax=ax,label='Row Z-score (VST)')
   else:
    data=pd.read_csv(p/(kind+'.csv')).sort_values('padj').head(int(cfg['top']))
    if data.empty:ax.text(.5,.5,'No eligible gene sets',ha='center',va='center');ax.axis('off')
    else:
     x=-np.log10(data.padj.clip(lower=1e-300)) if kind in ('go','pathway') else data.NES
     ax.barh(np.arange(len(data)),x,color=[cfg['up'] if a>=0 else cfg['down'] for a in x]);ax.set_yticks(range(len(data)),[str(s)[:65] for s in data.description],fontsize=8);ax.invert_yaxis();ax.set_xlabel('−log10 BH adjusted P' if kind in ('go','pathway') else 'Normalized enrichment score')
   ax.set_title(cfg['title'] or {'volcano':'Differential expression','pca':'Sample PCA','heatmap':'Top differential genes','go':'GO over-representation','gsea':'GO ranked GSEA','pathway':'KEGG / GMT enrichment','pathway_gsea':'KEGG / GMT GSEA'}[kind])
   if cfg['xlabel']:ax.set_xlabel(cfg['xlabel'])
   if cfg['ylabel']:ax.set_ylabel(cfg['ylabel'])
   for fmt in ('svg','pdf','png'):fig.savefig(p/(kind+'.'+fmt),dpi=int(cfg['dpi']))
   data.to_csv(p/(kind+'.plot-data.csv'),index=True);saved.write_text(json.dumps(cfg,indent=2))
  finally:plt.close(fig)
