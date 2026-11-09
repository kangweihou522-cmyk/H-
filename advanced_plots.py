"""Multi-comparison figures use explicit, consistent profiles and original statistics."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from dataio import read_frame

LABELS={'composition':'参考辅助RNA组成贡献', 'contrast_effects':'多比较基因效应热图','gene_forest':'跨比较基因效应与区间','pathway_comparison':'多比较通路效应','deg_upset':'跨比较DEG精确交集'}


def data(p,kind,cfg):
    if kind=='composition':
        r=read_frame(p/'composition-fractions.csv');return r,[c for c in r if c!='sample']
    targets=[q for q in sorted((p/'contrasts').glob('*')) if (q/'differential.csv').exists()]
    if kind=='deg_upset':
        sets={};profile=cfg['comparison_profile']
        for q in targets:
            f=q/'profiles'/profile/'differential.csv'
            if not f.exists():continue
            r=read_frame(f);mask=r.significant==True
            if cfg['direction']!='all':mask &= r.log2FoldChange.gt(0) if cfg['direction']=='up' else r.log2FoldChange.lt(0)
            sets[q.name]=set(r.loc[mask,'gene'])
        rows=[];union=set().union(*sets.values()) if sets else set();counts={}
        for gene in sorted(union):
            membership=tuple(name for name,genes in sets.items() if gene in genes);counts.setdefault(membership,[]).append(gene)
        for membership,genes in sorted(counts.items(),key=lambda x:(-len(x[1]),x[0]))[:int(cfg['top'])]:rows.append(dict(intersection=' & '.join(membership),count=len(genes),members=';'.join(membership),genes=';'.join(genes),profile=profile,direction=cfg['direction']))
        return pd.DataFrame(rows,columns=['intersection','count','members','genes','profile','direction']),list(sets)
    if kind in ('contrast_effects','gene_forest'):
        frames=[]
        for q in targets:
            r=read_frame(q/'differential.csv');r['comparison']=q.name;frames.append(r)
        r=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
        if r.empty:return r,[]
        from plots import resolve_genes
        genes=resolve_genes(p,None,cfg['genes'] or cfg['gene']) if cfg['genes'] or cfg['gene'] else r.groupby('gene').padj.min().nsmallest(int(cfg['top'])).index.tolist()
        r=r[r.gene.isin(genes)].copy()
        if kind=='gene_forest':r=r[r.gene==genes[0]] if genes else r.iloc[:0]
        return r,list(r.comparison.unique())
    metric=cfg['comparison_metric'];file={'gsea':'gsea.csv' if cfg['database']=='go' else 'pathway_gsea.csv','activity':'pathway-activity-differential.csv','tf':'tf-activity-differential.csv'}[metric]
    frames=[]
    for q in sorted((p/'contrasts').glob('*')):
        if not (q/file).exists():continue
        x=read_frame(q/file).rename(columns={'logFC':'NES','adj.P.Val':'padj'});x['comparison']=q.name;frames.append(x)
    r=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
    if len(r):
        terms=r.groupby('term').padj.min().nsmallest(int(cfg['top'])).index;r=r[r.term.isin(terms)]
    return r,list(r.comparison.unique()) if len(r) else []


def spec(p,kind,cfg):
    r,names=data(p,kind,cfg);traces=[];layout=dict(title=dict(text=cfg['title'] or {'contrast_effects':'Gene effects across comparisons','gene_forest':'Raw gene effect and Wald interval','pathway_comparison':'Pathway effects across comparisons','deg_upset':'Exact DEG intersections','composition':'Reference-assisted RNA contributions'}[kind]),font=dict(size=cfg['font']),margin=dict(l=135,r=70,t=65,b=100),paper_bgcolor='white',plot_bgcolor='white',showlegend=False)
    metadata={'kind':kind,'rows':json.loads(r.to_json(orient='records')),'comparisons':names}
    if r.empty:layout['annotations']=[dict(text='无符合条件的结果；请检查比较、方案及输入',xref='paper',yref='paper',x=.5,y=.5,showarrow=False)]
    elif kind=='composition':
        from plots import PALETTE
        for i,name in enumerate(names):traces.append(dict(type='bar',name=name,x=r['sample'].tolist(),y=r[name].tolist(),marker=dict(color=PALETTE[i%len(PALETTE)])))
        layout.update(barmode='stack',showlegend=True,yaxis=dict(title=dict(text='Conditional RNA contribution'),range=[0,1]))
    elif kind=='deg_upset':
        x=list(range(len(r)));traces.append(dict(type='bar',x=x,y=r['count'].tolist(),marker=dict(color=cfg['down']),customdata=r.genes.tolist(),hovertemplate='%{y} genes<br>%{customdata}<extra></extra>'))
        for j,row in enumerate(r.itertuples()):
            present=[i for i,n in enumerate(names) if n in row.members.split(';')]
            traces.append(dict(type='scatter',x=[j]*len(names),y=list(range(len(names))),mode='markers',marker=dict(color=[cfg['down'] if i in present else '#e0e5e2' for i in range(len(names))],size=9),yaxis='y2',hoverinfo='skip'))
            if len(present)>1:traces.append(dict(type='scatter',x=[j,j],y=[min(present),max(present)],mode='lines',line=dict(color=cfg['down'],width=2),yaxis='y2',hoverinfo='skip'))
        layout.update(xaxis=dict(tickvals=x,ticktext=[str(i+1) for i in x]),yaxis=dict(domain=[.38,1],title=dict(text='Exact intersection genes')),yaxis2=dict(domain=[0,.3],tickvals=list(range(len(names))),ticktext=names,autorange='reversed'),annotations=[dict(text=f"{cfg['comparison_profile']} · {cfg['direction']} · exact intersections",xref='paper',yref='paper',x=.5,y=1.08,showarrow=False)])
    elif kind=='gene_forest':
        lo=r.get('ci_low_raw',pd.Series(np.nan,index=r.index));hi=r.get('ci_high_raw',pd.Series(np.nan,index=r.index));lfc=r.get('log2FoldChange_raw',r.log2FoldChange)
        traces=[dict(type='scatter',mode='markers',x=lfc.tolist(),y=r.comparison.tolist(),customdata=r.gene.tolist(),marker=dict(color=cfg['down'],size=10),error_x=dict(type='data',symmetric=False,array=(hi-lfc).tolist(),arrayminus=(lfc-lo).tolist()),hovertemplate='%{y}<br>raw log2FC=%{x:.4g}<extra></extra>')]
        layout.update(xaxis=dict(title=dict(text='Raw log2 fold change; Wald 95% interval')),shapes=[dict(type='line',x0=0,x1=0,y0=0,y1=1,yref='paper',line=dict(color='#aaaaaa',dash='dash'))])
    else:
        key='gene' if kind=='contrast_effects' else 'term';col='log2FoldChange' if key=='gene' else 'NES';matrix=r.pivot(index=key,columns='comparison',values=col);q=r.pivot(index=key,columns='comparison',values='padj').reindex_like(matrix)
        traces=[dict(type='heatmap',x=matrix.columns.tolist(),y=matrix.index.tolist(),z=matrix.to_numpy().tolist(),zmid=0,colorscale=[[0,cfg['down']],[.5,'#ffffff'],[1,cfg['up']]],customdata=[[str(g)]*len(matrix.columns) for g in matrix.index],text=[[f'q={v:.3g}' for v in row] for row in q.to_numpy()],hovertemplate='%{y}<br>%{x}: %{z:.4g}<br>%{text}<extra></extra>',colorbar=dict(title=dict(text='log2FC' if key=='gene' else 'NES' if cfg['comparison_metric']=='gsea' else 'Score difference')))]
        layout['yaxis']=dict(autorange='reversed');metadata['matrix']=matrix.to_dict()
    return dict(data=traces,layout=layout,metadata=metadata,category='advanced',settings=cfg,config=dict(responsive=True,displaylogo=False),note='多比较结果来自同一研究，不是独立验证。UpSet使用同一指定方案和方向；效应热图保留全部显示基因/通路，不以颜色代替显著性。森林图使用原始Wald效应区间；旧任务没有原始SE时不伪造区间。')


def render(p,kind,cfg):
    from plots import plt
    from matplotlib.gridspec import GridSpec
    r,names=data(p,kind,cfg);fig=plt.figure(figsize=(cfg['width'],cfg['height']),layout='constrained')
    if r.empty:ax=fig.add_subplot();ax.text(.5,.5,'No eligible results',ha='center',transform=ax.transAxes)
    elif kind=='composition':
        from plots import PALETTE
        ax=fig.add_subplot();bottom=np.zeros(len(r))
        for i,name in enumerate(names):ax.bar(r['sample'],r[name],bottom=bottom,label=name,color=PALETTE[i%len(PALETTE)]);bottom+=r[name].to_numpy()
        ax.set_ylim(0,1);ax.set_ylabel('Conditional RNA contribution');ax.tick_params(axis='x',rotation=45);ax.legend(fontsize=8)
    elif kind=='deg_upset':
        gs=GridSpec(2,1,figure=fig,height_ratios=[3,1.5]);ax=fig.add_subplot(gs[0]);dot=fig.add_subplot(gs[1],sharex=ax);x=np.arange(len(r));ax.bar(x,r['count'],color=cfg['down']);ax.set_ylabel('Exact intersection genes')
        for j,row in enumerate(r.itertuples()):
            present=[i for i,n in enumerate(names) if n in row.members.split(';')];dot.scatter([j]*len(names),range(len(names)),c=[cfg['down'] if i in present else '#e0e5e2' for i in range(len(names))],s=25)
            if len(present)>1:dot.plot([j,j],[min(present),max(present)],color=cfg['down'])
        dot.set_yticks(range(len(names)),names);dot.invert_yaxis();dot.set_xticks(x,[str(i+1) for i in x]);dot.set_xlabel(cfg['comparison_profile']+' / '+cfg['direction']);ax.tick_params(labelbottom=False)
    elif kind=='gene_forest':
        ax=fig.add_subplot();lfc=r.get('log2FoldChange_raw',r.log2FoldChange);ax.scatter(lfc,r.comparison,color=cfg['down']);ax.axvline(0,color='#aaa',ls='--')
        if 'ci_low_raw' in r:ax.errorbar(lfc,r.comparison,xerr=np.stack([lfc-r.ci_low_raw,r.ci_high_raw-lfc]),fmt='none',ecolor=cfg['down'])
        ax.set_xlabel('Raw log2 fold change; Wald 95% interval')
    else:
        from plots import palette
        ax=fig.add_subplot();key='gene' if kind=='contrast_effects' else 'term';col='log2FoldChange' if key=='gene' else 'NES';x=r.pivot(index=key,columns='comparison',values=col);limit=max(.01,float(np.nanmax(np.abs(x.to_numpy()))));im=ax.imshow(x,cmap=palette(cfg)[0],vmin=-limit,vmax=limit,aspect='auto');fig.colorbar(im,ax=ax);ax.set_xticks(range(len(x.columns)),x.columns,rotation=40,ha='right');ax.set_yticks(range(len(x)),x.index)
    for axis in fig.axes:
        axis.tick_params(labelsize=cfg['label_font']);axis.xaxis.label.set_size(cfg['font']);axis.yaxis.label.set_size(cfg['font'])
    fig.suptitle(cfg['title'] or {'contrast_effects':'Gene effects across comparisons','gene_forest':'Raw gene effect and Wald interval','pathway_comparison':'Pathway effects across comparisons','deg_upset':'Exact DEG intersections','composition':'Reference-assisted RNA contributions'}[kind],fontsize=cfg['font']+3)
    from plots import font_manager
    from matplotlib.text import Text
    available={f.name for f in font_manager.fontManager.ttflist};family=next((name for name in ('Noto Sans CJK SC','Noto Sans CJK JP','Microsoft YaHei','WenQuanYi Zen Hei') if name in available),'DejaVu Sans')
    for text in fig.findobj(Text):text.set_fontfamily(family)
    with plt.rc_context({'svg.fonttype':'none','pdf.fonttype':42}):
        for fmt in ('svg','pdf','png'):fig.savefig(p/(kind+'.'+fmt),dpi=cfg['dpi'])
    plt.close(fig);r.to_csv(p/(kind+'.plot-data.csv'),index=False);(p/(kind+'.settings.json')).write_text(json.dumps(cfg,ensure_ascii=False,indent=2))
