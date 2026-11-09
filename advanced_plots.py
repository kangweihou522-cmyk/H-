"""Multi-comparison figures use explicit, consistent profiles and original statistics."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from dataio import read_frame

LABELS={'composition':'参考辅助RNA组成贡献', 'contrast_effects':'多比较基因效应热图','gene_forest':'跨比较基因效应与区间','pathway_comparison':'多比较通路效应','deg_upset':'跨比较DEG精确交集'}


def comparison_targets(p):
    found={q.name:q for q in (p/'contrasts').glob('*') if q.is_dir()}
    config=json.loads((p/'config.json').read_text()) if (p/'config.json').exists() else {}
    names=[c['name'] for c in config.get('comparisons',[]) if c.get('name') in found]
    names+=sorted(set(found)-set(names))
    return [found[name] for name in names]


def data(p,kind,cfg):
    if kind=='composition':
        r=read_frame(p/'composition-fractions.csv');return r,[c for c in r if c!='sample']
    targets=[q for q in comparison_targets(p) if (q/'differential.csv').exists()]
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
        if cfg['significant_only'] and kind=='contrast_effects':genes=[g for g in genes if r.loc[r.gene==g,'padj'].lt(cfg['alpha']).any()]
        r=r[r.gene.isin(genes)].copy();r.attrs['feature_order']=list(genes)
        if kind=='gene_forest':r=r[r.gene==genes[0]] if genes else r.iloc[:0]
        return r,list(r.comparison.unique())
    metric=cfg['comparison_metric'];file={'gsea':'gsea.csv' if cfg['database']=='go' else 'pathway_gsea.csv','activity':'pathway-activity-differential.csv','tf':'tf-activity-differential.csv','footprint':'footprint-activity-differential.csv'}[metric]
    frames=[]
    for q in comparison_targets(p):
        if not (q/file).exists():continue
        x=read_frame(q/file).rename(columns={'logFC':'NES','adj.P.Val':'padj'});x['comparison']=q.name;frames.append(x)
    r=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
    if len(r):
        minimum=r.groupby('term').padj.min();minimum=minimum[minimum<cfg['alpha']] if cfg['significant_only'] else minimum
        terms=minimum.nsmallest(int(cfg['top'])).index;r=r[r.term.isin(terms)].copy();r.attrs['feature_order']=terms.tolist()
    return r,list(r.comparison.unique()) if len(r) else []


TITLES={'contrast_effects':'Gene effects across comparisons','gene_forest':'Raw gene effect and Wald interval','pathway_comparison':'Pathway effects across comparisons','deg_upset':'Exact DEG intersections','composition':'Reference-assisted RNA contributions'}


def comparison_labels(p,names):
    from comparison_context import resolve,label
    cfg=json.loads((p/'config.json').read_text()) if (p/'config.json').exists() else {}
    return [label(resolve(cfg,n),cfg.get('interaction_factor') or cfg.get('display_factor','')) if resolve(cfg,n) else str(n) for n in names]


def matrix_data(r,kind,cfg):
    key='gene' if kind=='contrast_effects' else 'term';effect='log2FoldChange' if key=='gene' else 'NES'
    x=r.pivot(index=key,columns='comparison',values=effect).reindex(columns=r.comparison.drop_duplicates());q=r.pivot(index=key,columns='comparison',values='padj')
    order=[g for g in r.attrs.get('feature_order',[]) if g in x.index]
    if order:x=x.reindex(order)
    q=q.reindex_like(x)
    finite=x.to_numpy()[np.isfinite(x.to_numpy())];extent=max(.01,float(np.max(np.abs(finite)))) if len(finite) else 1
    lo=float(cfg['color_min']) if cfg['color_min']!='' else -extent;hi=float(cfg['color_max']) if cfg['color_max']!='' else extent
    if lo>=hi:raise ValueError('色阶范围不包含有效区间，请同时设置最小/最大值')
    return x,q,lo,hi


def axis_labels(p,matrix,kind,cfg):
    from plots import gene_labels
    y=gene_labels(p,list(matrix.index),cfg) if kind=='contrast_effects' else list(matrix.index)
    return comparison_labels(p,list(matrix.columns)),y


def forest_values(r):
    center=r.get('log2FoldChange_raw',r.log2FoldChange)
    lo=r.get('ci_low_raw',pd.Series(np.nan,index=r.index));hi=r.get('ci_high_raw',pd.Series(np.nan,index=r.index))
    if ((lo.notna()&center.notna()&(lo>center))|(hi.notna()&center.notna()&(hi<center))).any():raise ValueError('原始Wald区间与效应估计不一致，请核对导入结果')
    return center,lo,hi


def configure_layout(layout,cfg,reverse_y=False):
    for axis in ('x','y'):
        values=layout.setdefault(axis+'axis',{});values.update(showgrid=cfg['grid'],gridcolor='#e8ece9',tickfont=dict(size=cfg['label_font']),automargin=True)
        if cfg[axis+'label']:values['title']=dict(text=cfg[axis+'label'])
        low,high=cfg[axis+'_min'],cfg[axis+'_max']
        if low!='' or high!='':
            limits=[float(low) if low!='' else None,float(high) if high!='' else None]
            values['range']=limits[::-1] if axis=='y' and reverse_y else limits;values['autorange']=False


def spec(p,kind,cfg):
    from plots import palette,wrapped,PALETTE
    r,names=data(p,kind,cfg);traces=[];layout=dict(title=dict(text=cfg['title'] or TITLES[kind]),font=dict(size=cfg['font']),margin=dict(l=155,r=75,t=85,b=125),paper_bgcolor='white',plot_bgcolor='white',showlegend=False)
    metadata={'kind':kind,'rows':json.loads(r.to_json(orient='records')),'comparisons':names}
    marker=dict(size=np.sqrt(cfg['point_size'])*1.5,opacity=cfg['point_alpha'])
    if r.empty:layout['annotations']=[dict(text='无符合条件的结果；请检查比较、方案及输入',xref='paper',yref='paper',x=.5,y=.5,showarrow=False)]
    elif kind=='composition':
        for i,name in enumerate(names):traces.append(dict(type='bar',name=name,x=r['sample'].tolist(),y=r[name].tolist(),marker=dict(color=PALETTE[i%len(PALETTE)]),opacity=cfg['point_alpha']))
        layout.update(barmode='stack',showlegend=cfg['legend'],yaxis=dict(title=dict(text='Conditional RNA contribution'),range=[0,1]))
    elif kind=='deg_upset':
        x=list(range(len(r)));traces.append(dict(type='bar',x=x,y=r['count'].tolist(),marker=dict(color=cfg['down']),opacity=cfg['point_alpha'],customdata=r.genes.tolist(),text=r.intersection.tolist(),hovertemplate='%{text}<br>%{y} genes<extra></extra>'))
        for j,row in enumerate(r.itertuples()):
            present=[i for i,n in enumerate(names) if n in row.members.split(';')]
            traces.append(dict(type='scatter',x=[j]*len(names),y=list(range(len(names))),mode='markers',marker=dict(color=[cfg['down'] if i in present else '#e0e5e2' for i in range(len(names))],**marker),yaxis='y2',hoverinfo='skip'))
            if len(present)>1:traces.append(dict(type='scatter',x=[j,j],y=[min(present),max(present)],mode='lines',line=dict(color=cfg['down'],width=2),yaxis='y2',hoverinfo='skip'))
        layout.update(xaxis=dict(tickvals=x,ticktext=[str(i+1) for i in x]),yaxis=dict(domain=[.38,1],title=dict(text='Exact intersection genes')),yaxis2=dict(domain=[0,.3],tickvals=list(range(len(names))),ticktext=wrapped(comparison_labels(p,names),cfg,web=True),autorange='reversed',showgrid=False),annotations=[dict(text=f"{cfg['comparison_profile']} · {cfg['direction']} · top {len(r)} exact intersections",xref='paper',yref='paper',x=.5,y=1.08,showarrow=False)])
    elif kind=='gene_forest':
        lfc,lo,hi=forest_values(r);color=[cfg['up'] if v>=0 else cfg['down'] for v in lfc]
        traces=[dict(type='scatter',mode='markers',x=lfc.tolist(),y=r.comparison.tolist(),customdata=r.gene.tolist(),marker=dict(color=color,**marker),error_x=dict(type='data',symmetric=False,array=(hi-lfc).tolist(),arrayminus=(lfc-lo).tolist()),text=[f"p={p:.4g}; q={q:.4g}" for p,q in zip(r.get('pvalue',pd.Series(np.nan,index=r.index)),r.padj)],hovertemplate='%{y}<br>raw log2FC=%{x:.4g}<br>%{text}<extra></extra>')]
        layout.update(xaxis=dict(title=dict(text='Raw log2 fold change; Wald 95% interval')),yaxis=dict(type='category',tickvals=r.comparison.tolist(),ticktext=wrapped(comparison_labels(p,r.comparison.tolist()),cfg,web=True),autorange='reversed'),shapes=[dict(type='line',x0=0,x1=0,y0=0,y1=1,yref='paper',line=dict(color='#aaaaaa',dash='dash'))])
    else:
        matrix,q,lo,hi=matrix_data(r,kind,cfg);xlabels,ylabels=axis_labels(p,matrix,kind,cfg)
        trace=dict(type='heatmap',x=matrix.columns.tolist(),y=matrix.index.tolist(),z=matrix.to_numpy().tolist(),zmin=lo,zmax=hi,zauto=False,colorscale=palette(cfg)[1],customdata=[[str(g)]*len(matrix.columns) for g in matrix.index],hovertext=[[f'q={v:.3g}' for v in row] for row in q.to_numpy()],hovertemplate='%{y}<br>%{x}: %{z:.4g}<br>%{hovertext}<extra></extra>',showscale=cfg['legend'],colorbar=dict(title=dict(text='log2FC' if kind=='contrast_effects' else 'NES' if cfg['comparison_metric']=='gsea' else 'Score difference')))
        if cfg['show_values'] and matrix.size<=600:trace.update(text=[[f'{value:.2f}'+(' *' if np.isfinite(padj) and padj<cfg['alpha'] else '') if np.isfinite(value) else 'NA' for value,padj in zip(values,qs)] for values,qs in zip(matrix.to_numpy(),q.to_numpy())],texttemplate='%{text}',textfont=dict(size=min(11,cfg['label_font'])))
        traces=[trace];layout.update(xaxis=dict(type='category',tickvals=matrix.columns.tolist(),ticktext=wrapped(xlabels,cfg,web=True)),yaxis=dict(type='category',tickvals=matrix.index.tolist(),ticktext=wrapped(ylabels,cfg,web=True),autorange='reversed'));metadata['matrix']=matrix.to_dict()
    configure_layout(layout,cfg,reverse_y=kind in ('contrast_effects','pathway_comparison','gene_forest'))
    note='多比较结果来自同一研究；颜色编码效应，星号表示各比较内q<阈值，不是整体FDR。UpSet固定同一方案/方向并只显示前top个精确交集。森林图显示原始Wald区间；旧任务缺少原始SE时不伪造区间。'
    if kind=='composition':note='参考条件下的RNA混合贡献；请同时查看拟合诊断，不是已校准的细胞数量比例。'
    return dict(data=traces,layout=layout,metadata=metadata,category='advanced',settings=cfg,config=dict(responsive=True,displaylogo=False),note=note)


def render(p,kind,cfg,output=None):
    from plots import plt,palette,wrapped,font_manager,PALETTE
    from matplotlib.gridspec import GridSpec
    from matplotlib.text import Text
    output=Path(output) if output else p
    r,names=data(p,kind,cfg);fig=plt.figure(figsize=(cfg['width'],cfg['height']),layout='constrained');reverse=False
    try:
        if r.empty:ax=fig.add_subplot();ax.text(.5,.5,'No eligible results',ha='center',transform=ax.transAxes)
        elif kind=='composition':
            ax=fig.add_subplot();bottom=np.zeros(len(r))
            for i,name in enumerate(names):ax.bar(r['sample'],r[name],bottom=bottom,label=name,color=PALETTE[i%len(PALETTE)],alpha=cfg['point_alpha']);bottom+=r[name].to_numpy()
            ax.set_ylim(0,1);ax.set_ylabel('Conditional RNA contribution');ax.tick_params(axis='x',rotation=45)
            if cfg['legend']:ax.legend(fontsize=cfg['label_font'])
        elif kind=='deg_upset':
            gs=GridSpec(2,1,figure=fig,height_ratios=[3,1.5]);ax=fig.add_subplot(gs[0]);dot=fig.add_subplot(gs[1],sharex=ax);x=np.arange(len(r));ax.bar(x,r['count'],color=cfg['down'],alpha=cfg['point_alpha']);ax.set_ylabel('Exact intersection genes')
            for j,row in enumerate(r.itertuples()):
                present=[i for i,n in enumerate(names) if n in row.members.split(';')];dot.scatter([j]*len(names),range(len(names)),c=[cfg['down'] if i in present else '#e0e5e2' for i in range(len(names))],s=cfg['point_size'],alpha=cfg['point_alpha'])
                if len(present)>1:dot.plot([j,j],[min(present),max(present)],color=cfg['down'])
            dot.set_yticks(range(len(names)),wrapped(comparison_labels(p,names),cfg));dot.invert_yaxis();dot.set_xticks(x,[str(i+1) for i in x]);dot.set_xlabel(cfg['comparison_profile']+' / '+cfg['direction']);ax.tick_params(labelbottom=False)
        elif kind=='gene_forest':
            ax=fig.add_subplot();lfc,lo,hi=forest_values(r);y=np.arange(len(r));color=[cfg['up'] if v>=0 else cfg['down'] for v in lfc];ax.scatter(lfc,y,c=color,s=cfg['point_size'],alpha=cfg['point_alpha']);ax.axvline(0,color='#aaa',ls='--')
            if 'ci_low_raw' in r:ax.errorbar(lfc,y,xerr=np.stack([lfc-lo,hi-lfc]),fmt='none',ecolor=cfg['neutral'])
            ax.set_yticks(y,wrapped(comparison_labels(p,r.comparison.tolist()),cfg));ax.invert_yaxis();reverse=True;ax.set_xlabel('Raw log2 fold change; Wald 95% interval')
        else:
            ax=fig.add_subplot();x,q,lo,hi=matrix_data(r,kind,cfg);xlabels,ylabels=axis_labels(p,x,kind,cfg);cmap=palette(cfg)[0].copy();cmap.set_bad('white');im=ax.imshow(np.ma.masked_invalid(x.to_numpy()),cmap=cmap,vmin=lo,vmax=hi,aspect='auto')
            if cfg['legend']:fig.colorbar(im,ax=ax,label='log2FC' if kind=='contrast_effects' else 'NES' if cfg['comparison_metric']=='gsea' else 'Score difference')
            ax.set_xticks(range(len(x.columns)),wrapped(xlabels,cfg),rotation=35,ha='right');ax.set_yticks(range(len(x)),wrapped(ylabels,cfg));reverse=True
            if cfg['show_values'] and x.size<=600:
                for i in range(len(x)):
                    for j in range(len(x.columns)):
                        value=x.iloc[i,j];padj=q.iloc[i,j];text=f'{value:.2f}'+(' *' if np.isfinite(padj) and padj<cfg['alpha'] else '') if np.isfinite(value) else 'NA';ax.text(j,i,text,ha='center',va='center',fontsize=min(11,cfg['label_font']))
        for axis in fig.axes:
            axis.tick_params(labelsize=cfg['label_font']);axis.xaxis.label.set_size(cfg['font']);axis.yaxis.label.set_size(cfg['font'])
        if kind in ('gene_forest','deg_upset','composition') and cfg['grid']:ax.grid(True,color='#e8ece9',axis='y');ax.set_axisbelow(True)
        for name in ('x','y'):
            if cfg[name+'label']:getattr(ax,'set_'+name+'label')(cfg[name+'label'])
            low,high=cfg[name+'_min'],cfg[name+'_max']
            if low!='' or high!='':
                current=getattr(ax,'get_'+name+'lim')();a=float(low) if low!='' else min(current);b=float(high) if high!='' else max(current)
                getattr(ax,'set_'+name+'lim')(*((b,a) if name=='y' and reverse else (a,b)))
        fig.suptitle(cfg['title'] or TITLES[kind],fontsize=cfg['font']+3)
        available={f.name for f in font_manager.fontManager.ttflist};family=next((name for name in ('Noto Sans CJK SC','Noto Sans CJK JP','Microsoft YaHei','WenQuanYi Zen Hei') if name in available),'DejaVu Sans')
        for text in fig.findobj(Text):text.set_fontfamily(family)
        with plt.rc_context({'svg.fonttype':'none','pdf.fonttype':42}):
            for fmt in ('svg','pdf','png'):fig.savefig(output/(kind+'.'+fmt),dpi=cfg['dpi'])
        r.to_csv(output/(kind+'.plot-data.csv'),index=False);(output/(kind+'.settings.json')).write_text(json.dumps(cfg,ensure_ascii=False,indent=2))
    finally:plt.close(fig)
