"""Keep displayed samples, labels and clicked statistics tied to an explicit contrast."""
import pandas as pd


def resolve(cfg,name=None):
    comparisons=cfg.get('comparisons') or []
    if not comparisons and cfg.get('control') and cfg.get('treatment'):comparisons=[dict(name='legacy',control=cfg['control'],treatment=cfg['treatment'])]
    if name:
        return next((dict(c) for c in comparisons if c.get('name')==name),None)
    return dict(comparisons[0]) if comparisons else None


def label(c,factor=''):
    if not c:return '全部样本 · 未指定比较'
    text=f"{c.get('treatment','')} / {c.get('control','')}"
    if c.get('factor_level') is not None:
        text+=' · '+(f"差中之差 {c['factor_level']} − {c.get('reference_level','')}" if c.get('kind')=='interaction' else f"{factor or '因素'}={c['factor_level']}")
    return text+' · '+str(c.get('name',''))


def expression_samples(samples,cfg,contrast=None,scope='comparison'):
    frame=pd.DataFrame(samples).copy();cmp=resolve(cfg,contrast);factor=cfg.get('interaction_factor') or cfg.get('display_factor','')
    if scope not in ('comparison','all'):raise ValueError('表达图样本范围无效')
    if frame.empty:return [],dict(label=label(cmp,factor),scope=scope,samples=0)
    if scope=='comparison' and cmp:
        if factor and cmp.get('factor_level') is not None and factor not in frame:raise ValueError('所选条件比较缺少因素列 '+factor+'；请补充样本信息，或显式选择全部样本展示')
        frame=frame[frame.group.isin([cmp['control'],cmp['treatment']])].copy()
        if factor and factor in frame and cmp.get('factor_level') is not None:
            levels=[cmp['factor_level']]
            if cmp.get('kind')=='interaction':levels=[cmp['reference_level'],cmp['factor_level']]
            frame=frame[frame[factor].isin(levels)].copy()
    frame['biological_group']=frame.group
    groups=([cmp['control'],cmp['treatment']] if cmp else [])+sorted(set(frame.group)-set([cmp['control'],cmp['treatment']] if cmp else []))
    if factor and factor in frame:
        reference=cfg.get('interaction_reference');levels=sorted(frame[factor].unique());levels=([reference] if reference in levels else [])+[x for x in levels if x!=reference]
        frame['_level_order']=frame[factor].map(dict(zip(levels,range(len(levels)))))
        frame['_group_order']=frame.group.map(dict(zip(groups,range(len(groups)))))
        frame=frame.sort_values(['_level_order','_group_order'],kind='stable').drop(columns=['_level_order','_group_order'])
        frame['group']=frame[factor].astype(str)+' · '+frame.group.astype(str)
        frame['stratum']=frame[factor].astype(str)
    else:
        frame['_group_order']=frame.group.map(dict(zip(groups,range(len(groups)))))
        frame=frame.sort_values('_group_order',kind='stable').drop(columns='_group_order')
    return frame.to_dict('records'),dict(label=('全部样本 · ' if scope=='all' else '')+label(cmp,factor),scope=scope,samples=len(frame),factor=factor,kind=(cmp or {}).get('kind','pairwise'),note='未调整的逐样本表达；效应/P来自所选模型，不能由图上均值比反推。')
