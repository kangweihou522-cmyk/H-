"""Explicit estimands and sample-design validation shared by all entry points."""
import itertools,re
import numpy as np
import pandas as pd


def prepare_design(cfg,m,inferential=True):
    kind=cfg.get('design_type','independent')
    if kind not in ('independent','paired','interaction'):raise ValueError('实验设计必须为独立、配对或交互')
    cov=cfg.get('covariates')
    if cov is None:cov=[dict(name='batch',type='categorical')] if 'batch' in m and m.batch.nunique()>1 else []
    if not isinstance(cov,list):raise ValueError('协变量必须为列表')
    cov=[dict(c) for c in cov]
    subject=cfg.get('subject_column','subject')
    repeated=subject in m and m[subject].nunique()<len(m)
    if inferential and kind=='paired':
        if not repeated or (m[subject]=='').any():raise ValueError('配对设计需要完整且重复出现的对象ID')
        if m[[subject,'group']].duplicated().any():raise ValueError('同一对象在同组有多次观测；请合并技术重复或使用专门重复测量模型')
        if (m.groupby(subject).group.nunique()<2).any():raise ValueError('配对对象须在至少两个组中有观测；请核对配对或另用专门模型')
        if not any(c.get('name')==subject for c in cov):cov.append(dict(name=subject,type='categorical'))
    if inferential and repeated and not any(c.get('name')==subject and c.get('type','categorical')=='categorical' for c in cov):
        raise ValueError('检测到重复对象：请选择配对设计或明确加入subject分类协变量；不能按独立动物静默运行')
    cfg['design_type']=kind;cfg['subject_column']=subject;cfg['covariates']=cov
    cfg['omnibus']=cfg.get('omnibus',True)
    if not isinstance(cfg['omnibus'],bool):raise ValueError('总体检验选项必须为布尔值')
    cfg['trait_adjustment']=cfg.get('trait_adjustment','design')
    if cfg['trait_adjustment'] not in ('design','covariates'):raise ValueError('表型关联调整范围无效')
    cfg['ora_background']=cfg.get('ora_background','pvalue')
    if cfg['ora_background'] not in ('pvalue','stat'):raise ValueError('ORA背景策略无效')
    cfg['effect_test']=cfg.get('effect_test',False)
    if not isinstance(cfg['effect_test'],bool):raise ValueError('效应阈值检验选项无效')
    cfg['global_fdr']=cfg.get('global_fdr',False)
    if not isinstance(cfg['global_fdr'],bool):raise ValueError('跨比较BH选项无效')
    cfg['wgcna_correlation']=cfg.get('wgcna_correlation','pearson')
    if cfg['wgcna_correlation'] not in ('pearson','bicor'):raise ValueError('网络相关方法无效')
    cfg['wgcna_resamples']=int(cfg.get('wgcna_resamples',0))
    if not 0<=cfg['wgcna_resamples']<=50:raise ValueError('网络稳定性重抽样次数须为0–50')
    if kind=='interaction' and inferential:
        factor=cfg.get('interaction_factor','')
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',factor) or factor in ('group','sample') or factor not in m:raise ValueError('请选择有效的第二分类因素')
        if factor in [c.get('name') for c in cov]:raise ValueError('交互因素不能同时重复填入协变量')
        if (m[factor]=='').any() or m[factor].nunique()<2:raise ValueError('交互因素需要至少两个完整水平')
        cells=pd.crosstab(m.group,m[factor]);
        if (cells<2).any().any():raise ValueError('交互设计每个组×因素组合至少需要2个独立重复')
        cfg['interaction_reference']=cfg.get('interaction_reference') or sorted(m[factor].unique())[0]
        if cfg['interaction_reference'] not in set(m[factor]):raise ValueError('交互参考水平不存在')
        if cfg.get('shrinkage','none')!='none':raise ValueError('交互数字对比当前请使用原始LFC，不能套用normal收缩')
    else:cfg['interaction_factor']=''
    cfg['independent_subjects']=int(m[subject].nunique()) if subject in m else len(m)
    return cfg


def design_matrix(cfg,m):
    parts=[np.ones((len(m),1))];names=['Intercept']
    for c in cfg['covariates']:
        name=c['name']
        if c.get('type','categorical')=='numeric':parts.append(pd.to_numeric(m[name]).to_numpy(float)[:,None]);names.append(name)
        else:
            x=pd.get_dummies(m[name],drop_first=True,dtype=float);parts.append(x.to_numpy());names.extend(name+'='+str(x) for x in x.columns)
    group_levels=sorted(m.group.unique());reference=cfg.get('control',group_levels[0]);group_levels=[reference]+[x for x in group_levels if x!=reference]
    group=pd.get_dummies(pd.Categorical(m.group,categories=group_levels),drop_first=True,dtype=float);parts.append(group.to_numpy());names.extend('group='+str(x) for x in group.columns)
    if cfg.get('interaction_factor'):
        factor=cfg['interaction_factor'];levels=sorted(m[factor].unique());ref=cfg['interaction_reference'];levels=[ref]+[v for v in levels if v!=ref];x=pd.get_dummies(pd.Categorical(m[factor],categories=levels),drop_first=True,dtype=float);parts.append(x.to_numpy());names.extend(factor+'='+str(v) for v in x.columns)
        for g in group:
            for f in x:parts.append((group[g]*x[f]).to_numpy()[:,None]);names.append('group='+str(g)+':'+factor+'='+str(f))
    matrix=np.concatenate(parts,axis=1)
    return matrix,names


def generate_comparisons(groups,mode='control',control=None,factor_levels=None,reference=None):
    groups=list(dict.fromkeys(groups))
    if mode not in ('control','all'):raise ValueError('比较生成方式无效')
    if len(groups)<2:raise ValueError('至少两个组')
    if mode=='control' and control not in groups:raise ValueError('请选择样本中的统一对照')
    pairs=[(control,g) for g in groups if g!=control] if mode=='control' else list(itertools.combinations(groups,2))
    # Names do not depend on user text being valid as an R or file identifier.
    out=[]
    for i,(a,b) in enumerate(pairs,1):
        for j,level in enumerate(factor_levels or [None],1):
            out.append(dict(name=f'cmp{i}'+(f'_level{j}' if level is not None else ''),control=a,treatment=b,kind='pairwise',**({'factor_level':level} if level is not None else {})))
        for j,level in enumerate(factor_levels or [],1):
            if level!=reference:out.append(dict(name=f'cmp{i}_interaction{j}',control=a,treatment=b,kind='interaction',factor_level=level,reference_level=reference))
    if len(out)>30:raise ValueError(f'生成{len(out)}个比较，超过每任务30个；请选统一对照或减少组合')
    return out
