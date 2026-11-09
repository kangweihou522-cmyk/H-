"""Cross-comparison display context and recoverable figure export regressions."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from comparison_context import expression_samples
from plots import interactive_data,render,get_settings,chart_catalog
from results import gene_data

@pytest.fixture
def study(tmp_path):
    rows=[dict(sample=f'{g}{s}{i}',group=g,sex=s,mouse=f'{s}{i}') for g in ['C','T','R'] for s in ['F','M'] for i in range(2)]
    m=pd.DataFrame(rows);m.to_csv(tmp_path/'samples.csv',index=False)
    values=pd.DataFrame(np.arange(36).reshape(3,12)+1,index=['g1','g2','g3'],columns=m['sample'])
    for file in ['vst.csv','normalized.csv']:values.to_csv(tmp_path/file)
    comparisons=[dict(name='female',control='C',treatment='T',factor_level='F',kind='pairwise'),dict(name='male',control='C',treatment='T',factor_level='M',kind='pairwise'),dict(name='interaction',control='C',treatment='T',factor_level='M',reference_level='F',kind='interaction')]
    cfg=dict(comparisons=comparisons,interaction_factor='sex',interaction_reference='F',subject_column='mouse',mode='counts')
    (tmp_path/'config.json').write_text(json.dumps(cfg))
    for i,c in enumerate(comparisons):
        q=tmp_path/'contrasts'/c['name'];q.mkdir(parents=True)
        pd.DataFrame(dict(gene=['g1','g2','g3'],symbol=['G1','G2','G3'],log2FoldChange=[2+i,-2-i,.1],log2FoldChange_raw=[3+i,-3-i,.1],ci_low_raw=[2+i,-4-i,-.4],ci_high_raw=[4+i,-2-i,.6],pvalue=[.0001,.0002,.8],padj=[.001,.002,.8],baseMean=[10,20,30],significant=[True,True,False])).to_csv(q/'differential.csv',index=False)
        pd.read_csv(q/'differential.csv').to_csv(q/'deg-all.csv',index=False)
    return tmp_path,cfg,m


def test_gene_card_retains_full_data_but_focuses_conditional_samples(study):
    p,cfg,m=study
    d=gene_data(p,'G1','female')
    assert len(d['samples'])==12 and len(d['expression_samples'])==4
    assert {s['sex'] for s in d['expression_samples']}=={'F'}
    assert {s['group'] for s in d['expression_samples']}=={'F · C','F · T'}
    assert d['comparison']['name']=='female' and d['candidates'][0]['log2FoldChange']==2
    assert len(gene_data(p,'g1','interaction')['expression_samples'])==8
    assert len(gene_data(p,'g1','male')['expression_samples'])==4
    assert all(a['pvalue'] is None for a in d['associations'])


def test_expression_scope_and_pair_lines_agree_with_metadata(study):
    p,cfg,m=study
    payload=interactive_data(p,'gene_expression',dict(contrast='female'),dict(gene='g1',gene_scale='vst',connect_pairs=True))
    assert len(payload['metadata']['samples'])==4
    assert len([t for t in payload['data'] if t.get('mode')=='lines'])==2
    payload=interactive_data(p,'gene_expression',dict(contrast='interaction'),dict(gene='g1',gene_scale='vst',connect_pairs=True))
    lines=[t for t in payload['data'] if t.get('mode')=='lines']
    assert len(lines)==4 and all(len(t['x'])==2 for t in lines)
    payload=interactive_data(p,'gene_expression',dict(contrast='female'),dict(gene='g1',gene_scope='all'))
    assert len(payload['metadata']['samples'])==12 and '全部样本' in payload['metadata']['expression_context']['label']
    assert {r['biological_group'] for r in payload['metadata']['samples']}=={'C','T','R'}
    render(p/'contrasts'/'female','gene_expression',dict(gene='g1',gene_scope='all',dpi=72,width=5,height=4))
    saved=json.loads((p/'contrasts'/'female'/'gene_expression.plot-data.json').read_text())
    assert saved['samples']==payload['metadata']['samples']


def test_missing_factor_fails_clearly_with_explicit_all_samples_fallback(study):
    p,cfg,m=study;m.drop(columns='sex').to_csv(p/'samples.csv',index=False)
    with pytest.raises(ValueError,match='因素列'):gene_data(p,'g1','female')
    payload=interactive_data(p,'gene_expression',dict(contrast='female'),dict(gene='g1',gene_scope='all'))
    assert len(payload['metadata']['samples'])==12
    # Non-inferential handoffs retain the original stratification column.
    rows,context=expression_samples(m.to_dict('records'),dict(cfg,interaction_factor='',display_factor='sex'),'male')
    assert len(rows)==4 and context['factor']=='sex'


def test_advanced_palettes_limits_order_and_raw_wald_intervals(study):
    p,cfg,m=study
    options=dict(genes='g2,g1',show_values=True,color_map='purple_green',color_min='-5',color_max='5',xlabel='Contrasts',grid=False)
    payload=interactive_data(p,'contrast_effects',{},options);trace=payload['data'][0]
    assert trace['y']==['g2','g1'] and trace['zmin']==-5 and trace['zmax']==5 and trace['text'][0][0].endswith('*')
    assert 'q=0.002'==trace['hovertext'][0][0] and '%{hovertext}' in trace['hovertemplate']
    assert payload['layout']['xaxis']['title']['text']=='Contrasts' and not payload['layout']['xaxis']['showgrid']
    assert 'sex=F' in payload['layout']['xaxis']['ticktext'][0]
    render(p,'contrast_effects',dict(options,dpi=72))
    exported=pd.read_csv(p/'contrast_effects.plot-data.csv');assert set(exported.gene)=={'g1','g2'}
    forest=interactive_data(p,'gene_forest',{},dict(gene='g1'))
    assert forest['data'][0]['x']==[3,4,5] and forest['data'][0]['error_x']['array']==[1,1,1]
    assert forest['layout']['yaxis']['autorange']=='reversed'


def test_activity_only_and_footprint_comparisons_are_discoverable(study):
    p,cfg,m=study
    for c in cfg['comparisons']:
        pd.DataFrame(dict(term=['P1','P2'],logFC=[1,-1],t=[3,-3],**{'P.Value':[.01,.01],'adj.P.Val':[.02,.02]})).to_csv(p/'contrasts'/c['name']/'footprint-activity-differential.csv',index=False)
    assert any(c['kind']=='pathway_comparison' for c in chart_catalog(p))
    settings=get_settings(p,'pathway_comparison');assert settings['comparison_metric']=='footprint'
    payload=interactive_data(p,'pathway_comparison',{});assert payload['data'][0]['z']==[[1,1,1],[-1,-1,-1]]
    render(p,'pathway_comparison',dict(dpi=72))
    assert (p/'pathway_comparison.png').stat().st_size>1000


@pytest.mark.parametrize('failure',['png','publish'])
def test_failed_plot_export_preserves_complete_last_generation(study,monkeypatch,failure):
    import plots
    from matplotlib.figure import Figure
    p,cfg,m=study;target=p/'contrasts'/'female';options=dict(gene='g1',dpi=72,width=4,height=4)
    render(target,'gene_expression',options)
    files=list(target.glob('gene_expression.*'));before={f.name:f.read_bytes() for f in files}
    if failure=='png':
        original=Figure.savefig
        def save(self,path,*args,**kwargs):
            if Path(path).suffix=='.png':raise OSError('Synthetic PNG failure')
            return original(self,path,*args,**kwargs)
        monkeypatch.setattr(Figure,'savefig',save)
    else:
        original=plots.os.replace;failed=False
        def replace(source,target):
            nonlocal failed
            if Path(target).name=='gene_expression.png' and not failed:failed=True;raise OSError('Synthetic publish failure')
            return original(source,target)
        monkeypatch.setattr(plots.os,'replace',replace)
    with pytest.raises(OSError,match='Synthetic'):render(target,'gene_expression',dict(options,title='Changed'))
    assert {f.name:f.read_bytes() for f in files}==before
    assert not list(target.glob('.plots-*'))


def test_expression_handoff_keeps_factor_without_claiming_an_interaction_fit(study):
    from validation import validate
    p,cfg,m=study
    pd.DataFrame(np.arange(144).reshape(12,12),index=[f'g{i}' for i in range(12)],columns=m['sample']).to_csv(p/'vst.csv')
    valid,_=validate(dict(cfg,analysis_stage='expression',samples=str(p/'samples.csv'),expression=str(p/'vst.csv'),expression_unit='vst',design_type='interaction',covariates=[]))
    assert valid['display_factor']=='sex' and valid['interaction_factor']==''
    assert valid['analysis_stage']=='expression' and valid['comparisons']==[]



def test_imported_interaction_contrast_can_run_standalone_enrichment(study):
    from validation import validate
    p,cfg,m=study
    cmp=cfg['comparisons'][2]
    pd.DataFrame(dict(gene=[f'g{i}' for i in range(12)],log2FoldChange=np.linspace(-2,2,12),pvalue=np.linspace(.01,.9,12),padj=np.ones(12))).to_csv(p/'contrasts'/cmp['name']/'differential.csv',index=False)
    valid,_=validate(dict(cfg,analysis_stage='enrichment',differential_input=str(p/'contrasts'/cmp['name']/'differential.csv'),comparisons=[cmp],samples='',expression='',design_type='interaction'))
    assert valid['comparisons'][0]['kind']=='interaction' and valid['comparisons'][0]['reference_level']=='F'
    assert valid['analysis_stage']=='enrichment' and valid['run_gsea'] is False
    assert valid['display_factor']=='sex' and valid['interaction_factor']==''
