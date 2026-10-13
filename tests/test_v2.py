from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from validation import validate, local_path
from results import curve_data, gene_data, table_page
from scripts.fixture import generate

@pytest.fixture
def data(tmp_path):
    generate(tmp_path)
    cfg=dict(mode='counts',samples=str(tmp_path/'samples.csv'),counts=str(tmp_path/'counts.csv'),control='Control',treatment='Treatment')
    return cfg,tmp_path

def test_windows_path_conversion():
    assert str(local_path(r'D:\data\samples.csv'))=='/mnt/d/data/samples.csv'
    assert str(local_path('file:///D:/data/samples.csv'))=='/mnt/d/data/samples.csv'

def test_multiple_group_contrasts_and_numeric_covariate(data):
    cfg,p=data
    m=pd.read_csv(p/'samples.csv');m.loc[len(m)]=['R1','Rescue'];m.loc[len(m)]=['R2','Rescue'];m['age']=[10,12,11,12,10,13,12,10];m.to_csv(p/'samples.csv',index=False)
    a=pd.read_csv(p/'counts.csv',index_col=0);a['R1']=a.C1;a['R2']=a.C2;a.to_csv(p/'counts.csv')
    cfg.update(comparisons=[dict(name='TC',control='Control',treatment='Treatment'),dict(name='RT',control='Treatment',treatment='Rescue')],covariates=[dict(name='age',type='numeric')])
    result,m=validate(cfg)
    assert len(result['comparisons'])==2 and len(m)==8

def test_duplicate_contrast_and_missing_trait(data):
    cfg,_=data
    c=dict(name='same',control='Control',treatment='Treatment');cfg['comparisons']=[c.copy(),c.copy()]
    with pytest.raises(ValueError,match='唯一'):validate(cfg)
    cfg.pop('comparisons');cfg['traits']=['missing']
    with pytest.raises(ValueError,match='表型'):validate(cfg)

def test_tabular_xls_counts(data):
    cfg,p=data
    a=pd.read_csv(p/'counts.csv',index_col=0);a.to_csv(p/'gene_count.xls',sep='\t');cfg['counts']=str(p/'gene_count.xls')
    _,m=validate(cfg);assert len(m)==6

def test_exclusion_rechecks_replicates(data):
    cfg,_=data;cfg['excluded_samples']=['C1','C2']
    with pytest.raises(ValueError,match='重复'):validate(cfg)

def test_running_sum_and_leading_edge(tmp_path):
    (tmp_path/'gene-sets.json').write_text(json.dumps({'go':{'GO:0006955':['g1','g2']}}))
    pd.DataFrame({'gene':['g1','g2','g3','g4'],'stat':[4,2,-2,-4]}).to_csv(tmp_path/'gsea-ranks.csv',index=False)
    pd.DataFrame({'term':['GO:0006955'],'NES':[2.0],'padj':[.01],'leadingEdge':['g1;g2']}).to_csv(tmp_path/'gsea.csv',index=False)
    curve=curve_data(tmp_path,'go','GO:0006955')
    assert np.allclose(curve['running'],[2/3,1,.5,0])
    assert curve['leading']==['g1','g2'] and curve['hits']==[0,1]

def test_results_table_null_and_search(tmp_path):
    pd.DataFrame({'gene':['A','B'],'symbol':['Hello',None],'padj':[1e-25,np.nan]}).to_csv(tmp_path/'differential.csv',index=False)
    d=table_page(tmp_path,'differential',{'query':'hello'});assert d['total']==1 and d['rows'][0]['padj']==pytest.approx(1e-25)
    d=table_page(tmp_path,'differential',{});assert d['rows'][1]['padj'] is None
    with pytest.raises(ValueError):table_page(tmp_path,'../../config',{})

def test_gene_alias_with_missing_entrez(data):
    _,p=data
    counts=pd.read_csv(p/'counts.csv',index_col=0)
    counts.to_csv(p/'normalized.csv');np.log2(counts+1).to_csv(p/'vst.csv')
    genes=list(counts.index)
    pd.DataFrame({'gene':genes,'symbol':[f'Synthetic_{i}' for i in range(len(genes))],
                  'padj':[np.nan]*len(genes),'baseMean':counts.mean(axis=1).to_numpy()}).to_csv(p/'differential.csv',index=False)
    pd.DataFrame({'gene':genes[:2],'id':genes[:2],'entrez':['123',None],'symbol':['Synthetic_0','Synthetic_1']}).to_csv(p/'gene-annotation.csv',index=False)
    assert gene_data(p,'123')['gene']==genes[0]
    assert gene_data(p,'synthetic_0')['gene']==genes[0]
    from plots import load_plot,get_settings
    chart=load_plot(p,'gene_expression',get_settings(p,'gene_expression'))
    assert len(chart[1]['samples'])==6

def test_wgcna_r_module_color_exports(data):
    _,p=data
    counts=pd.read_csv(p/'counts.csv',index_col=0);np.log2(counts+1).to_csv(p/'vst.csv')
    pd.DataFrame({'gene':counts.index[:3],'module':['grey','lightcyan1','brown'],'kME':[np.nan,.8,.9]}).to_csv(p/'wgcna-modules.csv',index=False)
    from plots import render,interactive_data
    render(p,'wgcna_modules',{})
    assert (p/'wgcna_modules.svg').stat().st_size>100
    assert len(interactive_data(p,'wgcna_modules',{})['data'][0]['marker']['color'])==3
