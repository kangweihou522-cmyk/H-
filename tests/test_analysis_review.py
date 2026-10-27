"""Regression checks for reviewed statistics, references and shared plot controls."""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from discovery import adjusted_traits, bh
from scripts.prepare_decoy import prepare
from scripts.make_mining_data import create
from validation import validate


def test_partial_pearson_and_global_bh(tmp_path):
    rng=np.random.default_rng(12);n=20
    group=np.repeat([0,1],10);trait=group*5+rng.normal(size=n);batch=np.tile([0,1],10)
    design=np.stack([np.ones(n),group,batch],axis=1)
    expression=np.stack([5*group+rng.normal(size=n),2*trait+rng.normal(0,.3,size=n)])
    samples=pd.DataFrame(dict(sample=[f'S{i}' for i in range(n)],group=group,score=trait,confounded=group*3+batch))
    samples.to_csv(tmp_path/'samples.csv',index=False)
    (tmp_path/'config.json').write_text(json.dumps(dict(traits=['score','confounded'])))
    pd.DataFrame(design,index=samples['sample']).to_csv(tmp_path/'design-matrix.csv')
    pd.DataFrame(expression,index=['group_gene','trait_gene'],columns=samples['sample']).to_csv(tmp_path/'vst.csv')
    status=adjusted_traits(tmp_path);assert status['state']=='completed'
    assert '完全解释' in status['skipped'][0]['reason']
    result=pd.read_csv(tmp_path/'gene-trait-adjusted.csv')
    residual=trait-design@np.linalg.lstsq(design,trait,rcond=None)[0]
    for row,x in zip(result.itertuples(),expression):
        rx=x-design@np.linalg.lstsq(design,x,rcond=None)[0]
        r=stats.pearsonr(rx,residual).statistic
        t=abs(r)*np.sqrt(16/(1-r*r))
        assert row.df==16 and np.isclose(row.correlation,r)
        assert np.isclose(row.pvalue,2*stats.t.sf(t,16))
    assert np.allclose(result.padj,bh(result.pvalue))
    assert result.iloc[1].correlation>.95


def test_bh_and_missing_tests():
    values=[.04,np.nan,.001,.8,.02,.2]
    assert np.allclose(bh(values)[[0,2,3,4,5]],[.0666666667,.005,.8,.05,.25])
    assert np.isnan(bh(values)[1])


def test_reference_decoys_stream_gzip_and_reject_collision(tmp_path):
    tx=tmp_path/'tx.fa';tx.write_text('>tx1 description\nACGT\n>tx2\nAACC')
    genome=tmp_path/'genome.fa.gz'
    with gzip.open(genome,'wt') as f:f.write('>chr1 version\nAACCGG\n>chr2\nTTGG\n')
    combined,decoys=prepare(tx,genome,tmp_path/'reference')
    assert decoys.read_text()=='chr1\nchr2\n'
    assert combined.read_text().count('>')==4 and '\n>chr1' in combined.read_text()
    tx.write_text('>chr1\nTTGG\n')
    with pytest.raises(ValueError,match='冲突'):prepare(tx,genome,tmp_path/'bad')
    assert not (tmp_path/'bad/gentrome.fa').exists()


def test_single_level_covariate_and_library_type(tmp_path):
    cfg=create(tmp_path);m=pd.read_csv(tmp_path/'samples.csv');m['constant']='B1';m.to_csv(tmp_path/'samples.csv',index=False)
    cfg['covariates']=[dict(name='constant',type='categorical')]
    with pytest.raises(ValueError,match='一个水平'):validate(cfg)


def test_plot_controls_share_matrix_and_dimensions(tmp_path):
    from plots import get_settings,interactive_data,load_plot
    cfg=create(tmp_path);counts=pd.read_csv(tmp_path/'counts.csv',index_col=0);v=np.log2(counts+1);v.to_csv(tmp_path/'vst.csv')
    m=pd.read_csv(tmp_path/'samples.csv');v.corr().to_csv(tmp_path/'sample-correlation.csv')
    settings=get_settings(tmp_path,'volcano');assert settings['width']/settings['height']<1.1
    opts=dict(cluster=False,heatmap_scale='none',show_values=True,color_min='-3',color_max='3')
    v.iloc[:3].to_csv(tmp_path/'pathway-activity.csv')
    payload=interactive_data(tmp_path,'activity',{},opts)
    assert payload['data'][0]['zmin']==-3 and payload['data'][0]['zmax']==3
    assert np.allclose(payload['data'][0]['z'],load_plot(tmp_path,'activity',get_settings(tmp_path,'activity',opts))[1].to_numpy())
    result=interactive_data(tmp_path,'correlation',{},dict(cluster=False,correlation_method='spearman'))
    assert np.allclose(result['data'][0]['z'],v.corr(method='spearman'))
    with pytest.raises(ValueError,match='范围'):get_settings(tmp_path,'volcano',dict(point_size=float('nan')))


def test_preview_and_lazy_export_leave_statistics_intact(tmp_path,monkeypatch):
    import app as station
    job='a'*32;p=tmp_path/job;p.mkdir()
    samples=pd.DataFrame(dict(sample=['C1','C2','T1','T2'],group=['C','C','T','T']))
    samples.to_csv(p/'samples.csv',index=False)
    (p/'status.json').write_text(json.dumps(dict(state='completed')))
    result=pd.DataFrame(dict(gene=['g1','g2'],symbol=['A','B'],log2FoldChange=[2,-2],pvalue=[.001,.01],padj=[.01,.04],baseMean=[50,60]))
    result.to_csv(p/'differential.csv',index=False)
    pd.DataFrame([[1,2,4,5],[4,3,2,1]],index=['g1','g2'],columns=samples['sample']).to_csv(p/'vst.csv')
    before=(p/'differential.csv').read_bytes()
    monkeypatch.setattr(station,'RUNS',tmp_path)
    client=station.app.test_client()
    response=client.post('/api/jobs/'+job+'/preview',json=dict(kind='volcano',width=7,height=6.5))
    assert response.status_code==200 and not (p/'volcano.settings.json').exists()
    response=client.get('/api/jobs/'+job+'/chart-download/volcano/svg')
    assert response.status_code==200 and b'<svg' in response.data
    assert json.loads((p/'volcano.settings.json').read_text())['style_version']==3
    assert before==(p/'differential.csv').read_bytes()
    assert client.get('/api/jobs/'+job+'/chart-download/volcano/not-an-extension').status_code==400


def test_numeric_identifiers_keep_leading_zeros(tmp_path):
    from dataio import read_frame
    from results import gene_data
    from plots import load_plot,get_settings
    samples=pd.DataFrame(dict(sample=['001','002','003','004'],group=['01','01','02','02']))
    samples.to_csv(tmp_path/'samples.csv',index=False)
    values=pd.DataFrame([[10,11,30,40],[40,30,20,10]],index=['0001','0002'],columns=samples['sample'])
    values.to_csv(tmp_path/'normalized.csv');np.log2(values+1).to_csv(tmp_path/'vst.csv')
    pd.DataFrame(dict(gene=['0001','0002'],symbol=['A','B'],padj=[.01,.02],baseMean=[20,25],log2FoldChange=[2,-2])).to_csv(tmp_path/'differential.csv',index=False)
    assert read_frame(tmp_path/'normalized.csv',index_col=0).index.tolist()==['0001','0002']
    result=gene_data(tmp_path,'0001');assert result['gene']=='0001' and result['samples'][0]['sample']=='001'
    data=load_plot(tmp_path,'heatmap',get_settings(tmp_path,'heatmap',dict(genes='0001')))[1]
    assert data.index.tolist()==['0001'] and set(data.columns)==set(samples['sample'])


def test_quant_import_rejects_reused_samples_and_invalid_metadata(tmp_path):
    cfg=create(tmp_path)
    names=['tx1.1|gene1','tx2.1|gene2']
    mapping=tmp_path/'tx2gene.tsv';pd.DataFrame(dict(tx=['tx1','tx2'],gene=['g1','g2'])).to_csv(mapping,sep='\t',index=False)
    samples=pd.read_csv(tmp_path/'samples.csv',dtype=str)
    paths=[]
    for i in range(len(samples)):
        f=tmp_path/f'quant{i}.sf'
        pd.DataFrame(dict(Name=names,Length=[1000,1200],EffectiveLength=[800,1000],TPM=[4.,6.],NumReads=[10.,20.])).to_csv(f,sep='\t',index=False)
        paths.append(str(f))
    samples['quant']=paths;samples.to_csv(tmp_path/'samples.csv',index=False)
    cfg.update(mode='salmon',tx2gene=str(mapping))
    validate(cfg)
    samples.loc[1,'quant']=samples.loc[0,'quant'];samples.to_csv(tmp_path/'samples.csv',index=False)
    with pytest.raises(ValueError,match='重复当作'):validate(cfg)
    samples.loc[1,'quant']=paths[1];samples.to_csv(tmp_path/'samples.csv',index=False)
    q=pd.read_csv(paths[0],sep='\t');q.loc[0,'EffectiveLength']=-1;q.to_csv(paths[0],sep='\t',index=False)
    with pytest.raises(ValueError,match='定量无效'):validate(cfg)


def test_salmon_single_end_uses_explicit_fragment_distribution():
    from app import salmon_quant_args
    cfg=dict(threads=2,library='A',fragment_mean=250,fragment_sd=40)
    single=salmon_quant_args(cfg,'index','out',['R1.fastq.gz'])
    paired=salmon_quant_args(cfg,'index','out',['R1.fastq.gz','R2.fastq.gz'])
    assert single[-6:]==['-r','R1.fastq.gz','--fldMean','250','--fldSD','40']
    assert paired[-4:]==['-1','R1.fastq.gz','-2','R2.fastq.gz'] and '--fldMean' not in paired


def test_r_module_fallback_is_consistent_across_plots():
    from plots import module_colors
    assert module_colors(['lightcyan1','lightcyan1'],'#cccccc')[0]==module_colors(['lightcyan1'],'#cccccc')[0]
    assert module_colors(['lightcyan1','lightcyan1'],'#cccccc')[0]==module_colors(['lightcyan1','lightcyan1'],'#cccccc')[1]
