import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scripts.fixture import generate
from validation import validate
from scripts.resume_cache import reuse_completed,fingerprint,reference_signature

@pytest.fixture
def module_data(tmp_path):
    generate(tmp_path)
    counts=pd.read_csv(tmp_path/'counts.csv',index_col=0)
    np.log2(counts+1).to_csv(tmp_path/'expression.csv')
    result=pd.DataFrame(dict(gene=counts.index,log2FoldChange=np.linspace(-3,3,len(counts)),pvalue=np.linspace(.001,1,len(counts)),padj=np.minimum(1,np.linspace(.001,1,len(counts))*2),stat=np.linspace(-5,5,len(counts))))
    result.to_csv(tmp_path/'differential.csv',index=False)
    return dict(samples=str(tmp_path/'samples.csv'),expression=str(tmp_path/'expression.csv'),expression_unit='log2_normalized',control='Control',treatment='Treatment'),tmp_path

def test_independent_inputs_do_not_require_counts_or_fastq(module_data):
    cfg,p=module_data
    valid,m=validate(dict(cfg,analysis_stage='expression'))
    assert valid['mode']=='module' and valid['comparisons']==[] and len(m)==6
    with pytest.raises(ValueError,match='已转换'):validate(dict(cfg,analysis_stage='expression',expression_unit='tpm'))
    valid,m=validate(dict(analysis_stage='enrichment',differential_input=str(p/'differential.csv'),control='C',treatment='T'))
    assert valid['run_gsea'] is True and len(m)==0
    json.dumps(valid)
    r=pd.read_csv(p/'differential.csv').drop(columns='stat');r.to_csv(p/'differential.csv',index=False)
    valid,_=validate(dict(analysis_stage='enrichment',differential_input=str(p/'differential.csv'),control='C',treatment='T'))
    assert valid['run_gsea'] is False

@pytest.fixture
def cache(tmp_path):
    source=tmp_path/'source';source.mkdir();(source/'index').mkdir()
    ref=tmp_path/'transcripts.fa';ref.write_text('>tx1\nAACCGG\n')
    raw=tmp_path/'read.fastq.gz';raw.write_bytes(b'synthetic raw reads')
    cfg=dict(mode='fastq',transcripts=str(ref),preprocess='fastp',quality=20,length=30,library='A',threads=2,min_count=10)
    (source/'config.json').write_text(json.dumps(cfg));(source/'status.json').write_text(json.dumps(dict(state='failed')))
    (source/'index/versionInfo.json').write_text('{}');(source/'index-complete.json').write_text(json.dumps(reference_signature(cfg)))
    samples=pd.DataFrame(dict(sample=['S1'],group=['C'],fastq1=[str(raw)]))
    for directory in ['qc/S1','quant/S1']:
        path=source/directory;path.mkdir(parents=True);(path/'complete.json').write_text(json.dumps([fingerprint(raw)]))
    (source/'qc/S1/R1.fastq.gz').write_bytes(b'clean reads');(source/'qc/S1/fastp.json').write_text('{}')
    (source/'quant/S1/quant.sf').write_text('Name\tLength\tEffectiveLength\tTPM\tNumReads\ntx1\t100\t80\t1\t2\n')
    (source/'quant/S1/aux_info').mkdir();(source/'quant/S1/aux_info/meta_info.json').write_text('{}')
    return source,cfg,samples

def test_changed_statistics_reuses_only_completed_upstream_and_preserves_source(cache,tmp_path):
    source,cfg,m=cache;target=tmp_path/'new';target.mkdir()
    report=reuse_completed(source,target,dict(cfg,min_count=20),m)
    assert report['reused']==dict(qc=['S1'],index=True,quant=['S1'])
    (target/'qc/S1/R1.fastq.gz').write_bytes(b'new run modification')
    assert (source/'qc/S1/R1.fastq.gz').read_bytes()==b'clean reads'
    assert not (target/'differential.csv').exists()

@pytest.mark.parametrize('change,expected',[({'quality':25},dict(qc=[],index=True,quant=[])),({'library':'U'},dict(qc=['S1'],index=True,quant=[])),({'reuse_completed':False},dict(qc=[],index=False,quant=[]))])
def test_changed_processing_invalidates_dependent_cache(cache,tmp_path,change,expected):
    source,cfg,m=cache;target=tmp_path/'new';target.mkdir()
    assert reuse_completed(source,target,dict(cfg,**change),m)['reused']==expected

def test_changed_inputs_and_incomplete_quant_not_reused(cache,tmp_path):
    source,cfg,m=cache;target=tmp_path/'new';target.mkdir()
    Path(cfg['transcripts']).write_text('>tx2\nAATTCCGG\n')
    report=reuse_completed(source,target,cfg,m)
    assert report['reused']==dict(qc=['S1'],index=False,quant=[])
    (source/'quant/S1/complete.json').unlink()
    target=tmp_path/'new2';target.mkdir();assert reuse_completed(source,target,cfg,m)['reused']['quant']==[]

def test_saved_settings_and_handoff_are_explicit(module_data,tmp_path,monkeypatch):
    import app as station
    cfg,p=module_data;job='c'*32;run=tmp_path/job;run.mkdir()
    cfg.update(mode='counts',counts=str(p/'counts.csv'),comparisons=[dict(name='comparison',control='Control',treatment='Treatment')])
    (run/'config.json').write_text(json.dumps(cfg));(run/'status.json').write_text(json.dumps(dict(state='failed')))
    pd.read_csv(p/'samples.csv').to_csv(run/'samples.csv',index=False)
    pd.read_csv(p/'counts.csv').to_csv(run/'counts.csv',index=False)
    monkeypatch.setattr(station,'RUNS',tmp_path);client=station.app.test_client()
    loaded=client.get('/api/jobs/'+job+'/settings');assert loaded.status_code==200 and loaded.json['resume_from']==job
    response=client.post('/api/jobs/'+job+'/module-inputs',json=dict(analysis_stage='differential'))
    assert response.status_code==200 and response.json['mode']=='counts'
    assert response.json['counts']==str(run/'counts.csv') and response.json['upstream_from']==job
    assert client.post('/api/jobs/'+job+'/module-inputs',json=dict(analysis_stage='network')).status_code==400

def test_enrichment_can_compute_missing_bh_from_complete_pvalues(module_data):
    cfg,p=module_data
    result=pd.read_csv(p/'differential.csv').drop(columns=['padj','stat']);result.to_csv(p/'differential.csv',index=False)
    valid,m=validate(dict(analysis_stage='enrichment',differential_input=str(p/'differential.csv'),control='C',treatment='T'))
    assert valid['run_gsea'] is False and len(m)==0

def test_handoff_preserves_selected_deg_profile(module_data,tmp_path,monkeypatch):
    import app as station
    cfg,p=module_data;job='d'*32;run=tmp_path/job;run.mkdir()
    cfg.update(mode='counts',deg_profile='strict',comparisons=[dict(name='comparison',control='Control',treatment='Treatment')])
    (run/'config.json').write_text(json.dumps(cfg));(run/'status.json').write_text(json.dumps(dict(state='completed')))
    pd.read_csv(p/'samples.csv').to_csv(run/'samples.csv',index=False)
    q=run/'contrasts/comparison';q.mkdir(parents=True)
    pd.read_csv(p/'differential.csv').to_csv(q/'differential.csv',index=False)
    (q/'active-profile.json').write_text(json.dumps(dict(profile='nominal')))
    monkeypatch.setattr(station,'RUNS',tmp_path)
    data=station.app.test_client().post('/api/jobs/'+job+'/module-inputs',json=dict(analysis_stage='enrichment',contrast='comparison')).json
    assert data['deg_profile']=='nominal' and data['differential_input']==str(q/'differential.csv')
    profile=q/'profiles/nominal';profile.mkdir(parents=True);(q/'differential.csv').rename(profile/'differential.csv')
    data=station.app.test_client().post('/api/jobs/'+job+'/module-inputs',json=dict(analysis_stage='enrichment',contrast='comparison')).json
    assert data['differential_input']==str(profile/'differential.csv')
