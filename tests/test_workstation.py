import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import app as station
from scripts.fixture import generate

@pytest.fixture
def inputs(tmp_path):
 generate(tmp_path)
 return dict(mode='counts',samples=str(tmp_path/'samples.csv'),counts=str(tmp_path/'counts.csv'),control='Control',treatment='Treatment',permutations=100),tmp_path

def test_raw_counts_rejects_tpm(inputs):
 cfg,p=inputs;a=pd.read_csv(p/'counts.csv',index_col=0).astype(float);a.iloc[0,0]=.5;a.to_csv(p/'counts.csv')
 with pytest.raises(ValueError,match='整数'):station.validate(cfg)

def test_sample_names_prevent_shell_paths(inputs):
 cfg,p=inputs;m=pd.read_csv(p/'samples.csv');m.loc[0,'sample']='../../bad';m.to_csv(p/'samples.csv',index=False)
 with pytest.raises(ValueError,match='样本名'):station.validate(cfg)

def test_confounded_batch(inputs):
 cfg,p=inputs;m=pd.read_csv(p/'samples.csv');m['batch']=m.group;m.to_csv(p/'samples.csv',index=False)
 with pytest.raises(ValueError,match='混杂'):station.validate(cfg)

def test_partial_replicates(inputs):
 cfg,p=inputs;m=pd.read_csv(p/'samples.csv').iloc[:4];m.to_csv(p/'samples.csv',index=False)
 with pytest.raises(ValueError,match='重复'):station.validate(cfg)

def test_valid_and_batch(inputs):
 cfg,p=inputs;m=pd.read_csv(p/'samples.csv');m['batch']=['A','B','A','A','B','A'];m.to_csv(p/'samples.csv',index=False)
 valid,m=station.validate(cfg);assert valid['threads']==4 and len(m)==6

def test_gmt_validation(inputs):
 cfg,p=inputs;(p/'bad.gmt').write_text('path\tdesc\tg1\npath\tdesc\tg2\n');cfg['gmt']=str(p/'bad.gmt')
 with pytest.raises(ValueError,match='不能重复'):station.validate(cfg)

def test_api_origin_and_errors(inputs):
 c=station.app.test_client();assert c.get('/').status_code==200
 assert c.post('/api/demo',json={},headers={'Origin':'https://evil.test'}).status_code==403
 assert c.get('/api/jobs/not-a-job').status_code==404
 assert c.post('/api/jobs',json={'samples':'/nonexistent'}).status_code==400
 assert c.get('/',headers={'Host':'evil.test'}).status_code==403

def test_plot_exports_and_parameter_validation(tmp_path):
 from plots import render
 m=pd.DataFrame({'sample':['A','B','C','D'],'group':['C','C','T','T']});m.to_csv(tmp_path/'samples.csv',index=False)
 genes=['g1','g2','g3'];pd.DataFrame({'gene':genes,'symbol':genes,'padj':[.001,.01,.8],'log2FoldChange':[3,-2,.1]}).to_csv(tmp_path/'differential.csv',index=False)
 pd.DataFrame(np.arange(12).reshape(3,4),index=genes,columns=m['sample']).to_csv(tmp_path/'vst.csv')
 render(tmp_path,'pca',{})
 pca=pd.read_csv(tmp_path/'pca.plot-data.csv');assert len(pca)==4 and np.isfinite(pca[['PC1','PC2']].to_numpy()).all()
 render(tmp_path,'volcano',{'title':'Edited title','up':'#ee5555'})
 for fmt in ['svg','png','pdf','plot-data.csv']:assert (tmp_path/('volcano.'+fmt)).stat().st_size>50
 assert 'Edited title' in (tmp_path/'volcano.svg').read_text()
 with pytest.raises(ValueError):render(tmp_path,'volcano',{'dpi':9999})
