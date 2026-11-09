"""Meaningful design, IO and cross-comparison invariants for the upgraded workflow."""
import gzip,io,json
from pathlib import Path
import numpy as np,pandas as pd,pytest
from designs import generate_comparisons,design_matrix
from validation import validate,read_table
from scripts.fixture import generate
from workspace_io import scan_fastq,merge_lanes,save_sample_draft,excel_export

@pytest.fixture
def basic(tmp_path):
 generate(tmp_path)
 return dict(mode='counts',samples=str(tmp_path/'samples.csv'),counts=str(tmp_path/'counts.csv'),control='Control',treatment='Treatment'),tmp_path

def test_pairing_requires_explicit_model_and_paired_mode_builds_it(basic):
 cfg,p=basic;m=pd.read_csv(p/'samples.csv');m['subject']=['M1','M2','M3']*2;m.to_csv(p/'samples.csv',index=False)
 with pytest.raises(ValueError,match='重复对象'):validate(cfg)
 valid,m=validate(dict(cfg,design_type='paired'))
 assert valid['covariates']==[dict(name='subject',type='categorical')]
 matrix,names=design_matrix(valid,m);assert matrix.shape==(6,4) and np.linalg.matrix_rank(matrix)==4
 m.loc[1,'subject']='M1';m.to_csv(p/'samples.csv',index=False)
 with pytest.raises(ValueError,match='同组'):validate(dict(cfg,design_type='paired'))

def test_multigroup_and_interaction_hypotheses_are_explicit(basic):
 cfg,p=basic
 rows=[]
 for g in ['C','T']:
  for sex in ['F','M']:
   rows.extend([dict(sample=f'{g}{sex}{i}',group=g,sex=sex) for i in range(3)])
 m=pd.DataFrame(rows);m.to_csv(p/'samples.csv',index=False);pd.DataFrame(np.random.default_rng(8).integers(10,100,(50,12)),columns=m['sample']).to_csv(p/'counts.csv')
 cmps=generate_comparisons(['C','T'],'control','C',['F','M'],'F');assert len(cmps)==3 and cmps[-1]['kind']=='interaction'
 valid,m=validate(dict(cfg,comparisons=cmps,design_type='interaction',interaction_factor='sex',interaction_reference='F'))
 x,names=design_matrix(valid,m);assert x.shape==(12,4) and any(':' in n for n in names)
 assert len(generate_comparisons(['C','M','T'],'all'))==3
 with pytest.raises(ValueError,match='超过'):generate_comparisons([str(i) for i in range(9)],'all')
 with pytest.raises(ValueError,match='重复'):validate(dict(cfg,control='C',treatment='T',comparisons=[dict(name='a',control='C',treatment='T'),dict(name='b',control='C',treatment='T')]))

def test_excel_and_gzip_lane_drafts_preserve_inputs(tmp_path):
 for lane in ['L001','L002']:
  for end in [1,2]:
   with gzip.open(tmp_path/f'A_S1_{lane}_R{end}_001.fastq.gz','wb') as f:f.write(b'@synthetic\nAA\n+\nII\n')
 scan=scan_fastq(tmp_path);assert len(scan['samples'])==2 and scan['lanes_require_merge']
 for r in scan['samples']:r['group']='C'
 with pytest.raises(ValueError,match='lane'):save_sample_draft(tmp_path,scan['samples'])
 rows=merge_lanes(tmp_path,scan['samples']);assert len(rows)==1 and gzip.open(rows[0]['fastq1']).read().count(b'@synthetic')==2
 saved=save_sample_draft(tmp_path,rows);assert Path(saved).exists()
 table=pd.DataFrame(dict(sample=['001','002'],group=['C','T']));table.to_excel(tmp_path/'samples.xlsx',index=False)
 assert read_table(tmp_path/'samples.xlsx',dtype=str)['sample'].tolist()==['001','002']


def test_export_text_ids_formulas_and_full_precision(tmp_path):
 pd.DataFrame(dict(gene=['001','=SUM(1,2)'],pvalue=[1e-20,.0134567891234])).to_csv(tmp_path/'differential.csv',index=False)
 (tmp_path/'config.json').write_text(json.dumps(dict(mode='counts')))
 output=excel_export(tmp_path)
 from openpyxl import load_workbook
 w=load_workbook(output);s=w['differential'];assert s['A2'].value=='001' and s['A3'].data_type=='s'
 assert s['B2'].value==1e-20 and s['B3'].value==.0134567891234


def test_crosscomparison_profile_and_direction_are_not_mixed(tmp_path):
 from advanced_plots import data
 from plots import get_settings
 pd.DataFrame(dict(sample=['s1','s2','s3','s4'],group=['C','C','T','T'])).to_csv(tmp_path/'samples.csv',index=False)
 for name,genes in [('A',['g1','g2']),('B',['g2','g3'])]:
  q=tmp_path/'contrasts'/name/'profiles/strict';q.mkdir(parents=True)
  r=pd.DataFrame(dict(gene=genes,significant=[True,True],log2FoldChange=[2,-2],padj=[.01,.01]))
  r.to_csv(q/'differential.csv',index=False);r.to_csv(q.parent.parent/'differential.csv',index=False)
 cfg=get_settings(tmp_path,'deg_upset');r,names=data(tmp_path,'deg_upset',cfg)
 assert sorted(r['count'])==[1,1,1] and r['count'].sum()==3
 cfg['direction']='up';r,_=data(tmp_path,'deg_upset',cfg);assert r['count'].sum()==2


def test_composition_recovers_known_rna_mixture_and_rejects_collinear_reference():
 from composition import estimate
 rng=np.random.default_rng(15);reference=pd.DataFrame(rng.uniform(1,100,(100,3)),index=[f'g{i}' for i in range(100)],columns=['A','B','C'])
 reference=reference.div(reference.sum())*1e6;truth=np.array([[.2,.6],[.3,.3],[.5,.1]])
 expression=pd.DataFrame(reference.to_numpy()@truth,index=reference.index,columns=['S1','S2'])
 fractions,diagnostics=estimate(expression,reference)
 assert np.allclose(fractions[['A','B','C']].to_numpy(),truth.T,atol=1e-9) and diagnostics.cosine.min()>.99999
 reference['C']=reference['B']
 with pytest.raises(ValueError,match='区分'):estimate(expression,reference)


def test_resource_snapshots_normalize_real_excel_without_overwriting_source(basic,monkeypatch,tmp_path):
 import app as station
 cfg,p=basic
 annotation=pd.DataFrame(dict(gene=['ENSMUSG00000000001'],symbol=['G'],entrez=['001'],go=['GO:0006955'],ontology=['BP']))
 annotation.to_excel(p/'annotation.xlsx',index=False)
 runroot=p/'runs';runroot.mkdir();monkeypatch.setattr(station,'RUNS',runroot);monkeypatch.setattr(station.pool,'submit',lambda *a,**k:None)
 response=station.app.test_client().post('/api/jobs',json=dict(cfg,annotation=str(p/'annotation.xlsx')))
 assert response.status_code==202,response.json
 snapshot=runroot/response.json['id']/'annotation.tsv'
 assert snapshot.read_bytes()[:2]!=b'PK' and pd.read_csv(snapshot,sep='\t',dtype=str).entrez.iloc[0]=='001'
 assert (p/'annotation.xlsx').read_bytes()[:2]==b'PK'


def test_single_end_draft_drops_unused_r2_column(tmp_path):
 path=save_sample_draft(tmp_path,[dict(sample='S1',group='C',fastq1='/reads/R1.fastq.gz',fastq2='')])
 assert 'fastq2' not in pd.read_csv(path)


def test_empty_deg_heatmap_is_explained_without_annotation_crash(tmp_path):
 from plots import interactive_data
 m=pd.DataFrame(dict(sample=['C1','C2','T1','T2'],group=['C','C','T','T']));m.to_csv(tmp_path/'samples.csv',index=False)
 pd.DataFrame(np.arange(12).reshape(3,4),index=['g1','g2','g3'],columns=m['sample']).to_csv(tmp_path/'vst.csv')
 pd.DataFrame(dict(gene=['g1','g2','g3'],log2FoldChange=[0,.2,-.1],padj=[1,.5,.8],significant=[False]*3)).to_csv(tmp_path/'differential.csv',index=False)
 payload=interactive_data(tmp_path,'heatmap',{})
 assert payload['data'][0]['z']==[] and '没有符合' in payload['layout']['annotations'][0]['text']


def test_custom_pair_ids_keep_strings_and_avoid_independent_association_pvalues(tmp_path):
 from results import gene_data
 from plots import interactive_data
 (tmp_path/'config.json').write_text(json.dumps(dict(subject_column='mouse',mode='counts')))
 m=pd.DataFrame(dict(sample=['C1','C2','T1','T2'],group=['C','C','T','T'],mouse=['001','002','001','002'],score=[1,2,3,4]));m.to_csv(tmp_path/'samples.csv',index=False)
 pd.DataFrame([[1,2,3,4]],index=['g1'],columns=m['sample']).to_csv(tmp_path/'vst.csv')
 pd.DataFrame(dict(gene=['g1'],padj=[.01],baseMean=[10],log2FoldChange=[1])).to_csv(tmp_path/'differential.csv',index=False)
 data=gene_data(tmp_path,'g1');assert data['samples'][0]['mouse']=='001' and all(a['pvalue'] is None for a in data['associations'])
 d=interactive_data(tmp_path,'gene_expression',{},dict(gene='g1',connect_pairs=True,gene_scale='vst'))
 assert len([t for t in d['data'] if t.get('mode')=='lines'])==2


def test_named_plot_recipes_and_project_theme_keep_expression_unchanged(tmp_path,monkeypatch):
 import app as station
 job='e'*32;p=tmp_path/job;p.mkdir()
 m=pd.DataFrame(dict(sample=['C1','C2','T1','T2'],group=['C','C','T','T']));m.to_csv(p/'samples.csv',index=False)
 (p/'config.json').write_text(json.dumps(dict(analysis_stage='expression')));(p/'status.json').write_text(json.dumps(dict(state='completed')))
 pd.DataFrame(np.random.default_rng(7).normal(size=(15,4)),index=[f'g{i}' for i in range(15)],columns=m['sample']).to_csv(p/'vst.csv')
 source=(p/'vst.csv').read_bytes();monkeypatch.setattr(station,'RUNS',tmp_path);client=station.app.test_client()
 r=client.post(f'/api/jobs/{job}/recipes',json=dict(name='paper PCA',kind='pca',settings=dict(font=13)));assert r.status_code==200 and r.json['recipes']['paper PCA']['settings']['font']==13
 r=client.post(f'/api/jobs/{job}/theme',json=dict(font=14,up='#ee6633',apply_all=True));assert r.status_code==200 and not r.json['warnings']
 assert json.loads((p/'pca.settings.json').read_text())['font']==14 and (p/'vst.csv').read_bytes()==source
 assert client.post(f'/api/jobs/{job}/theme',json=dict(font=99)).status_code==400
