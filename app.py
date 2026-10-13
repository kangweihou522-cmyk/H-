"""Local mouse bulk RNA-seq workstation. Bind to loopback only."""
from pathlib import Path
import csv, json, os, re, shutil, subprocess, threading, uuid, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from flask import Flask, request, jsonify, send_from_directory, abort
ROOT=Path(__file__).resolve().parent
RUNS=Path(os.environ.get('RNASEQ_RUNS',ROOT/'runs')).resolve(); RUNS.mkdir(parents=True,exist_ok=True)
app=Flask(__name__,static_folder='static'); app.config['MAX_CONTENT_LENGTH']=2*1024*1024
render_lock=threading.Lock()
pool=ThreadPoolExecutor(max_workers=1); processes={}; cancelled=set(); lock=threading.Lock()

def save(path,obj):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2));tmp.replace(path)
def runpath(job):
 if not re.fullmatch('[a-f0-9]{32}',job):abort(404)
 p=RUNS/job
 if not p.is_dir():abort(404)
 return p
def status(p,**kw):
 with lock:
  f=p/'status.json';v=json.loads(f.read_text()) if f.exists() else {};v.update(kw);save(f,v)
def local_file(value):
 if not isinstance(value,str) or not value.strip():raise ValueError('文件路径不能为空')
 p=Path(value).expanduser().resolve()
 if not p.is_file():raise ValueError(f'文件不存在：{p}')
 return p

def validate(cfg):
 if not isinstance(cfg,dict):raise ValueError('请求必须为 JSON 对象')
 mode=cfg.get('mode','counts')
 if mode not in ('counts','fastq','salmon'):raise ValueError('未知输入类型')
 m=pd.read_csv(local_file(cfg.get('samples')),dtype=str).fillna('')
 if not {'sample','group'}.issubset(m.columns):raise ValueError('样本表需要 sample、group 列')
 if len(m)<4 or m['sample'].duplicated().any() or not all(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]{0,79}',s) for s in m['sample']):raise ValueError('至少4个样本，样本名须唯一且只含字母、数字、下划线、点或短横线')
 control=cfg.get('control'); treatment=cfg.get('treatment')
 if not control or not treatment or control==treatment:raise ValueError('请选择不同的对照组和处理组')
 if set(m.group)!={control,treatment} or min(m.groupby('group').size())<2:raise ValueError('本次对比仅允许两个分组，每组至少2个生物学重复')
 if 'batch' in m and m.batch.nunique()>1:
  if (m.batch=='').any():raise ValueError('batch 不可为空')
  design=pd.get_dummies(m[['batch','group']],drop_first=True,dtype=float);design.insert(0,'intercept',1.)
  if np.linalg.matrix_rank(design.to_numpy())<design.shape[1] or len(m)<=design.shape[1]:raise ValueError('批次与分组混杂或没有残差自由度')
 cfg={**cfg,'mode':mode,'threads':int(cfg.get('threads',4)),'alpha':float(cfg.get('alpha',.05)),'lfc':float(cfg.get('lfc',1)),'permutations':int(cfg.get('permutations',1000)),'gmt_id':cfg.get('gmt_id','ENSEMBL')}
 if not 1<=cfg['threads']<=8 or not 0<cfg['alpha']<1 or not 0<=cfg['lfc']<=10 or not 100<=cfg['permutations']<=10000:raise ValueError('参数超出允许范围')
 if mode=='counts':
  counts=pd.read_csv(local_file(cfg.get('counts')),index_col=0)
  if counts.index.has_duplicates or counts.index.isna().any() or not set(m['sample']).issubset(counts.columns):raise ValueError('计数矩阵基因须唯一，列须包含所有样本')
  a=counts[m['sample']].to_numpy(dtype=float)
  if not np.isfinite(a).all() or (a<0).any() or not np.equal(a,np.floor(a)).all():raise ValueError('DESeq2 需要非负整数原始计数；不能使用 TPM/FPKM')
  if (a.sum(axis=0)==0).any():raise ValueError('样本总计数不能为0')
 else:
  mapping=pd.read_csv(local_file(cfg.get('tx2gene')),sep='\t',dtype=str)
  if mapping.shape[1]!=2 or mapping.isna().any().any() or mapping.iloc[:,0].str.replace(r'\..*$','',regex=True).duplicated().any():raise ValueError('tx2gene 需要带表头的两列，转录本（去版本号后）须唯一')
  if mode=='fastq':
   local_file(cfg.get('transcripts'))
   if cfg.get('library','A') not in ['A','IU','ISF','ISR','U','SF','SR']:raise ValueError('不支持的链特异性类型')
   if 'fastq1' not in m:raise ValueError('FASTQ 模式需要 fastq1 列，双端添加 fastq2 列')
   for col in ('fastq1','fastq2'):
    if col in m:
     if (m[col]=='').any():raise ValueError('同一任务不能混合单端和双端')
     m[col]=[str(local_file(f)) for f in m[col]]
  else:
   if 'quant' not in m:raise ValueError('Salmon 导入需要 quant 列（quant.sf 路径）')
   for f in m['quant']:local_file(f)
 if cfg.get('annotation'):
  ann=pd.read_csv(local_file(cfg['annotation']),sep='\t',dtype=str).fillna('')
  if not {'gene','symbol','entrez','go','ontology'}.issubset(ann.columns):raise ValueError('注释表需要 gene symbol entrez go ontology 五列')
  if not set(ann.ontology).issubset({'BP','CC','MF',''}):raise ValueError('GO ontology 只能为 BP/CC/MF')
 if cfg.get('gmt'):
  gmt=local_file(cfg['gmt'])
  if cfg.get('gmt_id','ENSEMBL') not in ('ENSEMBL','ENTREZID','SYMBOL'):raise ValueError('GMT ID 类型无效')
  terms=set()
  for line in gmt.read_text().splitlines():
   parts=line.split('\t')
   if len(parts)<3 or parts[0] in terms:raise ValueError('GMT 每行至少3列，通路名不能重复')
   terms.add(parts[0])
  if not terms:raise ValueError('GMT 文件为空')
 for key in ('counts','transcripts','tx2gene','gmt','annotation'):
  if cfg.get(key):cfg[key]=str(local_file(cfg[key]))
 return cfg,m

def command(p,args,mounts=()):
 if p.name in cancelled:raise RuntimeError('任务已取消')
 exe=args[0]
 if shutil.which(exe):cmd=args
 else:
  cmd=['docker','run','--rm','--name','rnaseq-'+p.name,'--cpus',str(json.loads((p/'config.json').read_text())['threads']),'--memory','24g','--user',f'{os.getuid()}:{os.getgid()}','-e','HOME=/tmp','-v',f'{p}:{p}','-v',f'{ROOT}/scripts:{ROOT}/scripts:ro']
  for directory in sorted(set(str(Path(f).resolve().parent) for f in mounts)):
   if directory!=str(p):cmd+=['-v',f'{directory}:{directory}:ro']
  cmd+=['mouse-rnaseq:1']+args
 with (p/'execution.jsonl').open('a') as f:f.write(json.dumps({'time':time.time(),'command':cmd})+'\n')
 with (p/'log.txt').open('a') as log:
  proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);processes[p.name]=proc
  save(p/'runtime.json',{'pid':proc.pid,'argv':cmd})
  try:
   rc=proc.wait()
   if rc:raise RuntimeError(f'{exe} 退出码 {rc}；请查看日志')
  finally:
   processes.pop(p.name,None);(p/'runtime.json').unlink(missing_ok=True)
 if p.name in cancelled:raise RuntimeError('任务已取消')

def execute(p):
 try:
  cfg=json.loads((p/'config.json').read_text());m=pd.read_csv(p/'samples.csv',dtype=str)
  if cfg['mode']=='fastq':
   status(p,state='running',stage='Salmon 建索引')
   if not (p/'index'/'versionInfo.json').exists():command(p,['salmon','index','-t',cfg['transcripts'],'-i',str(p/'index'),'-p',str(cfg['threads'])], [cfg['transcripts']])
   for _,s in m.iterrows():
    if (p/'quant'/s['sample']/'quant.sf').exists() and (p/'quant'/s['sample']/'aux_info/meta_info.json').exists():continue
    status(p,stage='Salmon 定量 '+s['sample'])
    args=['salmon','quant','-i',str(p/'index'),'-l',cfg.get('library','A'),'-p',str(cfg['threads']),'--validateMappings','--seqBias','--gcBias','-o',str(p/'quant'/s['sample'])]
    files=[s.fastq1];args+=['-1',s.fastq1,'-2',s.fastq2] if 'fastq2' in m else ['-r',s.fastq1]
    if 'fastq2' in m:files.append(s.fastq2)
    command(p,args,files)
  if cfg['mode'] in ('fastq','salmon'):
   qc=[]
   for sample in m['sample']:
    q=p/'quant'/sample/'quant.sf';table=pd.read_csv(q,sep='\t')
    item={'sample':sample,'estimated_reads':float(table.NumReads.sum()),'transcripts':len(table)}
    info=p/'quant'/sample/'aux_info/meta_info.json'
    if info.exists():
     meta=json.loads(info.read_text());item.update({k:meta.get(k) for k in ('num_processed','num_mapped','percent_mapped','salmon_version')})
    qc.append(item)
   pd.DataFrame(qc).to_csv(p/'salmon-qc.csv',index=False)
  status(p,state='running',stage='DESeq2 / GO / GSEA')
  command(p,['Rscript',str(ROOT/'scripts/analyze.R'),str(p)])
  status(p,stage='生成图表')
  from plots import render
  for kind in ('pca','volcano','heatmap','go','gsea')+ (('pathway','pathway_gsea') if cfg.get('gmt') else ()):
   with render_lock:render(p,kind,{})
  status(p,state='completed',stage='完成（GO 未运行：缺少注释）' if (p/'annotation-warning.txt').exists() else '完成',finished=time.time(),warnings=['GO 未运行：请安装 org.Mm.eg.db 或提供注释 TSV'] if (p/'annotation-warning.txt').exists() else [])
 except Exception as e:status(p,state='cancelled' if p.name in cancelled else 'failed',stage=str(e),finished=time.time())

@app.before_request
def protect():
 # Loopback application: prevent cross-origin websites from launching local commands.
 if request.host.split(':')[0] not in ('127.0.0.1','localhost','[::1]'):abort(403)
 if request.method!='GET' and request.headers.get('Origin') and request.headers['Origin']!=request.host_url.rstrip('/'):abort(403)
 if request.method!='GET' and not request.is_json:abort(415)
@app.errorhandler(ValueError)
def invalid(e):return jsonify(error=str(e)),400
@app.errorhandler(FileNotFoundError)
def missing(e):return jsonify(error=str(e)),400
@app.get('/')
def index():return send_from_directory(ROOT/'static','index.html')
@app.get('/api/health')
def health():
 image=subprocess.run(['docker','image','inspect','mouse-rnaseq:1'],capture_output=True).returncode==0 if shutil.which('docker') else False
 return jsonify(python=True,salmon=bool(shutil.which('salmon')),Rscript=bool(shutil.which('Rscript')),docker_image=image,run_dir=str(RUNS))
@app.get('/api/jobs')
def jobs():return jsonify([{'id':p.name,**json.loads((p/'status.json').read_text())} for p in sorted(RUNS.iterdir(),key=lambda x:x.stat().st_mtime,reverse=True) if (p/'status.json').is_file()])
@app.post('/api/jobs')
def create():
 cfg,m=validate(request.get_json());p=RUNS/uuid.uuid4().hex;p.mkdir()
 m.to_csv(p/'samples.csv',index=False);save(p/'config.json',cfg)
 if cfg.get('annotation'):shutil.copyfile(cfg['annotation'],p/'annotation.tsv')
 if cfg.get('gmt'):shutil.copyfile(cfg['gmt'],p/'pathways.gmt')
 if cfg['mode']=='counts':shutil.copyfile(cfg['counts'],p/'counts.csv')
 else:
  shutil.copyfile(cfg['tx2gene'],p/'tx2gene.tsv')
  if cfg['mode']=='salmon':
   for _,s in m.iterrows():
    q=p/'quant'/s['sample'];q.mkdir(parents=True);shutil.copyfile(s['quant'],q/'quant.sf')
    original=Path(s['quant']).parent
    for name in ('aux_info/meta_info.json','lib_format_counts.json','cmd_info.json'):
     if (original/name).is_file():
      (q/name).parent.mkdir(exist_ok=True);shutil.copyfile(original/name,q/name)
 status(p,state='queued',stage='等待',name=cfg.get('name','小鼠 RNA-seq'),created=time.time())
 pool.submit(execute,p);return jsonify(id=p.name),202
@app.get('/api/jobs/<job>')
def detail(job):
 p=runpath(job);return jsonify(**json.loads((p/'status.json').read_text()),log=(p/'log.txt').read_text(errors='replace')[-24000:] if (p/'log.txt').exists() else '',files=[f.name for f in p.iterdir() if f.is_file()])
@app.post('/api/jobs/<job>/retry')
def retry(job):
 p=runpath(job)
 if json.loads((p/'status.json').read_text())['state'] not in ('failed','cancelled'):raise ValueError('仅能重试失败或取消的任务')
 cancelled.discard(job);status(p,state='queued',stage='重试，复用已完成定量')
 pool.submit(execute,p);return jsonify(id=job),202
@app.post('/api/jobs/<job>/cancel')
def cancel(job):
 p=runpath(job)
 if json.loads((p/'status.json').read_text())['state'] in ('completed','failed','cancelled'):raise ValueError('任务已结束')
 cancelled.add(job)
 proc=processes.get(job)
 if proc:
  import signal
  os.killpg(proc.pid,signal.SIGTERM)
  subprocess.run(['docker','stop','rnaseq-'+job],capture_output=True) if shutil.which('docker') else None
 return jsonify(ok=True)
@app.get('/api/jobs/<job>/files/<name>')
def download(job,name):
 p=runpath(job)
 if not re.fullmatch(r'[\w.-]+',name) or name.endswith('.tmp'):abort(404)
 return send_from_directory(p,name,as_attachment=request.args.get('download')=='1')
@app.post('/api/jobs/<job>/plot')
def plot(job):
 p=runpath(job)
 if json.loads((p/'status.json').read_text())['state']!='completed':raise ValueError('请等待分析完成')
 from plots import render
 cfg=request.get_json();kind=cfg.pop('kind','volcano')
 with render_lock:render(p,kind,cfg)
 return jsonify(file=kind+'.svg')
@app.post('/api/demo')
def demo():
 from scripts.fixture import generate
 p=RUNS/'demo_inputs';p.mkdir(exist_ok=True);generate(p)
 return jsonify(samples=str(p/'samples.csv'),counts=str(p/'counts.csv'),control='Control',treatment='Treatment',mode='counts')
def recover_interrupted():
 import signal
 for p in RUNS.iterdir():
  if not (p/'status.json').is_file():continue
  if json.loads((p/'status.json').read_text()).get('state') not in ('queued','running'):continue
  rt=p/'runtime.json'
  if rt.exists():
   record=json.loads(rt.read_text());procfile=Path('/proc')/str(record['pid'])/'cmdline'
   if procfile.exists() and procfile.read_bytes().decode().rstrip('\0').split('\0')==record['argv']:
    try:os.killpg(record['pid'],signal.SIGTERM)
    except ProcessLookupError:pass
   if record['argv'][0]=='docker':subprocess.run(['docker','stop','rnaseq-'+p.name],capture_output=True)
   rt.unlink(missing_ok=True)
  status(p,state='failed',stage='服务重启导致任务中断；可重试')
if __name__=='__main__':
 recover_interrupted()
 app.run(host='127.0.0.1',port=int(os.environ.get('PORT',8765)),threaded=True)
