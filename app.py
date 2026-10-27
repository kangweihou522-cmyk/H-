"""Local mouse bulk RNA-seq workstation, version 2."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json, os, re, shutil, signal, subprocess, threading, time, uuid, hashlib, sys
import numpy as np
import pandas as pd
from flask import Flask, request, jsonify, send_from_directory, abort, send_file
from validation import validate, local_file, local_path, read_table
ROOT = Path(__file__).resolve().parent
STATION_ROOT = Path(os.environ.get('RNASEQ_HOME', ROOT)).resolve()
RUNS = Path(os.environ.get('RNASEQ_RUNS', ROOT/'runs')).resolve()
RUNS.mkdir(parents=True, exist_ok=True)
app = Flask(__name__, static_folder='static')
app.config['MAX_CONTENT_LENGTH'] = 2*1024*1024
pool = ThreadPoolExecutor(max_workers=1)
processes, cancelled, pending = {}, set(), set()
lock, render_lock = threading.RLock(), threading.Lock()

def save(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    tmp.replace(path)

def runpath(job):
    if not re.fullmatch('[a-f0-9]{32}', job): abort(404)
    p = RUNS/job
    if not p.is_dir(): abort(404)
    return p

def status(p, **kwargs):
    with lock:
        f = p/'status.json'
        value = json.loads(f.read_text()) if f.exists() else {}
        value.update(kwargs)
        save(f, value)

def step(p, name, state, message=''):
    with lock:
        f = p/'steps.json'
        data = json.loads(f.read_text()) if f.exists() else {}
        data[name] = dict(state=state, message=message, updated=time.time())
        save(f, data)
    status(p, stage=message or name)

def contrast_path(p, name=None):
    if not name: return p
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', name): abort(404)
    q = p/'contrasts'/name
    if not q.is_dir(): abort(404)
    return q

def command(p, args, mounts=(), analysis_tool=True):
    if p.name in cancelled: raise RuntimeError('任务已取消')
    cfg = json.loads((p/'config.json').read_text())
    if analysis_tool and cfg.get('executor', 'docker') == 'docker':
        cmd = ['docker', 'run', '--rm', '--name', 'rnaseq-'+p.name,
               '--cpus', str(cfg['threads']), '--memory', str(cfg['memory'])+'g',
               '--user', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp',
               '-e', 'OPENBLAS_NUM_THREADS=1', '-e', 'OMP_NUM_THREADS=1',
               '-e', f'RNASEQ_SCRIPT_DIR={ROOT}/scripts',
               '-v', f'{p}:{p}', '-v', f'{ROOT}/scripts:{ROOT}/scripts:ro']
        for directory in sorted({str(Path(f).resolve().parent) for f in mounts}):
            if directory != str(p): cmd += ['-v', f'{directory}:{directory}:ro']
        cmd += ['mouse-rnaseq:2'] + args
    else:
        cmd = args
    with (p/'execution.jsonl').open('a') as f:
        f.write(json.dumps(dict(time=time.time(), command=cmd))+'\n')
    with (p/'log.txt').open('a') as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        with lock: processes[p.name] = proc
        try:
            save(p/'runtime.json', dict(pid=proc.pid, argv=cmd))
            rc = proc.wait()
            if rc: raise RuntimeError(f'{Path(args[0]).name} 退出码 {rc}；请查看日志')
        finally:
            with lock: processes.pop(p.name, None)
            (p/'runtime.json').unlink(missing_ok=True)
    if p.name in cancelled: raise RuntimeError('任务已取消')

def source_version():
    manifest=ROOT/'.source-manifest.json'
    if manifest.is_file():
        files=json.loads(manifest.read_text())
        return 'source-sha256:'+hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest()
    if not shutil.which('git'): return 'v2.0 release archive'
    result=subprocess.run(['git','-C',str(ROOT),'rev-parse','HEAD'],capture_output=True,text=True)
    return result.stdout.strip() if result.returncode==0 else 'v2.0 release archive'

def fastp_summary(p,m):
    rows=[]
    for sample in m['sample']:
        f=p/'qc'/sample/'fastp.json'
        if not f.exists():continue
        data=json.loads(f.read_text());raw=data['summary']['before_filtering'];clean=data['summary']['after_filtering']
        rows.append(dict(sample=sample,raw_reads=raw['total_reads'],clean_reads=clean['total_reads'],raw_bases=raw['total_bases'],clean_bases=clean['total_bases'],q20_rate=clean['q20_rate'],q30_rate=clean['q30_rate'],gc_content=clean['gc_content'],raw_q20_rate=raw['q20_rate'],raw_q30_rate=raw['q30_rate'],raw_gc_content=raw['gc_content'],retained_pct=100*clean['total_reads']/max(1,raw['total_reads']),adapter_trimmed_reads=data.get('adapter_cutting',{}).get('adapter_trimmed_reads',0)))
    if rows:pd.DataFrame(rows).to_csv(p/'fastp-qc.csv',index=False)

def fingerprint(file):
    path = Path(file); s = path.stat()
    return dict(path=str(path.resolve()), size=s.st_size, mtime_ns=s.st_mtime_ns)

def salmon_quant_args(cfg,index,output,files):
    args=['salmon','quant','-i',str(index),'-l',cfg.get('library','A'),'-p',str(cfg['threads']),
          '--validateMappings','--seqBias','--gcBias','-o',str(output)]
    if len(files)==2: return args+['-1',files[0],'-2',files[1]]
    return args+['-r',files[0],'--fldMean',str(cfg.get('fragment_mean',200)),'--fldSD',str(cfg.get('fragment_sd',80))]

def execute(p):
    try:
        if p.name in cancelled: raise RuntimeError('任务已取消')
        cfg = json.loads((p/'config.json').read_text())
        m = pd.read_csv(p/'samples.csv', dtype=str)
        status(p, state='running')
        step(p, 'inputs', 'completed', '输入检查完成')
        if cfg['mode'] == 'fastq':
            if cfg['preprocess'] == 'fastp':
                step(p, 'qc', 'running', 'fastp 接头/质量过滤')
                for _, s in m.iterrows():
                    q = p/'qc'/s['sample']; q.mkdir(parents=True, exist_ok=True)
                    marker = q/'complete.json'
                    files = [s.fastq1] + ([s.fastq2] if 'fastq2' in m else [])
                    signature = [fingerprint(f) for f in files]
                    if marker.exists() and json.loads(marker.read_text()) == signature: continue
                    args = ['fastp', '-i', s.fastq1, '-o', str(q/'R1.fastq.gz'),
                            '-w', str(min(8,cfg['threads'])), '-q', str(cfg['quality']),
                            '-l', str(cfg['length']), '--json', str(q/'fastp.json'), '--html', str(q/'fastp.html')]
                    if 'fastq2' in m:
                        args += ['-I', s.fastq2, '-O', str(q/'R2.fastq.gz'), '--detect_adapter_for_pe']
                    command(p, args, files)
                    save(marker, signature)
                fastp_summary(p,m)
                step(p, 'qc', 'completed', 'fastp 质量过滤完成')
            else: step(p, 'qc', 'skipped', '输入为已处理 FASTQ，跳过 fastp')
            step(p, 'quant', 'running', 'Salmon 建索引与定量')
            signature = dict(transcripts=fingerprint(cfg['transcripts']),genome=fingerprint(cfg['genome']) if cfg.get('genome') else None)
            marker = p/'index-complete.json'
            if not marker.exists() or json.loads(marker.read_text()) != signature:
                transcript=cfg['transcripts'];extra=[]
                if cfg.get('genome'):
                    command(p,[sys.executable,str(ROOT/'scripts/prepare_decoy.py'),transcript,cfg['genome'],str(p/'reference')],analysis_tool=False)
                    transcript=str(p/'reference/gentrome.fa');extra=['--decoys',str(p/'reference/decoys.txt')]
                command(p, ['salmon','index','-t',transcript,'-i',str(p/'index'),'-p',str(cfg['threads'])]+extra, [transcript])
                save(marker, signature)
            for _, s in m.iterrows():
                q = p/'quant'/s['sample']; done = q/'complete.json'
                files = [s.fastq1] + ([s.fastq2] if 'fastq2' in m else [])
                signature = [fingerprint(f) for f in files]
                if done.exists() and json.loads(done.read_text()) == signature: continue
                status(p, stage='Salmon 定量 '+s['sample'])
                if cfg['preprocess'] == 'fastp':
                    qcp = p/'qc'/s['sample']
                    files = [str(qcp/'R1.fastq.gz')] + ([str(qcp/'R2.fastq.gz')] if 'fastq2' in m else [])
                args = salmon_quant_args(cfg,p/'index',q,files)
                command(p,args,files)
                save(done, signature)
            step(p,'quant','completed','Salmon 定量完成')
        else:
            step(p,'qc','skipped','导入模式：原始 reads 质控未运行')
            step(p,'quant','completed','导入计数/定量结果')
        if cfg['mode'] in ('fastq','salmon'):
            qc=[]
            for sample in m['sample']:
                q=p/'quant'/sample; table=pd.read_csv(q/'quant.sf',sep='\t')
                item=dict(sample=sample,estimated_reads=float(table.NumReads.sum()),transcripts=len(table))
                info=q/'aux_info/meta_info.json'
                if info.exists():
                    meta=json.loads(info.read_text());item.update({k:meta.get(k) for k in ('num_processed','num_mapped','percent_mapped','salmon_version')})
                lib=q/'lib_format_counts.json'
                if lib.exists():
                    info=json.loads(lib.read_text());item['expected_format']=info.get('expected_format');item['compatible_ratio']=info.get('compatible_fragment_ratio')
                qc.append(item)
            pd.DataFrame(qc).to_csv(p/'salmon-qc.csv',index=False)
        step(p,'statistics','running','DESeq2 / fgsea / GSVA / WGCNA')
        command(p,['Rscript',str(ROOT/'scripts/analyze.R'),str(p)])
        from discovery import run_discovery
        run_discovery(p)
        step(p,'statistics','completed','统计与数据挖掘完成')
        step(p,'visuals','running','生成诊断图与分析报告')
        from plots import render, chart_catalog
        plot_errors=[]
        for entry in chart_catalog(p):
            target=contrast_path(p,entry.get('contrast'))
            try:
                with render_lock: render(target,entry['kind'],{})
            except (ValueError,KeyError,TypeError,IndexError) as error:
                plot_errors.append(dict(kind=entry['kind'],contrast=entry.get('contrast'),error=str(error)))
        save(p/'plot-warnings.json',plot_errors)
        multiqc=Path(sys.executable).parent/'multiqc'
        if cfg['mode']=='fastq' and multiqc.is_file():
            command(p,[str(multiqc),str(p/'qc'),str(p/'quant'),'-o',str(p/'multiqc'),'--force','--no-ansi','--quiet'],analysis_tool=False)
        from results import build_report
        build_report(p)
        step(p,'visuals','completed','图表与报告已生成')
        warnings=[]
        if plot_errors: warnings.append(f'{len(plot_errors)} 张图未成功导出；统计表已保留，详情见 plot-warnings.json，可调整后重绘。')
        if (p/'annotation-warning.txt').exists(): warnings.append((p/'annotation-warning.txt').read_text().strip())
        if cfg['significance']=='pvalue': warnings.append('采用未校正 P 值筛选；没有控制差异基因的 FDR。')
        if (p/'salmon-qc.csv').exists():
            qc=pd.read_csv(p/'salmon-qc.csv')
            if 'percent_mapped' in qc and (qc.percent_mapped<60).any(): warnings.append('存在 Salmon 映射率低于60%的样本，请检查参考版本、文库方向和污染。')
        status(p,state='completed',stage='完成（有注意事项）' if warnings else '完成',finished=time.time(),warnings=warnings)
    except Exception as e:
        status(p,state='cancelled' if p.name in cancelled else 'failed',stage=str(e),finished=time.time())
        for key,value in json.loads((p/'steps.json').read_text()).items() if (p/'steps.json').exists() else []:
            if value['state']=='running': step(p,key,'failed',str(e))
    finally:
        with lock: pending.discard(p.name)

@app.before_request
def protect():
    if request.host.split(':')[0] not in ('127.0.0.1','localhost','[::1]'): abort(403)
    if request.method!='GET' and request.headers.get('Origin') and request.headers['Origin']!=request.host_url.rstrip('/'): abort(403)
    if request.method!='GET' and not request.is_json: abort(415)

@app.errorhandler(ValueError)
@app.errorhandler(FileNotFoundError)
def invalid(e): return jsonify(error=str(e)),400

@app.get('/')
def index(): return send_from_directory(ROOT/'static','index.html')

@app.get('/api/manual')
def manual_status():
    manifest=ROOT/'docs/guide-manifest.json'
    metadata={}
    if manifest.is_file():
        try: metadata=json.loads(manifest.read_text(encoding='utf-8'))
        except (OSError,ValueError): pass
    return jsonify(html=(ROOT/'static/guide/index.html').is_file(),
                   word=(ROOT/'docs/USER_GUIDE.docx').is_file(),
                   source=(ROOT/'docs/USER_GUIDE.md').is_file(),metadata=metadata)

@app.get('/api/manual/word')
def manual_word():
    file=ROOT/'docs/USER_GUIDE.docx'
    if not file.is_file(): return jsonify(error='当前版本尚未包含Word使用说明，请更新工作站后重试。'),404
    return send_file(file,as_attachment=True,download_name='小鼠bulk-RNA-seq工作站使用说明.docx',
                     mimetype='application/vnd.openxmlformats-officedocument.wordprocessingml.document')

@app.get('/api/manual/source')
def manual_source():
    file=ROOT/'docs/USER_GUIDE.md'
    if not file.is_file(): return jsonify(error='当前版本尚未包含使用说明原稿，请更新工作站后重试。'),404
    return send_file(file,as_attachment=True,download_name='小鼠bulk-RNA-seq工作站使用说明.md',
                     mimetype='text/markdown; charset=utf-8')

@app.get('/api/health')
def health():
    image=subprocess.run(['docker','image','inspect','mouse-rnaseq:2'],capture_output=True).returncode==0 if shutil.which('docker') else False
    return jsonify(service_id='mouse-bulk-rnaseq',version='2.0',browse_root=str(STATION_ROOT if os.environ.get('RNASEQ_HOME') else ROOT.parent),python=True,docker_image=image,salmon=bool(shutil.which('salmon')),Rscript=bool(shutil.which('Rscript')),run_dir=str(RUNS),program_dir=str(ROOT),station_dir=str(STATION_ROOT),windows_dir=os.environ.get('RNASEQ_WINDOWS_HOME',''),python_environment=sys.prefix,pid=os.getpid())

@app.get('/api/files')
def browse():
    roots=[Path.home(),STATION_ROOT,ROOT.parent,Path('/mnt'),Path('/data')]
    p=local_path(request.args.get('path',str(STATION_ROOT if os.environ.get('RNASEQ_HOME') else ROOT.parent)))
    if not any(p==r or p.is_relative_to(r) for r in roots): raise ValueError('请浏览 home、workspace、mnt 或 data 目录')
    if not p.is_dir(): raise ValueError('目录不存在')
    entries=[]
    for item in sorted(p.iterdir(),key=lambda x:(not x.is_dir(),x.name.lower())):
        if item.name.startswith('.'): continue
        try: entries.append(dict(name=item.name,path=str(item),directory=item.is_dir(),size=item.stat().st_size))
        except OSError: continue
        if len(entries)>=1000: break
    return jsonify(path=str(p),parent=str(p.parent),entries=entries)

@app.post('/api/preflight')
def preflight():
    cfg,m=validate(request.get_json())
    warnings=[]
    if min(m.groupby('group').size())<3: warnings.append('存在每组仅2个生物学重复的情况，推荐≥3。')
    if cfg['significance']=='pvalue': warnings.append('未校正P值不能控制FDR，默认建议使用BH padj。')
    if cfg['mode']=='fastq' and not cfg.get('genome'): warnings.append('未提供匹配版本基因组：使用无 decoy 的转录组索引。正式分析建议 full-decoy；本工作站支持可选基因组 FASTA。')
    if cfg['wgcna'] and len(m)<15: warnings.append('WGCNA 将跳过：独立样本不足15个，推荐≥20。')
    if not cfg.get('annotation'): warnings.append('未指定离线注释；GO/Entrez/Symbol映射将尝试容器中的 org.Mm.eg.db。')
    return jsonify(samples=m.to_dict('records'),groups=m.groupby('group').size().to_dict(),comparisons=cfg['comparisons'],covariates=cfg['covariates'],warnings=warnings,output_free_gb=round(shutil.disk_usage(RUNS).free/2**30,1))

@app.post('/api/metadata')
def metadata():
    m=read_table(local_file(request.get_json().get('samples')),dtype=str).fillna('')
    if not {'sample','group'}.issubset(m): raise ValueError('样本表需要 sample 和 group')
    numeric=[c for c in m if c not in ('sample','fastq1','fastq2','quant') and pd.to_numeric(m[c],errors='coerce').notna().all()]
    return jsonify(samples=m.to_dict('records'),groups=list(m.group.unique()),columns=list(m.columns),numeric=numeric)

@app.get('/api/jobs')
def jobs():
    return jsonify([dict(id=p.name,**json.loads((p/'status.json').read_text())) for p in sorted(RUNS.iterdir(),key=lambda x:x.stat().st_mtime,reverse=True) if (p/'status.json').is_file()])

@app.post('/api/jobs')
def create():
    cfg,m=validate(request.get_json());p=RUNS/uuid.uuid4().hex;p.mkdir()
    m.to_csv(p/'samples.csv',index=False);save(p/'config.json',cfg)
    inputs=[fingerprint(cfg['samples'])]
    for key,name in [('annotation','annotation.tsv'),('gmt','pathways.gmt'),('tf_network','regulatory-network.tsv')]:
        if cfg.get(key): shutil.copyfile(cfg[key],p/name);inputs.append(fingerprint(cfg[key]))
    if cfg['mode']=='counts': read_table(cfg['counts'],index_col=0,converters={0:str})[m['sample']].to_csv(p/'counts.csv');inputs.append(fingerprint(cfg['counts']))
    else:
        read_table(cfg['tx2gene'],dtype=str).to_csv(p/'tx2gene.tsv',sep='\t',index=False)
        inputs.append(fingerprint(cfg['tx2gene']))
        if cfg['mode']=='fastq':
            inputs.append(fingerprint(cfg['transcripts']))
            if cfg.get('genome'): inputs.append(fingerprint(cfg['genome']))
        if cfg['mode']=='salmon':
            for _,s in m.iterrows():
                inputs.append(fingerprint(s['quant']))
                q=p/'quant'/s['sample'];q.mkdir(parents=True)
                shutil.copyfile(s['quant'],q/'quant.sf'); original=Path(s['quant']).parent
                for name in ('aux_info/meta_info.json','lib_format_counts.json','cmd_info.json'):
                    if (original/name).is_file():
                        (q/name).parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(original/name,q/name)
    save(p/'input-provenance.json',dict(files=inputs,source_commit=source_version(),created=time.time()))
    status(p,state='queued',stage='等待',input_mode=cfg['mode'],name=cfg.get('name','小鼠 RNA-seq'),created=time.time())
    with lock: pending.add(p.name)
    pool.submit(execute,p)
    return jsonify(id=p.name),202

@app.get('/api/jobs/<job>')
def detail(job):
    p=runpath(job)
    return jsonify(**json.loads((p/'status.json').read_text()),log=(p/'log.txt').read_text(errors='replace')[-24000:] if (p/'log.txt').exists() else '',steps=json.loads((p/'steps.json').read_text()) if (p/'steps.json').exists() else {},files=[f.name for f in p.iterdir() if f.is_file()])

@app.post('/api/jobs/<job>/retry')
def retry(job):
    p=runpath(job)
    with lock:
        if p.name in pending or json.loads((p/'status.json').read_text())['state'] not in ('failed','cancelled'): raise ValueError('仅能重试失败或取消的任务')
        cfg=json.loads((p/'config.json').read_text())
        if 'comparisons' not in cfg: raise ValueError('旧版本任务请导入原结果或创建新的v2任务')
        cancelled.discard(job);pending.add(job);status(p,state='queued',stage='重试，复用已完成步骤')
    pool.submit(execute,p)
    return jsonify(id=job),202

@app.post('/api/jobs/<job>/cancel')
def cancel(job):
    p=runpath(job)
    if json.loads((p/'status.json').read_text())['state'] in ('completed','failed','cancelled'): raise ValueError('任务已结束')
    cancelled.add(job)
    with lock: proc=processes.get(job)
    if proc:
        try: os.killpg(proc.pid,signal.SIGTERM)
        except ProcessLookupError: pass
        if shutil.which('docker'): subprocess.run(['docker','stop','rnaseq-'+job],capture_output=True)
    status(p,state='cancelled',stage='已取消',finished=time.time())
    return jsonify(ok=True)

@app.get('/api/jobs/<job>/files/<path:name>')
def download(job,name):
    p=runpath(job)
    target=(p/name).resolve()
    if not target.is_relative_to(p) or not target.is_file() or target.name.endswith('.tmp'): abort(404)
    return send_from_directory(p,name,as_attachment=request.args.get('download')=='1')

@app.get('/api/jobs/<job>/manifest')
def manifest(job):
    from results import manifest
    return jsonify(manifest(runpath(job)))

@app.get('/api/jobs/<job>/table/<name>')
def table(job,name):
    from results import table_page
    return jsonify(table_page(runpath(job),name,request.args))

@app.get('/api/jobs/<job>/gene')
def gene(job):
    from results import gene_data
    return jsonify(gene_data(runpath(job),request.args.get('gene',''),request.args.get('contrast')))

@app.get('/api/jobs/<job>/curve')
def curve(job):
    from results import curve_data
    return jsonify(curve_data(runpath(job),request.args.get('database','go'),request.args.get('term',''),request.args.get('contrast')))

@app.get('/api/jobs/<job>/chart/<kind>')
def chart(job,kind):
    from plots import interactive_data
    return jsonify(interactive_data(runpath(job),kind,request.args))

@app.post('/api/jobs/<job>/profile')
def choose_profile(job):
    p=runpath(job)
    if json.loads((p/'status.json').read_text())['state']!='completed':raise ValueError('请等待分析完成')
    cfg=request.get_json()
    from results import select_profile
    with render_lock:result=select_profile(p,cfg.get('contrast'),cfg.get('profile'))
    return jsonify(result)

@app.post('/api/jobs/<job>/plot')
def plot(job):
    p=runpath(job)
    if json.loads((p/'status.json').read_text())['state']!='completed': raise ValueError('请等待分析完成')
    from plots import render
    cfg=dict(request.get_json());kind=cfg.pop('kind','volcano');contrast=cfg.pop('contrast',None)
    with render_lock: render(contrast_path(p,contrast),kind,cfg)
    from results import build_report
    build_report(p)
    return jsonify(file=(f'contrasts/{contrast}/' if contrast else '')+kind+'.svg')

@app.post('/api/jobs/<job>/preview')
def preview(job):
    p=runpath(job)
    if json.loads((p/'status.json').read_text())['state']!='completed': raise ValueError('请等待分析完成')
    from plots import interactive_data
    cfg=dict(request.get_json());kind=cfg.pop('kind','volcano');contrast=cfg.pop('contrast',None)
    return jsonify(interactive_data(p,kind,dict(contrast=contrast),cfg))

@app.get('/api/jobs/<job>/chart-download/<kind>/<fmt>')
def chart_download(job,kind,fmt):
    p=runpath(job)
    if json.loads((p/'status.json').read_text())['state']!='completed': raise ValueError('请等待分析完成')
    from plots import LABELS,render
    if kind not in LABELS or fmt not in ('svg','pdf','png','plot-data.csv','plot-data.json','settings.json'):
        raise ValueError('未知图表或导出格式')
    q=contrast_path(p,request.args.get('contrast'));file=q/(kind+'.'+fmt)
    settings=q/(kind+'.settings.json')
    current=json.loads(settings.read_text()) if settings.exists() else {}
    if not file.exists() or current.get('style_version',0)<3:
        with render_lock: render(q,kind,{})
        from results import build_report
        build_report(p)
    if not file.is_file(): raise ValueError('此图使用另一种绘图数据格式，请选择对应 CSV 或 JSON')
    return send_file(file,as_attachment=True,download_name=file.name)

@app.get('/api/jobs/<job>/bundle')
def bundle(job):
    from results import bundle
    p=runpath(job)
    return send_file(bundle(p),as_attachment=True,download_name='rnaseq-'+job[:8]+'.zip')

@app.post('/api/demo')
def demo():
    from scripts.make_smoke_data import create
    p=RUNS/'demo_inputs_v2'
    if not (p/'request.json').exists(): create(p)
    cfg=json.loads((p/'request.json').read_text());cfg.update(mode='counts',counts=str(p/'counts.csv'),wgcna=True,pathway_activity=True,preprocess='fastp',name='合成演示 · 非真实实验数据')
    return jsonify(cfg)

def recover_interrupted():
    for p in RUNS.iterdir():
        if not (p/'status.json').is_file(): continue
        if json.loads((p/'status.json').read_text()).get('state') not in ('queued','running'): continue
        rt=p/'runtime.json'
        if rt.exists():
            rec=json.loads(rt.read_text());f=Path('/proc')/str(rec['pid'])/'cmdline'
            if f.exists() and f.read_bytes().decode().rstrip('\0').split('\0')==rec['argv']:
                try: os.killpg(rec['pid'],signal.SIGTERM)
                except ProcessLookupError: pass
            if rec['argv'][0]=='docker': subprocess.run(['docker','stop','rnaseq-'+p.name],capture_output=True)
            rt.unlink(missing_ok=True)
        status(p,state='failed',stage='服务重启导致任务中断；可重试')

if __name__=='__main__':
    recover_interrupted()
    app.run(host='127.0.0.1',port=int(os.environ.get('PORT',8765)),threaded=True)
