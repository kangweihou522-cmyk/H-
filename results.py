"""Typed result exploration and a portable offline report, without serving source files."""
from pathlib import Path
import base64, html, json, re, zipfile
import numpy as np
import pandas as pd
from scipy import stats

TABLES = {'fastp-qc','differential','go','gsea','pathway','pathway_gsea','pathway-activity-differential','tf-activity-differential','tf-activity-details','sample-qc','salmon-qc','dispersion','wgcna-modules','wgcna-module-trait','wgcna-soft-threshold','pathway-trait'}
CONTRAST_TABLES={'differential','go','gsea','pathway','pathway_gsea','pathway-activity-differential','tf-activity-differential'}

def contrast_dir(p, name):
    if not name: return p
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', name): raise ValueError('比较名无效')
    q=p/'contrasts'/name
    if not q.is_dir(): raise ValueError('比较不存在')
    return q

def read_json(p, default=None):
    return json.loads(p.read_text()) if p.is_file() else default

def records(df):
    return json.loads(df.to_json(orient='records',double_precision=15))

def manifest(p):
    from plots import chart_catalog
    cfg=read_json(p/'config.json',{})
    summary=read_json(p/'analysis-summary.json',{})
    contrasts=summary.get('contrasts',{})
    if not contrasts:
        contrasts={'legacy':dict(name='legacy',control=cfg.get('control'),treatment=cfg.get('treatment'))}
    charts=chart_catalog(p)
    for item in charts:
        q=contrast_dir(p,item.get('contrast'))
        item['settings']=read_json(q/(item['kind']+'.settings.json'),{})
        item['file']=(str(q.relative_to(p))+'/' if q!=p else '')+item['kind']+'.svg'
    return dict(config=cfg,summary=summary,contrasts=contrasts,charts=charts,mining=read_json(p/'mining-status.json',{}),samples=records(pd.read_csv(p/'samples.csv')),report=(p/'report.html').exists(),multiqc=(p/'multiqc/multiqc_report.html').exists())

def table_page(p,name,args):
    if name not in TABLES: raise ValueError('未知结果表')
    q=contrast_dir(p,args.get('contrast')) if name in CONTRAST_TABLES else p
    f=q/(name+'.csv')
    if not f.exists(): return dict(columns=[],rows=[],total=0,available=False)
    frame=pd.read_csv(f)
    if name in ('go','pathway') and 'direction' in frame:
        frame=frame[frame.direction==args.get('direction','all')]
    query=args.get('query','')[:200]
    if query:
        mask=np.zeros(len(frame),dtype=bool)
        for c in frame.select_dtypes(include=['str','object']).columns:
            mask |= frame[c].astype(str).str.contains(query,case=False,regex=False).to_numpy()
        frame=frame[mask]
    significant=args.get('significant')=='1'
    if significant:
        if 'significant' in frame: frame=frame[frame.significant==True]
        elif 'padj' in frame: frame=frame[frame.padj<=float(args.get('alpha',.05))]
        elif 'adj.P.Val' in frame: frame=frame[frame['adj.P.Val']<=float(args.get('alpha',.05))]
    column=args.get('sort','padj' if 'padj' in frame else frame.columns[0] if len(frame.columns) else '')
    if column in frame: frame=frame.sort_values(column,ascending=args.get('order','asc')=='asc',na_position='last')
    offset=max(0,int(args.get('offset',0)));limit=min(500,max(1,int(args.get('limit',25))))
    return dict(columns=list(frame.columns),rows=records(frame.iloc[offset:offset+limit]),total=len(frame),available=True)

def gene_data(p,query,contrast=None):
    if not query or len(query)>160: raise ValueError('请输入基因 ID 或 Symbol')
    q=contrast_dir(p,contrast);r=pd.read_csv(q/'differential.csv')
    mask=r.gene.astype(str).str.casefold()==query.casefold()
    if 'symbol' in r: mask|=r.symbol.fillna('').astype(str).str.casefold()==query.casefold()
    annotation=p/'gene-annotation.csv'
    if annotation.exists():
        ann=pd.read_csv(annotation,dtype={'gene':str,'id':str,'entrez':str,'symbol':str})
        matching=np.zeros(len(ann),dtype=bool)
        for key in ('id','entrez','symbol'):
            if key in ann:matching|=ann[key].fillna('').astype(str).str.casefold().eq(query.casefold()).to_numpy()
        mask|=r.gene.isin(ann.loc[matching,'gene'])
    candidates=r[mask]
    if candidates.empty: raise ValueError('本次保留的表达基因中未找到该 ID/Symbol')
    gene=str(candidates.iloc[0].gene)
    m=pd.read_csv(p/'samples.csv');norm=pd.read_csv(p/'normalized.csv',index_col=0);vst=pd.read_csv(p/'vst.csv',index_col=0)
    values=m.copy();values['normalized']=norm.loc[gene,m['sample']].to_numpy();values['vst']=vst.loc[gene,m['sample']].to_numpy()
    tpm=p/'tpm.csv'
    if tpm.exists():
        t=pd.read_csv(tpm,index_col=0)
        if gene in t.index: values['TPM']=t.loc[gene,m['sample']].to_numpy()
    associations=[]
    # Exploratory gene/phenotype association, not a substitute for the covariate-adjusted DE model.
    for col in m:
        if col in ('sample','group','fastq1','fastq2','quant'): continue
        numeric=pd.to_numeric(m[col],errors='coerce')
        if numeric.notna().all() and numeric.nunique()>1:
            estimate,pvalue=stats.spearmanr(numeric,values['vst'])
            associations.append(dict(trait=col,spearman=float(estimate) if np.isfinite(estimate) else None,pvalue=float(pvalue) if np.isfinite(pvalue) else None,n=len(m)))
    return dict(gene=gene,candidates=records(candidates),samples=records(values),associations=associations,association_note='探索性 Spearman 关联：未调整协变量，未作全基因多重检验，不代表因果或验证结果。')

def curve_data(p,database,term,contrast=None):
    if database not in ('go','pathway'): raise ValueError('未知基因集数据库')
    q=contrast_dir(p,contrast);name='gsea' if database=='go' else 'pathway_gsea'
    frame=pd.read_csv(q/(name+'.csv'))
    if not term and len(frame.dropna(subset=['padj','NES'])): term=str(frame.dropna(subset=['padj','NES']).sort_values('padj').iloc[0].term)
    result=frame[frame.term==term]
    if result.empty: raise ValueError('该基因集没有可用的 GSEA 结果')
    sets=read_json(p/'gene-sets.json',{})
    members=sets.get(database,{}).get(term,[])
    if isinstance(members,str): members=[members]
    rank=pd.read_csv(q/(name+'-ranks.csv'));hit=rank.gene.isin(members).to_numpy();scores=rank.stat.to_numpy(float)
    n,m=len(scores),int(hit.sum())
    if not 0<m<n: raise ValueError('基因集大小无法计算富集曲线')
    weight=np.abs(scores);increments=np.where(hit,weight/(weight[hit].sum() or 1),-1/(n-m))
    if weight[hit].sum()==0: increments=np.where(hit,1/m,-1/(n-m))
    running=np.cumsum(increments);leading=set(str(result.iloc[0].get('leadingEdge','')).split(';'))
    return dict(term=term,result=records(result)[0],rank=list(range(1,n+1)),running=running.tolist(),stat=scores.tolist(),gene=rank.gene.tolist(),hits=np.flatnonzero(hit).tolist(),leading=[g for g in rank.gene if g in leading])

def build_report(p):
    from plots import chart_catalog
    cfg=read_json(p/'config.json',{});summary=read_json(p/'analysis-summary.json',{});mining=read_json(p/'mining-status.json',{})
    esc=lambda v:html.escape(str(v),quote=True)
    content=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>RNA-seq 分析报告</title><style>body{font:15px/1.7 system-ui;margin:auto;max-width:1100px;padding:30px;color:#183531}h1,h2{color:#126657}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border-bottom:1px solid #d4e0d9;padding:8px;text-align:left}img{max-width:100%}section{margin:30px 0;padding:20px;border:1px solid #d4e0d9;border-radius:10px}.note{background:#fff6e6;padding:15px}pre{white-space:pre-wrap}</style>']
    content += ['<h1>'+esc(cfg.get('name','小鼠 RNA-seq 分析'))+'</h1>','<p>报告由本地工作站生成；此文件可离线打开。图表编辑会更新报告。统计结果由保存的分析参数决定。</p>']
    if (p/'annotation-warning.txt').exists(): content.append('<p class="note">'+esc((p/'annotation-warning.txt').read_text())+'</p>')
    if cfg.get('deg_profile')=='nominal' or any(c.get('active_profile')=='nominal' for c in summary.get('contrasts',{}).values()): content.append('<p class="note">本次采用未校正 P 值筛选差异基因，没有控制 FDR。</p>')
    sections=[('样本与实验设计',pd.read_csv(p/'samples.csv').drop(columns=['fastq1','fastq2','quant'],errors='ignore'))]
    for file,label in [('fastp-qc.csv','测序质控'),('salmon-qc.csv','Salmon 定量质控'),('sample-qc.csv','表达量与文库诊断')]:
        if (p/file).exists(): sections.append((label,pd.read_csv(p/file)))
    if summary.get('contrasts'): sections.append(('比较组合与实际统计阈值',pd.DataFrame([{k:v for k,v in c.items() if k!='profiles'} for c in summary['contrasts'].values()])))
    if summary.get('contrasts'):
        rows=[dict(comparison=name,**spec) for name,c in summary['contrasts'].items() for spec in c.get('profiles',{}).values()]
        if rows:sections.append(('两种差异筛选方案的数量',pd.DataFrame(rows)))
    for title,frame in sections: content.append('<section><h2>'+esc(title)+'</h2>'+frame.to_html(index=False,escape=True)+'</section>')
    content.append('<section><h2>分析方法与来源</h2><pre>'+esc(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre><p>Salmon selective alignment → tximport → DESeq2；GSEA 使用 fgsea multilevel 和 Wald stat；样本诊断使用 VST；通路活性使用 GSVA ssGSEA 和 limma。相关性/模块关联属于探索性统计，不等同于因果证据。</p></section>')
    for chart in chart_catalog(p):
        q=contrast_dir(p,chart.get('contrast'));file=q/(chart['kind']+'.svg')
        if not file.exists(): continue
        image=base64.b64encode(file.read_bytes()).decode()
        content.append('<section><h2>'+esc(chart['title']+(' · '+chart['contrast'] if chart.get('contrast') else ''))+'</h2><img src="data:image/svg+xml;base64,'+image+'" alt="'+esc(chart['title'])+'"></section>')
    for name,state in mining.items(): content.append('<section><h2>'+esc(name)+'</h2><pre>'+esc(json.dumps(state,ensure_ascii=False,indent=2))+'</pre></section>')
    content.append('<section><h2>分析边界</h2><p>Salmon 转录本定量不产生基因组 BAM，未运行的基因体覆盖、内含子/外显子分布、新转录本、剪接或 SNP 分析不能从此结果推算。需要这些结果时使用 nf-core/rnaseq 等基因组比对流程，并保留独立结果。</p></section></html>')
    (p/'report.html').write_text('\n'.join(content),encoding='utf-8')

def bundle(p):
    build_report(p);target=p/'analysis-results.zip';temporary=p/'analysis-results.tmp'
    with zipfile.ZipFile(temporary,'w',zipfile.ZIP_DEFLATED) as z:
        for file in p.rglob('*'):
            if not file.is_file() or file.name in ('analysis-results.zip','analysis-results.tmp','runtime.json'): continue
            if file.relative_to(p).parts[0] in ('index','quant','qc'): continue
            z.write(file,str(file.relative_to(p)))
    temporary.replace(target);return target

def select_profile(p,contrast,profile):
    """Select precomputed ORA/DE results; preserve both profiles and the fitted model."""
    import shutil
    from plots import chart_catalog, render
    if profile not in ('nominal','strict'): raise ValueError('未知差异筛选方案')
    q=contrast_dir(p,contrast);source=q/'profiles'/profile
    if not source.is_dir(): raise ValueError('旧任务没有双方案结果，请创建新任务')
    spec=dict(profile=profile,column='pvalue' if profile=='nominal' else 'padj',alpha=.05,lfc=0 if profile=='nominal' else 1)
    r=pd.read_csv(source/'differential.csv');spec.update(total=int(r.significant.sum()),up=int((r.direction=='up').sum()),down=int((r.direction=='down').sum()))
    for file in ('differential.csv','deg-all.csv','deg-up.csv','deg-down.csv','go.csv','pathway.csv'):
        shutil.copyfile(source/file,q/file)
    (q/'active-profile.json').write_text(json.dumps(spec,indent=2))
    f=p/'analysis-summary.json';summary=read_json(f,{})
    value=summary['contrasts'][contrast];value.update(active_profile=profile,significance=spec['column'],alpha=.05,lfc=spec['lfc'],up=spec['up'],down=spec['down'])
    f.write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    sf=p/'status.json';status=read_json(sf,{})
    warnings=[w for w in status.get('warnings',[]) if '未校正 P 值' not in w]
    if any(v.get('active_profile')=='nominal' for v in summary['contrasts'].values()):warnings.append('采用未校正 P 值筛选的比较未控制 FDR；属于探索性结果。')
    status['warnings']=warnings;sf.write_text(json.dumps(status,ensure_ascii=False,indent=2))
    if list(summary['contrasts'])[0]==contrast:
        for file in ('differential.csv','deg-all.csv','deg-up.csv','deg-down.csv','go.csv','pathway.csv'):shutil.copyfile(q/file,p/file)
    for entry in chart_catalog(p):
        if entry.get('contrast')==contrast and entry['kind'] in ('volcano','ma','heatmap','go','pathway'):
            options={'p_column':spec['column'],'alpha':.05,'lfc':spec['lfc']} if entry['kind'] in ('volcano','ma','heatmap') else {}
            render(q,entry['kind'],options)
    build_report(p)
    return spec
