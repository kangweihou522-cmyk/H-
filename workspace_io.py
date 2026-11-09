"""Local data preparation, resource bookmarks and curated exports; no read uploads."""
from pathlib import Path
import hashlib,io,json,re
import pandas as pd
from validation import local_file,local_path,read_table


def scan_fastq(directory):
    root=local_path(str(directory))
    if not root.is_dir():raise ValueError('FASTQ目录不存在')
    rows={};unmatched=[]
    for f in sorted(root.iterdir()):
        if not f.is_file() or not re.search(r'\.(?:fastq|fq)(?:\.gz)?$',f.name,re.I):continue
        name=re.sub(r'\.(?:fastq|fq)(?:\.gz)?$','',f.name,flags=re.I)
        match=re.fullmatch(r'(.+?)(?:_S\d+)?(?:_(L\d{3}))?[_\.-]R([12])(?:[_\.-]001)?',name,re.I)
        if not match:unmatched.append(f.name);continue
        sample,lane,end=match.groups();key=(sample,lane or '')
        row=rows.setdefault(key,dict(sample=sample,group='',lane=lane or '',fastq1='',fastq2=''))
        if row['fastq'+end]:raise ValueError('同一sample/lane测序端不唯一：'+sample)
        row['fastq'+end]=str(f)
    values=list(rows.values())
    for row in values:
        if not row['fastq1']:raise ValueError('缺少R1：'+row['sample'])
    if not values:raise ValueError('未发现可识别的_R1/_R2命名FASTQ，请使用样本表')
    if len({bool(r['fastq2']) for r in values})>1:raise ValueError('目录混合单端/双端或缺失R2')
    return dict(samples=values,unmatched=unmatched,lanes_require_merge=any(sum(r['sample']==x['sample'] for x in values)>1 for r in values))


def save_sample_draft(root,rows,column_map=None):
    if not isinstance(rows,list) or not rows or len(rows)>5000:raise ValueError('样本草稿需要1–5000行')
    frame=pd.DataFrame(rows).fillna('').astype(str)
    if column_map:frame=frame.rename(columns=column_map)
    if frame.columns.duplicated().any():raise ValueError('列映射产生重复列名，请更换映射或移除重复列')
    if 'fastq2' in frame and (frame.fastq2.str.strip()=='').all():frame=frame.drop(columns='fastq2')
    if not {'sample','group'}.issubset(frame):raise ValueError('样本草稿需要sample/group列')
    if frame['sample'].duplicated().any():raise ValueError('同一生物样本跨lane需先合并；不能把lane当成独立重复')
    payload=frame.to_csv(index=False).encode('utf-8');key=hashlib.sha256(payload).hexdigest()[:20]
    folder=Path(root)/'_inputs';folder.mkdir(exist_ok=True)
    path=folder/(key+'-samples.csv');path.write_bytes(payload);return str(path)


def excel_export(p):
    """Full precision Excel with identifiers as strings, plus a methods/index sheet."""
    from results import TABLES,read_json
    from dataio import read_frame
    output=io.BytesIO();index=[];used={'索引','分析参数'}
    with pd.ExcelWriter(output,engine='openpyxl') as writer:
        targets=sorted((p/'contrasts').glob('*'));profiles=[q/'profiles'/kind for q in targets for kind in ('strict','nominal') if (q/'profiles'/kind).exists()]
        for root in [p]+targets+profiles:
            for name in sorted(TABLES):
                file=root/(name+'.csv')
                if not file.exists():continue
                frame=read_frame(file)
                if len(frame)>1048575:continue
                label=(root.name[:10]+'_' if root!=p else '')+name
                label=re.sub(r'[\\/*?:\[\]]','_',label)[:28];sheet=label;counter=1
                while sheet in used:counter+=1;sheet=label[:27]+str(counter)
                used.add(sheet);frame.to_excel(writer,sheet_name=sheet,index=False)
                index.append(dict(sheet=sheet,file=str(file.relative_to(p)),rows=len(frame)))
        pd.DataFrame(index).to_excel(writer,sheet_name='索引',index=False)
        cfg=read_json(p/'config.json',{});pd.DataFrame([dict(parameter=k,value=json.dumps(v,ensure_ascii=False)) for k,v in cfg.items()]).to_excel(writer,sheet_name='分析参数',index=False)
        for sheet in writer.book:
            sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
            for row in sheet:
                for cell in row:
                    # IDs remain text; do not turn arbitrary annotation strings into Excel formulas.
                    if isinstance(cell.value,str):cell.data_type='s'
    output.seek(0);return output


def merge_lanes(root,rows):
    """Concatenate gzip members without recompressing; preserves original lane files."""
    import shutil
    groups={}
    for row in rows:groups.setdefault(row['sample'],[]).append(row)
    merged=[]
    for sample,lanes in groups.items():
        row=dict(lanes[0]);row.pop('lane',None)
        if len(lanes)>1:
            files=[local_file(x[k]) for x in lanes for k in ('fastq1','fastq2') if x.get(k)]
            if any(not str(f).lower().endswith('.gz') for f in files):raise ValueError('自动跨lane合并需gzip FASTQ；原始文本请先压缩')
            signatures=[dict(path=str(f),size=f.stat().st_size,mtime=f.stat().st_mtime_ns) for f in files]
            key=hashlib.sha256(json.dumps(signatures,sort_keys=True).encode()).hexdigest()[:20]
            folder=Path(root)/'_inputs'/('lanes-'+key);folder.mkdir(parents=True,exist_ok=True)
            for end in ('fastq1','fastq2'):
                if not row.get(end):continue
                target=folder/(end+'.fastq.gz');tmp=folder/(end+'.tmp')
                if not target.exists():
                    with tmp.open('wb') as dst:
                        for lane in lanes:
                            with local_file(lane[end]).open('rb') as src:shutil.copyfileobj(src,dst,1024*1024)
                    tmp.replace(target)
                row[end]=str(target)
        merged.append(row)
    return merged
