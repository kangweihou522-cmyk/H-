"""Copy verified completed FASTQ stages into a new, independently recorded run."""
from pathlib import Path
import json
import shutil
import uuid


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None


def fingerprint(path):
    path=Path(path);stat=path.stat()
    return dict(path=str(path.resolve()),size=stat.st_size,mtime_ns=stat.st_mtime_ns)


def reference_signature(cfg):
    return dict(transcripts=fingerprint(cfg['transcripts']),genome=fingerprint(cfg['genome']) if cfg.get('genome') else None)


def reuse_completed(source,target,cfg,samples,check_cancel=lambda:None):
    """Never alter the source run or reuse statistics under a changed design.

    Complete markers are input path/size/mtime signatures, as in the original
    executor. FASTQ processing and quantification parameters are checked in
    addition. Files are copied, with cancellation checks between chunks, so a
    subsequent retry cannot write through a link into the earlier experiment.
    """
    source,target=Path(source),Path(target)
    report=dict(source_job=source.name,reused=dict(qc=[],index=False,quant=[]),copied_bytes=0,notes=[])
    old=read_json(source/'config.json') or {}
    provenance=read_json(source/'input-provenance.json') or {}
    report['source_commit']=provenance.get('source_commit')
    report['changed_parameters']=[key for key in sorted(set(old)|set(cfg)) if key not in ('resume_from','reuse_completed') and old.get(key)!=cfg.get(key)]
    if not cfg.get('reuse_completed',True):
        report['notes'].append('用户选择重新运行全部步骤。');return report
    if cfg['mode']!='fastq' or old.get('mode')!='fastq':
        report['notes'].append('当前为计数/已有定量导入，使用本次导入输入重新执行统计。');return report
    if (read_json(source/'status.json') or {}).get('state') not in ('failed','cancelled'):
        report['notes'].append('来源任务已重新启动，未复制其缓存。');return report
    same_preprocess=old.get('preprocess','fastp')==cfg.get('preprocess','fastp')
    same_qc=same_preprocess and all(old.get(key,default)==cfg.get(key,default) for key,default in [('quality',20),('length',30)])
    refs=reference_signature(cfg)
    old_refs=read_json(source/'index-complete.json')
    # Legacy transcript-only runs recorded a single FASTA signature.
    if isinstance(old_refs,dict) and 'path' in old_refs:
        old_refs=dict(transcripts=old_refs,genome=None)
    same_index=old_refs==refs and (source/'index/versionInfo.json').is_file()
    quant_defaults=dict(library='A',fragment_mean=200,fragment_sd=80)
    same_quant=same_preprocess and same_index and all(old.get(k,v)==cfg.get(k,v) for k,v in quant_defaults.items())
    if cfg.get('preprocess','fastp')=='fastp':same_quant=same_quant and same_qc
    if not same_qc:report['notes'].append('预处理方式或 fastp 参数改变，相关质控和定量重新运行。')
    if not same_index:report['notes'].append('索引未完成或参考签名改变，重新建立索引并定量。')
    elif not same_quant:report['notes'].append('定量或预处理参数改变，重新执行 Salmon 定量。')
    staging=target/'.resume-cache'

    def copy_directory(original,destination):
        check_cancel()
        if destination.exists():return False
        required=sum(f.stat().st_size for f in original.rglob('*') if f.is_file())
        if shutil.disk_usage(target).free < required+16*1024**2:
            raise ValueError('续跑缓存复制空间不足；可关闭“复用已完成步骤”后新建任务，或释放输出盘空间。')
        stage=staging/uuid.uuid4().hex;stage.parent.mkdir(exist_ok=True)
        def copy_file(src,dst):
            with open(src,'rb') as incoming,open(dst,'wb') as outgoing:
                while True:
                    check_cancel();chunk=incoming.read(1024**2)
                    if not chunk:break
                    outgoing.write(chunk)
            shutil.copystat(src,dst)
            return dst
        try:
            shutil.copytree(original,stage,copy_function=copy_file)
            check_cancel();destination.parent.mkdir(parents=True,exist_ok=True);stage.rename(destination)
            report['copied_bytes']+=required
            return True
        finally:
            if stage.exists():shutil.rmtree(stage)

    try:
        if same_index and copy_directory(source/'index',target/'index'):
            (target/'index-complete.json').write_text(json.dumps(refs))
            report['reused']['index']=True
        for _,sample in samples.iterrows():
            check_cancel();name=sample['sample']
            files=[sample.fastq1]+([sample.fastq2] if 'fastq2' in samples else [])
            signature=[fingerprint(path) for path in files]
            qc=source/'qc'/name;quant=source/'quant'/name
            clean_files=['R1.fastq.gz','fastp.json']+(['R2.fastq.gz'] if len(files)==2 else [])
            qc_ok=same_qc and read_json(qc/'complete.json')==signature and all((qc/f).is_file() and (qc/f).stat().st_size for f in clean_files)
            if cfg.get('preprocess','fastp')=='fastp' and qc_ok:
                if copy_directory(qc,target/'qc'/name):report['reused']['qc'].append(name)
            quant_ok=same_quant and read_json(quant/'complete.json')==signature and (quant/'quant.sf').is_file() and (quant/'quant.sf').stat().st_size>0 and read_json(quant/'aux_info/meta_info.json') is not None
            if cfg.get('preprocess','fastp')=='fastp':quant_ok=quant_ok and qc_ok
            if quant_ok and copy_directory(quant,target/'quant'/name):report['reused']['quant'].append(name)
    finally:
        if staging.exists():shutil.rmtree(staging)
    report['notes'].append('统计、富集、挖掘和报告按本次参数重新计算；原任务保留。缓存依据输入路径/大小/修改时间与参数校验，不是完整文件哈希校验。')
    return report
