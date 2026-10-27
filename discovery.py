"""Exploratory evidence integration; no invented biomarker confidence scores."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from dataio import read_frame
from scipy import stats


def bh(values):
    """BH over all finite tests in the supplied family, preserving missing tests."""
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    order = valid[np.argsort(values[valid], kind='stable')]
    if len(order):
        result[order] = np.minimum(1, np.minimum.accumulate((values[order]*len(order)/np.arange(1, len(order)+1))[::-1])[::-1])
    return result


def adjusted_traits(p):
    cfg = json.loads((p/'config.json').read_text())
    samples = read_frame(p/'samples.csv', dtype=str)
    selected = cfg.get('traits', [])
    info = dict(state='skipped', reason='没有选定可检验的连续表型', skipped=[])
    if 'subject' in samples and samples.subject.nunique() < len(samples):
        info['reason'] = '存在重复 subject；当前独立样本相关检验不适用于重复测量'
        return info
    if len(samples) < 8:
        info['reason'] = '调整后基因–表型探索至少需要 8 个独立样本；较小样本仍可查看表达'
        return info
    design_file = p/'design-matrix.csv'
    if not design_file.exists():
        info['reason'] = '旧任务没有保存实际模型矩阵；需新建任务后计算调整后关联'
        return info
    design = read_frame(design_file, index_col=0).loc[samples['sample']].to_numpy(float)
    rank = np.linalg.matrix_rank(design)
    df = len(samples)-rank-1
    if df < 3:
        info['reason'] = '加入分组、协变量与表型后，残差自由度小于 3'
        return info
    projection = np.eye(len(samples))-design@np.linalg.pinv(design)
    expression = read_frame(p/'vst.csv', index_col=0)[samples['sample']]
    residuals = expression.to_numpy(float)@projection
    norms = np.linalg.norm(residuals, axis=1)
    frames = []
    for name in selected:
        trait = pd.to_numeric(samples[name], errors='coerce')
        if not trait.notna().all() or trait.nunique() < 3:
            info['skipped'].append(dict(trait=name, reason='仅检验至少三个取值的完整连续数值表型'))
            continue
        values = trait.to_numpy(float)
        residual = projection@values
        norm = np.linalg.norm(residual)
        if norm <= max(1, np.linalg.norm(values))*1e-10:
            info['skipped'].append(dict(trait=name, reason='表型被分组/协变量完全解释，不能独立估计关联'))
            continue
        correlation = np.divide(residuals@residual, norms*norm, out=np.full(len(norms), np.nan), where=norms>1e-10)
        correlation = np.clip(correlation, -1, 1)
        with np.errstate(divide='ignore', invalid='ignore'):
            statistic = np.abs(correlation)*np.sqrt(df/np.maximum(1-correlation**2, 1e-30))
        frames.append(pd.DataFrame(dict(gene=expression.index, trait=name, correlation=correlation,
                                       pvalue=2*stats.t.sf(statistic, df), n=len(samples), df=df)))
    if not frames:
        return info
    result = pd.concat(frames, ignore_index=True)
    result['padj'] = bh(result.pvalue)
    result.to_csv(p/'gene-trait-adjusted.csv', index=False)
    info.update(state='completed', reason='', tests=int(result.pvalue.notna().sum()),
                method='Partial Pearson on declared transformed expression; residualize gene and trait against the full DESeq2 design (including group); t test; BH across all finite gene × selected-trait tests',
                caution='独立样本线性模型探索；检查线性关系、离群点和残差，不能据此推断因果')
    return info


def build_candidates(p):
    cfg = json.loads((p/'config.json').read_text())
    ann = read_frame(p/'gene-annotation.csv', dtype=str).fillna('') if (p/'gene-annotation.csv').exists() else None
    module = read_frame(p/'wgcna-modules.csv') if (p/'wgcna-modules.csv').exists() else None
    trait = read_frame(p/'gene-trait-adjusted.csv') if (p/'gene-trait-adjusted.csv').exists() else None
    for q in sorted((p/'contrasts').glob('*')):
        if not (q/'differential.csv').exists():
            continue
        frame = read_frame(q/'differential.csv')
        counts = {gene: set() for gene in frame.gene}
        for filename, key in [('gsea.csv', 'id'), ('pathway_gsea.csv', {'ENSEMBL':'id','SYMBOL':'symbol','ENTREZID':'entrez'}.get(cfg.get('gmt_id'), 'id'))]:
            if not (q/filename).exists() or ann is None:
                continue
            enrichment = read_frame(q/filename)
            mapping = ann.groupby(key).gene.apply(list).to_dict()
            for row in enrichment[enrichment.padj < .05].itertuples():
                if pd.isna(row.leadingEdge):
                    continue
                for identifier in str(row.leadingEdge).split(';'):
                    for gene in mapping.get(identifier, []):
                        if gene in counts:
                            counts[gene].add(filename+':'+row.term)
        frame['leading_edge_sets'] = frame.gene.map(lambda g: ';'.join(sorted(counts[g])))
        frame['leading_edge_count'] = frame.gene.map(lambda g: len(counts[g]))
        if module is not None:
            frame = frame.merge(module[['gene', 'module', 'kME']], how='left', on='gene', validate='one_to_one')
        if trait is not None and not trait.empty:
            best = trait.dropna(subset=['padj']).sort_values(['padj', 'trait']).drop_duplicates('gene')
            best = best[['gene', 'trait', 'correlation', 'padj']].rename(columns={'trait':'best_trait', 'correlation':'trait_partial_r', 'padj':'trait_global_padj'})
            frame = frame.merge(best, how='left', on='gene', validate='one_to_one')
        frame.to_csv(q/'candidate-evidence.csv', index=False)


def run_discovery(p):
    p = Path(p)
    info = adjusted_traits(p)
    build_candidates(p)
    status_file = p/'mining-status.json'
    status = json.loads(status_file.read_text()) if status_file.exists() else {}
    status['gene_trait'] = info
    status['candidate_evidence'] = dict(state='completed', method='DESeq2 + significant GSEA leading-edge membership + optional WGCNA kME and adjusted trait association; evidence columns, no composite validation score')
    if not (p/'contrasts').is_dir() or not any((p/'contrasts').glob('*/differential.csv')):status['candidate_evidence']=dict(state='skipped',reason='此单步任务未拟合/导入差异基因结果')
    status_file.write_text(json.dumps(status, ensure_ascii=False, indent=2))
