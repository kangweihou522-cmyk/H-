"""Reference-assisted RNA mixture exploration using non-negative least squares."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from validation import read_table
from dataio import read_frame


def estimate(expression,reference):
    for x in (expression,reference):
        if x.index.has_duplicates or not np.isfinite(x.to_numpy(float)).all() or (x.to_numpy(float)<0).any() or (x.sum(axis=0)<=0).any():raise ValueError('组成矩阵须为唯一基因、非负有限线性TPM，列总量大于0')
    genes=expression.index.intersection(reference.index)
    if len(genes)<50 or reference.shape[1]<2:raise ValueError('至少50个匹配基因及2种参考细胞类型')
    expr=expression.div(expression.sum(axis=0),axis=1)*1e6;ref=reference.div(reference.sum(axis=0),axis=1)*1e6
    A=ref.loc[genes].to_numpy(float);norm=np.linalg.norm(A,axis=0)
    if np.any(norm==0) or np.linalg.matrix_rank(A)<A.shape[1] or np.linalg.cond(A/norm)>1e6:raise ValueError('参考签名无法区分细胞类型：零覆盖、秩不足或严重共线')
    fractions=[];diagnostics=[]
    for sample in expr:
        b=expr.loc[genes,sample].to_numpy(float);weights,residual=nnls(A,b);total=weights.sum()
        if total<=0:raise ValueError('无法拟合非零参考贡献：'+str(sample))
        fitted=A@weights;den=np.linalg.norm(fitted)*np.linalg.norm(b)
        fractions.append(dict(sample=sample,**dict(zip(ref.columns,weights/total))))
        diagnostics.append(dict(sample=sample,matched_genes=len(genes),input_genes=len(expr),reference_genes=len(ref),coefficient_sum=float(total),rmse=float(np.sqrt(np.mean((b-fitted)**2))),cosine=float(np.dot(b,fitted)/den) if den else np.nan))
    return pd.DataFrame(fractions),pd.DataFrame(diagnostics)


def run(p):
    p=Path(p);cfg=json.loads((p/'config.json').read_text());status=dict(state='skipped',reason='未提供匹配组织/物种的线性TPM参考签名')
    if not cfg.get('composition_reference'):return status
    try:
        if (p/'tpm.csv').exists():expr=read_frame(p/'tpm.csv',index_col=0)
        elif cfg.get('expression_unit')=='log2_tpm':expr=np.expm1(read_frame(p/'vst.csv',index_col=0)*np.log(2))
        else:raise ValueError('需Salmon基因TPM或明确的log2(TPM+1)；counts/VST不能当作TPM')
        ref=read_table(p/'composition-reference.csv',index_col=0,converters={0:str});fractions,diagnostics=estimate(expr,ref)
        fractions.to_csv(p/'composition-fractions.csv',index=False);diagnostics.to_csv(p/'composition-diagnostics.csv',index=False)
        status=dict(state='completed',method='SciPy NNLS; exact gene IDs; column-normalized linear TPM; conditional RNA mixture contributions, not calibrated cell counts or MuSiC',reference_source=cfg.get('composition_source','unspecified'),reference_condition_number=float(np.linalg.cond(ref.loc[expr.index.intersection(ref.index)].to_numpy(float))),caution='依赖参考组织/状态、缺失细胞类型与RNA含量；不是实际细胞比例，也不提供未经验证的P值')
    except (ValueError,KeyError,TypeError) as e:status=dict(state='failed',reason=str(e))
    (p/'composition-status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));return status
