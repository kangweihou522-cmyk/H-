"""Local, artificial 18-sample fixture for multi-group/pathway/TF/WGCNA validation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd


def create(directory):
    p = Path(directory).resolve()
    p.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2026)
    n = 18
    genes = [f'ENSMUSG{i:011d}' for i in range(1, 601)]
    groups = ['Control']*6 + ['Treatment']*6 + ['Rescue']*6
    trait = np.array([-1]*6 + [1]*6 + [0]*6) + rng.normal(0, .1, n)
    latent = np.stack([trait, rng.normal(size=n), rng.normal(size=n)])
    mu = np.stack([np.exp(4.7 + .55*latent[i//200])*rng.uniform(.7, 1.3) for i in range(600)])
    counts = rng.negative_binomial(60, 60/(60+mu))
    metadata = pd.DataFrame(dict(sample=[f'S{i+1}' for i in range(n)], group=groups,
                                 batch=['B1', 'B2']*9, disease_score=trait))
    metadata.to_csv(p/'samples.csv', index=False)
    pd.DataFrame(counts, index=genes, columns=metadata['sample']).rename_axis('gene').to_csv(p/'counts.csv')
    # Artificial GO memberships and regulatory edges are NOT mouse biology.
    pd.DataFrame([dict(gene=g, symbol=f'Synthetic_{i}', entrez=str(i+1),
                       go='GO:0006955' if i < 200 else 'GO:0008152', ontology='BP')
                  for i, g in enumerate(genes)]).to_csv(p/'annotation.tsv', sep='\t', index=False)
    (p/'pathways.gmt').write_text('\n'.join(f'module_{i}\tSynthetic program {i}\t'+
                                          '\t'.join(genes[i*200:(i+1)*200]) for i in range(3))+'\n')
    pd.DataFrame([dict(source=f'Synthetic_TF{i//200+1}', target=g, weight=1)
                  for i, g in enumerate(genes)]).to_csv(p/'regulatory.tsv', sep='\t', index=False)
    cfg = dict(name='合成 18 样本 · 全方向验证', mode='counts', samples=str(p/'samples.csv'),
               counts=str(p/'counts.csv'), annotation=str(p/'annotation.tsv'), gmt=str(p/'pathways.gmt'),
               tf_network=str(p/'regulatory.tsv'), tf_id='ENSEMBL', gmt_id='ENSEMBL',
               comparisons=[dict(name='T_C', control='Control', treatment='Treatment'),
                            dict(name='R_T', control='Treatment', treatment='Rescue')],
               covariates=[dict(name='batch', type='categorical')], traits=['disease_score'],
               wgcna=True, wgcna_genes=600, wgcna_power=6, pathway_activity=True, min_set=15,
               max_set=500, permutations=100, threads=2, shrinkage='normal', deg_profile='strict')
    (p/'request.json').write_text(json.dumps(cfg, ensure_ascii=False, indent=2))
    return cfg


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    args = parser.parse_args()
    cfg = create(args.directory)
    print('Created synthetic input and request.json:', Path(args.directory).resolve())
