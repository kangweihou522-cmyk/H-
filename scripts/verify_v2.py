"""Independently check a completed v2 run's profiles, GSEA ranks and mining outputs."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def verify(p):
    p = Path(p)
    assert json.loads((p/'status.json').read_text())['state'] == 'completed'
    summary = json.loads((p/'analysis-summary.json').read_text())
    samples = pd.read_csv(p/'samples.csv')['sample'].tolist()
    for name, contrast in summary['contrasts'].items():
        q = p/'contrasts'/name
        for profile, column, lfc in [('nominal', 'pvalue', 0), ('strict', 'padj', 1)]:
            s = q/'profiles'/profile
            table = pd.read_csv(s/'differential.csv')
            significant = (table[column] < .05) & (table.log2FoldChange.abs() > lfc)
            selected = set(table.loc[significant, 'gene'])
            assert set(pd.read_csv(s/'deg-all.csv').gene) == selected
            assert set(pd.read_csv(s/'deg-up.csv').gene) == set(table.loc[significant & (table.log2FoldChange > 0), 'gene'])
            assert set(pd.read_csv(s/'deg-down.csv').gene) == set(table.loc[significant & (table.log2FoldChange < 0), 'gene'])
            assert contrast['profiles'][profile]['total'] == len(selected)
        for prefix in ['gsea', 'pathway_gsea']:
            file=q/(prefix+'-ranks.csv')
            if not file.exists():continue
            ranks = pd.read_csv(file)
            assert ranks.gene.is_unique and np.isfinite(ranks.stat).all()
            assert np.all(np.diff(ranks.stat.to_numpy()) <= 0)
        print(name, {key: spec['total'] for key, spec in contrast['profiles'].items()})
    mining = json.loads((p/'mining-status.json').read_text())
    for key, file in [('pathway_activity', 'pathway-activity.csv'), ('regulatory_activity', 'tf-activity.csv')]:
        if mining.get(key, {}).get('state') == 'completed':
            values = pd.read_csv(p/file, index_col=0)
            assert list(values.columns) == samples and np.isfinite(values.to_numpy()).all()
    if mining['wgcna']['state'] == 'completed':
        eig = pd.read_csv(p/'wgcna-eigengenes.csv', index_col=0)
        traits = pd.read_csv(p/'traits.csv', index_col=0).loc[eig.index]
        associations = pd.read_csv(p/'wgcna-module-trait.csv')
        for row in associations.itertuples():
            correlation, pvalue = pearsonr(eig[row.module], traits[row.trait])
            assert np.isclose(row.correlation, correlation, atol=1e-12)
            assert np.isclose(row.pvalue, pvalue, atol=1e-12)
        assert (p/'wgcna-soft-threshold.csv').is_file()
        edges = pd.read_csv(p/'wgcna-network.csv')
        assert len(edges) <= 500
    assert '双方案' in (p/'report.html').read_text() or 'nominal' in (p/'report.html').read_text()
    print('PASS: exact profile conditions/DEG directions, sorted finite GSEA ranks, completed activity matrices, independent module–trait correlations, report')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory')
    verify(parser.parse_args().run_directory)
