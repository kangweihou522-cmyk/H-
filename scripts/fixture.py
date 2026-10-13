"""Deterministic synthetic negative-binomial counts, explicitly not biological evidence."""
from pathlib import Path
import numpy as np
import pandas as pd

def generate(p):
 rng=np.random.default_rng(42)
 # Mouse Ensembl IDs supplied from installed org.Mm.eg.db in the end-to-end fixture.
 ids_file=Path(p)/'mouse-ids.txt'
 ids=ids_file.read_text().splitlines()[:800] if ids_file.exists() else [f'ENSMUSG{i:011d}' for i in range(1,801)]
 names=['C1','C2','C3','T1','T2','T3'];mu=rng.uniform(50,500,(len(ids),6));mu[:80,3:]*=4;mu[80:160,3:]*=.25
 counts=rng.negative_binomial(20,20/(20+mu));pd.DataFrame(counts,index=ids,columns=names).rename_axis('gene').to_csv(Path(p)/'counts.csv')
 pd.DataFrame({'sample':names,'group':['Control']*3+['Treatment']*3}).to_csv(Path(p)/'samples.csv',index=False)
if __name__=='__main__':
 import sys
 generate(Path(sys.argv[1]))
