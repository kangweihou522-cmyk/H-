"""Assert end-to-end outputs and known synthetic fold-change directions."""
from pathlib import Path
import pandas as pd
import json, sys
p=Path(sys.argv[1]);r=pd.read_csv(p/'differential.csv');assert len(r)==240
up=r[r.gene.isin([f'ENSMUSG{i:011d}' for i in range(1,41)])];down=r[r.gene.isin([f'ENSMUSG{i:011d}' for i in range(41,81)])]
assert (up.log2FoldChange>1).mean()>.9;assert (down.log2FoldChange<-1).mean()>.9
assert (up.padj<.05).mean()>.8;assert (down.padj<.05).mean()>.8
for sample in ['C1','C2','C3','T1','T2','T3']:
 q=pd.read_csv(p/'quant'/sample/'quant.sf',sep='\t');assert len(q)==240 and q.NumReads.sum()>10000
 meta=json.loads((p/'quant'/sample/'aux_info/meta_info.json').read_text());assert meta['percent_mapped']>99
for kind in ['go','gsea','pathway','pathway_gsea']:
 a=pd.read_csv(p/(kind+'.csv'));assert len(a)>0 and a.padj.between(0,1).all()
g=pd.read_csv(p/'pathway_gsea.csv').set_index('term');assert g.loc['synthetic_up','NES']>0 and g.loc['synthetic_down','NES']<0
for kind in ['pca','volcano','heatmap','go','gsea','pathway','pathway_gsea']:
 for fmt in ['svg','png','pdf','plot-data.csv']:
  assert (p/(kind+'.'+fmt)).stat().st_size>100
print(json.dumps({'genes_tested':len(r),'up_significant':int((up.padj<.05).sum()),'down_significant':int((down.padj<.05).sum()),'samples_quantified':6,'charts':7,'formats':['SVG','PNG','PDF'],'GO_terms':len(pd.read_csv(p/'go.csv')),'pathway_sets':len(g),'GSEA_up_NES':float(g.loc['synthetic_up','NES']),'GSEA_down_NES':float(g.loc['synthetic_down','NES'])},indent=2))
