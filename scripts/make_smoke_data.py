"""Local FASTQ + counts + synthetic GO/GMT fixture; no remote reads required."""
from pathlib import Path
import gzip, csv, json
import numpy as np
import pandas as pd

def create(p):
 p=Path(p);p.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(1234);n=240
 genes=[f'ENSMUSG{i:011d}' for i in range(1,n+1)];tx=[f'ENSMUST{i:011d}' for i in range(1,n+1)];seq=[''.join(rng.choice(list('ACGT'),1000)) for _ in range(n)]
 (p/'transcripts.fa').write_text(''.join(f'>{t}\n{s}\n' for t,s in zip(tx,seq)))
 pd.DataFrame({'transcript':tx,'gene':genes}).to_csv(p/'tx2gene.tsv',sep='\t',index=False)
 samples=['C1','C2','C3','T1','T2','T3'];counts=[];rows=[];rc=str.maketrans('ACGT','TGCA')
 for si,sample in enumerate(samples):
  mu=np.full(n,80.);mu[:40]*=4 if si>=3 else 1;mu[40:80]*=.25 if si>=3 else 1
  c=rng.negative_binomial(30,30/(30+mu));counts.append(c)
  f1=p/(sample+'_R1.fastq.gz');f2=p/(sample+'_R2.fastq.gz')
  with gzip.open(f1,'wt') as a,gzip.open(f2,'wt') as b:
   for j,total in enumerate(c):
    for k in range(total):
     start=int(rng.integers(0,650));r1=seq[j][start:start+100];r2=seq[j][start+150:start+250].translate(rc)[::-1]
     a.write(f'@{sample}_{j}_{k}/1\n{r1}\n+\n'+('I'*100)+'\n');b.write(f'@{sample}_{j}_{k}/2\n{r2}\n+\n'+('I'*100)+'\n')
  rows.append(dict(sample=sample,group='Control' if si<3 else 'Treatment',fastq1=str(f1),fastq2=str(f2)))
 pd.DataFrame(rows).to_csv(p/'samples.csv',index=False);pd.DataFrame(np.array(counts).T,index=genes,columns=samples).rename_axis('gene').to_csv(p/'counts.csv')
 # Artificial memberships deliberately labelled as a fixture, not real mouse biology.
 pd.DataFrame([{'gene':g,'symbol':'Synthetic_'+str(i),'entrez':str(i+1),'go':'GO:0006955' if i<40 else 'GO:0008152','ontology':'BP'} for i,g in enumerate(genes)]).to_csv(p/'annotation.tsv',sep='\t',index=False)
 (p/'pathways.gmt').write_text('synthetic_up\tSynthetic induced genes\t'+'\t'.join(genes[:40])+'\nsynthetic_down\tSynthetic suppressed genes\t'+'\t'.join(genes[40:80])+'\nsynthetic_other\tSynthetic unchanged genes\t'+'\t'.join(genes[80:])+'\n')
 cfg=dict(mode='fastq',samples=str(p/'samples.csv'),transcripts=str(p/'transcripts.fa'),tx2gene=str(p/'tx2gene.tsv'),annotation=str(p/'annotation.tsv'),gmt=str(p/'pathways.gmt'),gmt_id='ENSEMBL',control='Control',treatment='Treatment',threads=2,alpha=.05,lfc=1,permutations=100,library='IU')
 (p/'request.json').write_text(json.dumps(cfg,indent=2));return cfg
if __name__=='__main__':
 import sys
 create(Path(sys.argv[1]).resolve());print('Created local synthetic FASTQ fixture')
