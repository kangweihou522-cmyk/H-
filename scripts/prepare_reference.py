"""Build tx2gene from an Ensembl GTF (same release as the transcript FASTA)."""
import argparse,gzip,re,csv
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('gtf',type=Path);p.add_argument('output',type=Path);a=p.parse_args();pairs={}
with (gzip.open(a.gtf,'rt') if a.gtf.suffix=='.gz' else a.gtf.open()) as f:
 for line in f:
  if line.startswith('#'):continue
  fields=line.rstrip().split('\t')
  if len(fields)!=9:continue
  attrs=dict(re.findall(r'(\w+) "([^"]+)"',fields[8]))
  if 'transcript_id' in attrs and 'gene_id' in attrs:
   t=attrs['transcript_id'];g=attrs['gene_id']
   if t in pairs and pairs[t]!=g:raise ValueError('ambiguous transcript '+t)
   pairs[t]=g
if not pairs:raise ValueError('GTF 没有 transcript_id / gene_id')
with a.output.open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['transcript','gene']);w.writerows(sorted(pairs.items()))
print(f'{len(pairs)} transcript-gene mappings')
