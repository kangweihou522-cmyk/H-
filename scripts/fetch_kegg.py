"""Download mouse KEGG pathways via public REST; respect KEGG terms of use."""
import urllib.request, hashlib, json, datetime, argparse
from pathlib import Path

def fetch(url):
 with urllib.request.urlopen(url,timeout=60) as r:return r.read().decode()
def main():
 parser=argparse.ArgumentParser(description='建立 KEGG mouse GMT（Entrez ID）');parser.add_argument('output',type=Path);args=parser.parse_args()
 names=dict(line.split('\t',1) for line in fetch('https://rest.kegg.jp/list/pathway/mmu').splitlines());sets={}
 for line in fetch('https://rest.kegg.jp/link/pathway/mmu').splitlines():
  gene,path=line.split('\t');sets.setdefault(path.removeprefix('path:'),set()).add(gene.removeprefix('mmu:'))
 text='\n'.join('\t'.join([p,names.get(p,p),*sorted(g)]) for p,g in sorted(sets.items()))+'\n'
 args.output.write_text(text)
 args.output.with_suffix('.source.json').write_text(json.dumps({'organism':'mmu','id_type':'ENTREZID','source':'https://rest.kegg.jp','retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(text.encode()).hexdigest()},indent=2))
 print(f'{len(sets)} pathways -> {args.output}; 界面中选择 ENTREZID')
if __name__=='__main__':main()
