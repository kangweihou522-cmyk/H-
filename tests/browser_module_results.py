"""Standalone completed projects open their own results, with truthful inputs."""
import json,shutil,argparse
from pathlib import Path
from playwright.sync_api import sync_playwright
parser=argparse.ArgumentParser();parser.add_argument('--base',default='http://127.0.0.1:18882');parser.add_argument('--ids',type=Path,default=Path('/tmp/rnaseq-modules-validation/job-ids.json'));args=parser.parse_args();ids=json.loads(args.ids.read_text());errors=[]
with sync_playwright() as p:
 b=p.chromium.launch(executable_path=shutil.which('chromium'),headless=True,args=['--no-sandbox']);page=b.new_page(viewport={'width':1440,'height':1000});page.on('pageerror',lambda e:errors.append(str(e)));page.goto(args.base)
 for key,view in [('enrichment','enrichment'),('activity_final','mining'),('network','network'),('log_expression','expression'),('differential','differential')]:
  page.locator('#project').select_option(ids[key]);page.wait_for_function("(id)=>state.job===id&&!!state.manifest",arg=ids[key]);page.wait_for_function('(view)=>state.view===view',arg=view);page.wait_for_function("!!document.getElementById('plot').data")
  assert not page.locator('#chartpanel').is_hidden()
 page.locator('#project').select_option(ids['log_expression']);page.wait_for_function("(id)=>state.job===id&&!!state.manifest&&state.view==='expression'",arg=ids['log_expression'])
 page.locator('#chartselect').select_option('distribution|');page.wait_for_function("document.getElementById('plot').layout.yaxis.title.text.includes('log2(TPM + 1)')")
 page.locator('[data-view="qc"]').click();assert '单步分析' in page.locator('#qcsummary').inner_text() and '已处理 FASTQ' not in page.locator('#qcsummary').inner_text()
 page.locator('[data-view="mining"]').click();page.locator('#genequery').fill('Synthetic_0');page.locator('#searchgene').click();page.wait_for_function("!!document.getElementById('geneplot').data");assert 'log2(TPM + 1)' in page.evaluate("document.getElementById('geneplot').layout.yaxis.title.text")
 page.locator('#project').select_option(ids['ora_only']);page.wait_for_function("(id)=>state.job===id&&!!state.manifest&&state.view==='enrichment'",arg=ids['ora_only'])
 assert page.locator('#chartselect option[value^="gsea_curve"]').count()==0
 page.locator('#chartselect').select_option('volcano|comparison');page.wait_for_function("state.chart?.kind==='volcano'&&document.getElementById('plot').data[0].customdata?.[0]?.startsWith('ENSMUSG')&&document.getElementById('plot')._ev?.listenerCount('plotly_click')>0")
 page.evaluate("document.getElementById('plot').emit('plotly_click',{points:[{customdata:document.getElementById('plot').data[0].customdata[0]}]})");page.wait_for_selector('#pointinspector strong');assert page.locator('#pointinspector button').filter(has_text='未导入表达矩阵').is_disabled()
 assert not errors,errors;b.close()
print('PASS: standalone project navigation, correct matrix units, no false FASTQ QC, gene expression without DE fit, no fabricated GSEA or sample expression')
