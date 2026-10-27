"""Optional Chromium/Playwright end-to-end check on the synthetic fixture."""
from pathlib import Path
from playwright.sync_api import sync_playwright
import json,urllib.request,zipfile,io,argparse,shutil
parser=argparse.ArgumentParser(description='Only use the completed six-sample synthetic smoke job; this edits its plot/profile.')
parser.add_argument('--job',required=True)
parser.add_argument('--base',default='http://127.0.0.1:8765')
parser.add_argument('--chromium',default=shutil.which('chromium') or shutil.which('chromium-browser'))
parser.add_argument('--output',type=Path,default=Path('runs/browser-validation'))
args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
j=args.job;base=args.base;errors=[]
manifest=json.load(urllib.request.urlopen(base+'/api/jobs/'+j+'/manifest'))
assert len(manifest['samples'])==6 and 'comparison' in manifest['contrasts'], 'Use the six-sample synthetic smoke fixture'
probe=json.load(urllib.request.urlopen(base+'/api/jobs/'+j+'/gene?gene=Synthetic_0'))
assert probe['candidates'][0]['symbol']=='Synthetic_0'
req=urllib.request.Request(base+'/api/jobs/'+j+'/profile',data=json.dumps(dict(contrast='comparison',profile='strict')).encode(),headers={'Content-Type':'application/json'})
json.load(urllib.request.urlopen(req,timeout=90))
with sync_playwright() as p:
 b=p.chromium.launch(executable_path=args.chromium,headless=True,args=['--no-sandbox']);page=b.new_page(viewport={'width':1500,'height':1000});page.on('pageerror',lambda e:errors.append(str(e)))
 page.goto(base);page.locator('#project').select_option(j);page.wait_for_function("document.querySelector('#projectstatus').textContent.includes('已完成')")
 page.locator('[data-view="differential"]').click();page.wait_for_selector('#profilecomparison .profile-choice');assert page.locator('#profilecomparison .profile-choice').count()==2
 page.wait_for_function("!!document.getElementById('plot').data && document.getElementById('plot').data[0].customdata")
 # Exercise Plotly click/relayout handlers using the events emitted by the library.
 page.evaluate("document.getElementById('plot').emit('plotly_click',{points:[{customdata:document.getElementById('plot').data[0].customdata[0]}]})")
 page.wait_for_selector('#pointinspector strong');assert 'ENSMUSG' in page.locator('#pointinspector').inner_text()
 page.locator('#markgene').fill('Synthetic_0');page.locator('#markpoint').click();page.wait_for_function("document.getElementById('chartstatus').textContent.includes('基因已定位并标注') && (document.getElementById('plot').layout.annotations||[]).some(a=>a.name==='ENSMUSG00000000001')")
 page.evaluate("async()=>{const el=document.getElementById('plot');const i=el.layout.annotations.findIndex(a=>a.name==='ENSMUSG00000000001');await Plotly.relayout(el,{['annotations['+i+'].ax']:70})}")
 page.locator('#editor').get_by_label('标题',exact=True).fill('Advanced labelled volcano');page.locator('#redraw').click();page.wait_for_function("document.getElementById('chartstatus').textContent.includes('报告已更新')",timeout=30000)
 cfg=json.load(urllib.request.urlopen(base+'/api/jobs/'+j+'/files/contrasts/comparison/volcano.settings.json'));assert cfg['label_positions']['ENSMUSG00000000001']['ax']==70
 svg=urllib.request.urlopen(base+'/api/jobs/'+j+'/files/contrasts/comparison/volcano.svg').read().decode();assert 'Synthetic_0' in svg and 'Advanced labelled volcano' in svg
 before=json.load(urllib.request.urlopen(base+'/api/jobs/'+j+'/curve?contrast=comparison&database=pathway&term=synthetic_up'))
 page.locator('#profilecomparison .profile-choice').filter(has_text='方案一').get_by_role('button').click();page.wait_for_function("document.querySelector('#profilecomparison .profile-choice.selected').textContent.includes('方案一')",timeout=30000)
 after=json.load(urllib.request.urlopen(base+'/api/jobs/'+j+'/curve?contrast=comparison&database=pathway&term=synthetic_up'));assert before['running']==after['running']
 manifest=json.load(urllib.request.urlopen(base+'/api/jobs/'+j+'/manifest'));c=manifest['contrasts']['comparison'];assert c['active_profile']=='nominal' and c['profiles']['nominal']['total']>=c['profiles']['strict']['total']
 page.wait_for_function("document.getElementById('plot').layout.yaxis.title.text.includes('unadjusted')",timeout=30000)
 assert page.evaluate("document.getElementById('plot').layout.yaxis.title.text.includes('unadjusted')")
 page.screenshot(path=str(args.output/'volcano.png'),full_page=True)
 page.locator('[data-view="enrichment"]').click();page.locator('#tablekind').select_option('pathway_gsea');page.wait_for_selector('#resulttable tbody tr td button');page.locator('#resulttable tbody tr td button').first.click();page.wait_for_function("document.querySelector('#chartselect').value.startsWith('gsea_curve')")
 page.wait_for_function("(document.getElementById('plot').data||[]).some(d=>d.name==='Gene-set hits')")
 page.locator('[data-view="report"]').click()
 with page.expect_download() as d:page.get_by_text('下载完整结果 ZIP',exact=True).click()
 d.value.save_as(str(args.output/'results.zip'))
 with zipfile.ZipFile(str(args.output/'results.zip')) as z:assert z.testzip() is None and 'report.html' in z.namelist()
 assert not errors,errors;b.close()
 print('PASS: point click handler, Symbol location, label relayout events/persistence, SVG, dual profile switch, unchanged GSEA, curve drilldown, ZIP, no JS errors')
