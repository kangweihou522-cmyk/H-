"""Browser regression for reviewed canvas, preview, controls and advanced plots."""
from pathlib import Path
import argparse,json,urllib.request,shutil
from playwright.sync_api import sync_playwright
parser=argparse.ArgumentParser();parser.add_argument('--base',default='http://127.0.0.1:18880');parser.add_argument('--job',required=True);parser.add_argument('--output',type=Path,default=Path('runs/review-browser'));args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
base=args.base;job=args.job
manifest=json.load(urllib.request.urlopen(base+'/api/jobs/'+job+'/manifest'))
assert len(manifest['samples'])==18 and 'T_C' in manifest['contrasts'],'Use the synthetic 18-sample fixture'
errors=[]
with sync_playwright() as p:
 browser=p.chromium.launch(executable_path=shutil.which('chromium') or shutil.which('chromium-browser'),headless=True,args=['--no-sandbox']);page=browser.new_page(viewport={'width':1500,'height':1000});page.on('pageerror',lambda e:errors.append(str(e)))
 page.goto(base);page.locator('#project').select_option(job);page.wait_for_function("document.querySelector('#projectstatus').textContent.includes('已完成')")
 page.locator('#contrast').select_option('T_C');page.locator('[data-view="differential"]').click();page.wait_for_function("document.querySelector('#chartselect').value==='volcano|T_C' && !!document.getElementById('plot').data")
 page.wait_for_timeout(500)
 size=page.locator('#plot').bounding_box();assert 1<size['width']/size['height']<1.12,size
 page.locator('#editoradvanced').evaluate('(e)=>e.open=true');page.locator('[data-setting="annotated_genes"]').fill('Synthetic_0');page.locator('[data-setting="point_size"]').fill('70')
 original=urllib.request.urlopen(base+'/api/jobs/'+job+'/files/contrasts/T_C/volcano.settings.json').read()
 page.locator('#previewchart').click();page.wait_for_function("document.querySelector('#chartstatus').textContent.includes('临时预览')")
 assert urllib.request.urlopen(base+'/api/jobs/'+job+'/files/contrasts/T_C/volcano.settings.json').read()==original
 assert page.evaluate("document.getElementById('plot').layout.annotations.some(a=>a.name==='ENSMUSG00000000001')")
 page.locator('#editor').get_by_label('标题',exact=True).fill('Reviewed volcano');page.locator('#redraw').click();page.wait_for_function("document.querySelector('#chartstatus').textContent.includes('报告已更新')",timeout=60000)
 page.locator('#plot').scroll_into_view_if_needed();page.screenshot(path=str(args.output/'volcano-desktop.png'),full_page=True)
 # Emitted click is the same Plotly event used by a real point pick.
 page.evaluate("document.getElementById('plot').emit('plotly_click',{points:[{customdata:'ENSMUSG00000000001'}]})");page.wait_for_selector('#pointinspector strong');assert 'Synthetic_0' in page.locator('#pointinspector').inner_text()
 page.locator('[data-view="enrichment"]').click();page.locator('#chartselect').select_option('leading_edge|T_C');page.wait_for_function("document.querySelector('#chartselect').value==='leading_edge|T_C' && (document.getElementById('plot').data||[])[0]?.type==='heatmap'")
 page.locator('#plot').scroll_into_view_if_needed();page.screenshot(path=str(args.output/'leading-edge.png'),full_page=True)
 page.locator('#chartselect').select_option('enrichment_network|T_C');page.wait_for_function("document.querySelector('#chartnote').textContent.includes('Jaccard')")
 assert page.locator('[data-setting="network_threshold"]').count()==1
 page.locator('[data-view="mining"]').click();page.locator('#tablekind').select_option('candidate-evidence');page.wait_for_function("document.querySelector('#resulttable').textContent.includes('leading_edge_count')")
 page.locator('#tablekind').select_option('gene-trait-adjusted');page.wait_for_function("document.querySelector('#resulttable').textContent.includes('correlation')")
 page.locator('#chartselect').select_option('gene_trait|T_C');page.wait_for_function("document.querySelector('#chartnote').textContent.includes('偏 Pearson')")
 page.locator('#plot').scroll_into_view_if_needed();page.screenshot(path=str(args.output/'adjusted-trait.png'),full_page=True)
 page.locator('[data-view="network"]').click();page.locator('#chartselect').select_option('wgcna_network|');page.wait_for_function("document.querySelector('#chartnote').textContent.includes('signed adjacency')")
 page.locator('#plot').scroll_into_view_if_needed();page.screenshot(path=str(args.output/'network.png'),full_page=True)
 page.set_viewport_size({'width':390,'height':844});page.locator('[data-view="differential"]').click();page.wait_for_function("document.querySelector('#chartselect').value==='volcano|T_C'");page.wait_for_timeout(800)
 assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth+1')
 page.locator('#plot').scroll_into_view_if_needed();page.screenshot(path=str(args.output/'volcano-mobile.png'),full_page=True)
 assert not errors,errors
 browser.close()
print('PASS: square canvas, non-writing preview, gene labels/click, synchronized save, leading-edge/network/adjusted trait/candidate views, mobile overflow, no JS errors')
