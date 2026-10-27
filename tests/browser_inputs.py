"""Fresh forms, explicit recovery and upstream handoff without running fake jobs."""
import argparse,json,shutil
from playwright.sync_api import sync_playwright
p=argparse.ArgumentParser();p.add_argument('--base',default='http://127.0.0.1:18882');args=p.parse_args()
with sync_playwright() as pw:
 browser=pw.chromium.launch(executable_path=shutil.which('chromium'),headless=True,args=['--no-sandbox']);context=browser.new_context(viewport={'width':1440,'height':1000});page=context.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
 page.add_init_script("localStorage.setItem('mouse-rnaseq-v2',JSON.stringify({samples:'/legacy/samples.csv',counts:'/legacy/counts.csv',threads:8}));localStorage.setItem('mouse-rnaseq-current','aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa');localStorage.setItem('unrelated-setting','keep');")
 page.goto(args.base);page.wait_for_function('state.formEpoch>0')
 assert page.locator('#samples').input_value()=='' and page.locator('#threads').input_value()=='4'
 assert page.evaluate("localStorage.getItem('mouse-rnaseq-v2')") is None
 assert page.evaluate("localStorage.getItem('unrelated-setting')")=='keep'
 def dirty():
  page.locator('#samples').fill('/old/samples.csv');page.locator('#mode').select_option('counts');page.locator('#counts').fill('/old/counts.csv');page.locator('#threads').fill('8');page.locator('#covariates').fill('batch');page.locator('#addcontrast').click()
 def clean():
  assert page.locator('#samples').input_value()=='' and page.locator('#counts').input_value()==''
  assert page.locator('#covariates').input_value()=='' and page.locator('#threads').input_value()=='4'
  assert page.locator('#comparisonrows .comparison').count()==1
  assert page.locator('#samplepreview').is_hidden() and page.locator('#preflightresult').inner_text()==''
  assert page.evaluate('state.job') is None
 dirty();page.locator('#newproject').click();clean()
 dirty();page.reload();page.wait_for_function('state.formEpoch>0');clean()
 dirty();page.evaluate("window.dispatchEvent(new PageTransitionEvent('pageshow',{persisted:true}))");clean()
 held=[];page.route('**/api/metadata',lambda route:held.append(route))
 page.locator('#samples').fill('/pending/samples.csv');page.locator('#readmetadata').click();page.wait_for_timeout(150);assert held
 page.locator('#newproject').click();held.pop().fulfill(json=dict(samples=[dict(sample='OLD1',group='old')],groups=['old'],columns=['sample','group'],numeric=[]));page.wait_for_timeout(200);clean();assert page.evaluate('state.metadata') is None
 page.unroute('**/api/metadata')
 held=[];page.route('**/api/preflight',lambda route:held.append(route))
 page.locator('#preflight').click();page.wait_for_timeout(100);assert held;page.locator('#newproject').click();held.pop().fulfill(json=dict(samples=[{}]*6,comparisons=[{}],warnings=['STALE WARNING'],output_free_gb=100));page.wait_for_timeout(200);clean();assert 'STALE' not in page.locator('body').inner_text();page.unroute('**/api/preflight')
 # A job submission already sent can finish, while a new form remains blank.
 held=[]
 page.route('**/api/preflight',lambda route:route.fulfill(json=dict(samples=[{}]*6,comparisons=[{}],warnings=[],output_free_gb=100)))
 page.route('**/api/jobs',lambda route:held.append(route) if route.request.method=='POST' else route.continue_())
 page.locator('#samples').fill('/submit/samples.csv');page.locator('#mode').select_option('counts');page.locator('#counts').fill('/submit/counts.csv');page.locator('#submit').click();page.wait_for_timeout(150);assert held
 page.locator('#newproject').click();held.pop().fulfill(status=202,json=dict(id='a'*32));page.wait_for_timeout(250);clean();assert page.locator('#submit').is_enabled()
 page.unroute('**/api/preflight');page.unroute('**/api/jobs')
 # Explicit recovery is allowed, while refresh returns to a fresh form.
 page.locator('#project').select_option('a'*32);page.wait_for_function("document.getElementById('projectstatus').textContent.includes('失败')");page.locator('#loadsettings').click();page.wait_for_function("state.resumeFrom==='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'")
 assert page.locator('#samples').input_value() and page.locator('#reusecompleted').is_checked()
 page.locator('#newproject').click();clean();assert page.evaluate('state.resumeFrom') is None
 page.locator('#project').select_option('b'*32);page.wait_for_function("document.getElementById('projectstatus').textContent.includes('失败')");page.locator('#nextmodule').select_option('activity');page.locator('#useoutputs').click();page.wait_for_function("document.getElementById('analysis_stage').value==='activity'")
 assert page.locator('#expression').input_value().endswith('/vst.csv') and page.locator('#samples').input_value().endswith('/samples.csv')
 assert page.evaluate('state.upstreamFrom')=='b'*32
 page.locator('#newproject').click();clean()
 page.locator('#analysis_stage').select_option('enrichment');assert page.locator('#differential_input').is_visible() and page.locator('#samples').get_attribute('required') is None
 page.set_viewport_size({'width':390,'height':844});assert page.evaluate('document.documentElement.scrollWidth<=innerWidth')
 assert not errors,errors;browser.close()
print('PASS: initial/reload/new/BFCache reset, old-cache removal, stale response rejection, explicit failure recovery, upstream handoff, standalone form and mobile layout')
