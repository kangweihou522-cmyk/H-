"""Run with the installed system Chromium; no browser download required."""
from playwright.sync_api import sync_playwright
from pathlib import Path
import json
errors=[]
with sync_playwright() as p:
 browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
 page=browser.new_page(viewport={'width':1440,'height':1000});page.on('pageerror',lambda e:errors.append(str(e)))
 page.goto('http://127.0.0.1:8765');page.get_by_text('填入合成演示数据',exact=True).click();page.wait_for_function("document.getElementById('mode').value === 'counts'")
 assert page.locator('#samples').input_value().endswith('samples.csv')
 assert page.locator('#transcripts').is_hidden()
 page.locator('#mode').select_option('fastq');assert page.locator('#transcripts').is_visible()
 page.screenshot(path='/workspace/H-/runs/ui-desktop.png',full_page=True)
 page.set_viewport_size({'width':390,'height':844});page.screenshot(path='/workspace/H-/runs/ui-mobile.png',full_page=True)
 assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
 assert not errors,errors
 browser.close();print('Browser desktop/mobile, input mode switching, demo generation, no JS errors: PASS')
