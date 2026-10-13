from playwright.sync_api import sync_playwright
from pathlib import Path
import json,urllib.request
j=Path('/tmp/rnaseq-smoke-job').read_text()
with sync_playwright() as p:
 b=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
 page=b.new_page(viewport={'width':1440,'height':1000});errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
 page.goto('http://127.0.0.1:8765');page.locator('.job').filter(has_text=j[:8]).locator('button').click();page.wait_for_selector('#gallery .figure')
 page.wait_for_function("document.querySelectorAll('#gallery .figure').length===7")
 card=page.locator('#gallery .figure').nth(1);card.get_by_label('标题',exact=True).fill('Edited volcano from browser');card.get_by_label('上调/组1颜色').fill('#dd3355');card.get_by_text('重绘并保存',exact=True).click()
 card.locator('.figure-status').filter(has_text='已保存').wait_for(timeout=30000)
 text=urllib.request.urlopen('http://127.0.0.1:8765/api/jobs/'+j+'/files/volcano.svg').read().decode();assert 'Edited volcano from browser' in text and '#dd3355' in text
 with page.expect_download() as d:card.get_by_text('下载 PDF',exact=True).click()
 d.value.save_as('/tmp/browser-volcano.pdf');assert Path('/tmp/browser-volcano.pdf').read_bytes().startswith(b'%PDF')
 card.screenshot(path='/workspace/H-/runs/ui-editor.png');assert not errors,errors;b.close()
 print('7 result cards, browser title/color redraw, SVG contents and PDF download: PASS')
