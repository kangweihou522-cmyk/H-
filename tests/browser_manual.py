"""Real-browser check for documentation with no analysis project selected."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from playwright.sync_api import sync_playwright

parser = argparse.ArgumentParser()
parser.add_argument('--base', default='http://127.0.0.1:18879')
parser.add_argument('--chromium', default=shutil.which('chromium'))
parser.add_argument('--output', type=Path, default=Path('runs/manual-browser-validation'))
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
errors = []
with sync_playwright() as p:
    browser = p.chromium.launch(executable_path=args.chromium, headless=True, args=['--no-sandbox'])
    page = browser.new_page(viewport={'width': 1440, 'height': 1000})
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.goto(args.base)
    assert page.request.get(args.base+'/api/jobs').json() == []
    assert page.evaluate('state.job') is None
    manifest = page.request.get(args.base+'/api/manual').json()['metadata']
    page.locator('nav [data-view="manual"]').click()
    page.wait_for_function('(version)=>document.getElementById("manualstatus").textContent.includes(version)', arg=manifest['guide_version'])
    frame = page.frame_locator('#manualframe')
    frame.locator('article h2').first.wait_for()
    assert frame.locator('article h2').count() >= 20
    frame.locator('#guide-search').fill('KEGG')
    assert frame.locator('#search-results a').count() > 0
    frame.locator('#search-results a').first.click()
    assert page.locator('#chartpanel').is_hidden() and page.locator('#tablepanel').is_hidden()
    # The periodic project refresh must not re-expose analysis-only controls.
    page.wait_for_timeout(4500)
    assert page.locator('#projectbanner').is_hidden() and page.locator('#steps').is_hidden()
    with page.expect_download() as downloaded:
        page.locator('#manualword').click()
    word = downloaded.value
    assert word.suggested_filename.endswith('.docx')
    assert hashlib.sha256(Path(word.path()).read_bytes()).hexdigest() == manifest['word_sha256']
    page.screenshot(path=str(args.output/'manual-desktop.png'), full_page=True)
    page.set_viewport_size({'width': 390, 'height': 844})
    assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
    assert frame.locator('body').evaluate('(body)=>body.scrollWidth <= innerWidth')
    page.screenshot(path=str(args.output/'manual-mobile.png'), full_page=True)
    # The standalone guide is independently readable without the workstation shell.
    page.goto(args.base+'/static/guide/index.html')
    page.locator('#guide-search').fill('火山图')
    assert page.locator('#search-results a').count() > 0
    assert not errors, errors
    browser.close()
print('PASS: no-project guide, chapter search, Word download, refresh, desktop/mobile, standalone page')
