"""Optional real-browser smoke check; supply installed Chromium or Playwright's browser."""
import argparse
import shutil
from pathlib import Path
from playwright.sync_api import sync_playwright

parser = argparse.ArgumentParser()
parser.add_argument('--base', default='http://127.0.0.1:8765')
parser.add_argument('--chromium', default=shutil.which('chromium') or shutil.which('chromium-browser'))
parser.add_argument('--output', type=Path, default=Path('runs/browser-validation'))
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
errors = []
with sync_playwright() as p:
    browser = p.chromium.launch(executable_path=args.chromium, headless=True, args=['--no-sandbox'])
    page = browser.new_page(viewport={'width': 1440, 'height': 1000})
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto(args.base)
    page.locator('#demo').click()
    page.wait_for_function("document.getElementById('mode').value === 'counts'")
    page.wait_for_function("!document.getElementById('samplepreview').hidden")
    assert page.locator('#samples').input_value().endswith('samples.csv')
    assert page.locator('#transcripts').is_hidden()
    page.locator('#mode').select_option('fastq')
    assert page.locator('#transcripts').is_visible()
    page.screenshot(path=str(args.output/'desktop.png'), full_page=True)
    page.set_viewport_size({'width': 390, 'height': 844})
    page.screenshot(path=str(args.output/'mobile.png'), full_page=True)
    assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
    assert not errors, errors
    browser.close()
print('PASS: demo metadata, input modes, desktop/mobile layout, no JavaScript errors')
