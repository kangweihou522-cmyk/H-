"""Manual routes, portable Word contents and single-source synchronization."""
import io
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import pytest
import app as station

ROOT = Path(__file__).resolve().parent.parent
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}


def test_manual_available_without_analysis_and_downloads_match():
    client = station.app.test_client()
    manifest = json.loads((ROOT/'docs/guide-manifest.json').read_text())
    status = client.get('/api/manual')
    assert status.status_code == 200
    assert status.json['word'] and status.json['source'] and status.json['html']
    assert status.json['metadata']['manual_sha256'] == manifest['manual_sha256']
    word = client.get('/api/manual/word')
    assert word.status_code == 200 and word.data == (ROOT/'docs/USER_GUIDE.docx').read_bytes()
    assert word.mimetype == 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
    assert 'attachment;' in word.headers['Content-Disposition']
    source = client.get('/api/manual/source')
    assert source.status_code == 200 and source.data == (ROOT/'docs/USER_GUIDE.md').read_bytes()
    assert client.get('/static/guide/index.html').status_code == 200


def test_missing_manual_has_clear_response(tmp_path, monkeypatch):
    monkeypatch.setattr(station, 'ROOT', tmp_path)
    client = station.app.test_client()
    assert client.get('/api/manual').json['word'] is False
    assert client.get('/api/manual/word').status_code == 404
    assert '更新' in client.get('/api/manual/source').json['error']


def test_word_contains_all_chapters_tables_links_and_utf8_content():
    source = (ROOT/'docs/USER_GUIDE.md').read_text()
    chapters = re.findall(r'^## (.+)$', source, flags=re.M)
    with ZipFile(ROOT/'docs/USER_GUIDE.docx') as package:
        assert package.testzip() is None
        xml = ET.fromstring(package.read('word/document.xml'))
        text = ''.join(xml.itertext())
        for title in chapters:
            assert title in text
        assert len(xml.findall('.//w:tbl', NS)) >= 5
        bookmarks = {b.attrib['{'+NS['w']+'}name'] for b in xml.findall('.//w:bookmarkStart', NS)}
        anchors = {a.attrib['{'+NS['w']+'}anchor'] for a in xml.findall('.//w:hyperlink', NS) if '{'+NS['w']+'}anchor' in a.attrib}
        assert anchors and anchors <= bookmarks
        relationships = package.read('word/_rels/document.xml.rels').decode()
        assert 'https://doi.org/' in relationships and 'https://github.com/' in relationships
        assets = re.findall(r'!\[[^\]]*\]\((assets/[^)]+)\)', source)
        embedded = {package.read(name) for name in package.namelist() if name.startswith('word/media/')}
        assert len(assets) == 4
        for asset in assets:
            original = (ROOT/'docs'/asset).read_bytes()
            assert original in embedded
            assert original == (ROOT/'static/guide'/asset).read_bytes()


def test_teaching_example_respects_bh_order_and_both_thresholds():
    source = (ROOT/'docs/USER_GUIDE.md').read_text()
    rows = re.findall(r'^\| ([A-Z]) \| (-?[\d.]+) \| ([\d.]+) \| ([\d.]+|NA) \|', source, flags=re.M)
    assert len(rows) == 7
    rows = [(gene, float(lfc), float(p), None if adj == 'NA' else float(adj)) for gene, lfc, p, adj in rows]
    finite = sorted((p, adj) for _, _, p, adj in rows if adj is not None)
    assert all(0 <= p <= adj <= 1 for p, adj in finite)
    assert all(a[1] <= b[1] for a, b in zip(finite, finite[1:]))
    nominal = {g for g, lfc, p, _ in rows if p < .05 and abs(lfc) > 0}
    strict = {g for g, lfc, _, adj in rows if adj is not None and adj < .05 and abs(lfc) > 1}
    assert nominal == {'A', 'B', 'C', 'E', 'H'} and strict == {'B', 'C'}


def test_generation_is_repeatable_and_detects_changed_runtime(tmp_path):
    pytest.importorskip('docx')
    from scripts.build_user_guide import build, OUTPUTS
    docs = tmp_path/'docs'; docs.mkdir()
    (docs/'USER_GUIDE.md').write_text('# 使用说明\n\n说明书版本：2026.10.26-r2\n\n## 输入与设计\n\n中文与 `FASTQ`。\n\n| A | B |\n|---|---|\n| 原始计数 | 2 |\n')
    runtime = tmp_path/'app.py'; runtime.write_text('version = 1\n')
    first = build(tmp_path)
    assert first == build(tmp_path)
    with ZipFile(io.BytesIO(first[OUTPUTS[1]])) as word:
        assert word.testzip() is None and '中文' in word.read('word/document.xml').decode()
    runtime.write_text('version = 2\n')
    changed = build(tmp_path)
    assert first[OUTPUTS[0]] != changed[OUTPUTS[0]] and first[OUTPUTS[1]] != changed[OUTPUTS[1]]
    old, new = (json.loads(item[OUTPUTS[2]]) for item in (first, changed))
    assert old['guide_version'] == '2026.10.26-r2' and old['updated'] == '2026-10-26'
    assert old['manual_sha256'] == new['manual_sha256']
    assert old['code_fingerprint'] != new['code_fingerprint']
