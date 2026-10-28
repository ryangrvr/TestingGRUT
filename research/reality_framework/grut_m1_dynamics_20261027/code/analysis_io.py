"""Lossless compressed catalogs; keep the full numerical measurements."""
from pathlib import Path
import gzip,json

ROOT=Path(__file__).resolve().parents[1]

def read_analysis():
    compressed=ROOT/'ANALYSIS.json.gz'
    if compressed.exists():
        with gzip.open(compressed,'rt',encoding='utf-8') as f:return json.load(f)
    return json.loads((ROOT/'ANALYSIS.json').read_text())

def write_analysis(data):
    raw=json.dumps(data,separators=(',',':'),allow_nan=False).encode()
    (ROOT/'ANALYSIS.json.gz').write_bytes(gzip.compress(raw,mtime=0))

