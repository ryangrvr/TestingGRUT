"""Post-computation packaging checks; does not change the frozen estimator."""
import hashlib
import json
from pathlib import Path
import re
import sys
import audit

P = Path(__file__).resolve().parent


def main():
    audit.guard()
    manifest = json.loads((P / 'PACKET_MANIFEST.json').read_text())['files']
    actual_files = {p.name for p in P.iterdir() if p.is_file() and p.suffix in {'.md', '.py', '.json', '.csv'} and p.name != 'PACKET_MANIFEST.json'}
    if set(manifest) != actual_files:
        raise RuntimeError('Packet manifest membership differs')
    for name, expected in manifest.items():
        data = (P / name).read_bytes()
        if {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()} != expected:
            raise RuntimeError('Packet file changed: '+name)
    for path in P.glob('*.md'):
        for url in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            if '://' in url or url.startswith('#'):
                continue
            if not (path.parent / url.split('#')[0]).exists():
                raise RuntimeError('Broken relative link: '+path.name+' '+url)
    expected_prediction = json.loads((P / 'PREDICTION.json').read_text())
    read_original = audit.read
    test = (P / audit.inputs()[0]['test_file']).resolve()

    def reject_test_values(path):
        if path.resolve() == test:
            raise RuntimeError('Numerical prediction consumed test file')
        return read_original(path)

    audit.read = reject_test_values
    try:
        if audit.prediction() != expected_prediction:
            raise RuntimeError('Prediction changed with blocked test-file values')
    finally:
        audit.read = read_original
    result = {'frozen_source_and_public_identity': True,
              'packet_membership_and_hashes': True,
              'relative_document_links': True,
              'numerical_prediction_with_test_value_reads_blocked': True,
              'scope': 'Same-author posthoc packaging/dependency checks, not blind analysis or scientific replication.'}
    if '--write' in sys.argv:
        (P / 'QUALITY_VERIFICATION.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
