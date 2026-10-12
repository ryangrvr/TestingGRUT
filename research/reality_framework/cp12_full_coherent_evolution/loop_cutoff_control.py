"""Reproduce exact L=3 versus L=2 perturbation-coefficient comparison."""
from loop_series import series
from pathlib import Path
import json
P=Path(__file__).resolve().parent
a=series(3,L=3,order=14);b=series(6,L=3,order=14)
old=json.loads((P/'loop_series_results.json').read_text())
assert all(str(b[k]-2*a[k])==old[k]['ring6_minus_two_triangles'] for k in range(15))
(P/'loop_series_cutoff_control.json').write_text(json.dumps(dict(L=3,order=14,exact_match_L2=True,differences=[str(b[k]-2*a[k]) for k in range(15)]),indent=2)+'\n')
print('EXACT_CUTOFF_MATCH_PASS')
