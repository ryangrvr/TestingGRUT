#!/usr/bin/env python3
"""Selected witnesses in the original 1000-toy historical XENONnT release."""
import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from widths import ROOT, Widths, population
from scan import Object, PUB, detector

PARENT = ROOT.parent/'grut_portal_empirical_test_20261022'
spec = importlib.util.spec_from_file_location('model', PARENT/'code/model.py')
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
sys.modules['model'] = model
spec = importlib.util.spec_from_file_location('window_public_detector', PARENT/'code/detector.py')
original = importlib.util.module_from_spec(spec)
spec.loader.exec_module(original)


def main(cache):
    path = Path(cache)/'xe_likelihood/data/xent_nr_20ty_projection/XENONnT_binwise_NR_LL_toyMC.json.gz'
    provenance = json.loads((PARENT/'SOURCE_INPUTS.json').read_text())
    expected = provenance['sources']['xenon/'+str(path.relative_to(cache))]['sha256']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    with gzip.open(path, 'rt') as f:
        data = json.load(f)
    axis = data['Yedges']
    edges = np.round(np.linspace(axis['low_edge'], axis['high_edge'], axis['N']), 6)
    ys = (edges[1:]+edges[:-1])/2
    z = data['Z']['Data']
    toys = np.array([z[str(i)] for i in range(len(z))])
    del z, data
    scan = json.loads((ROOT/'SCAN_RESULTS.json').read_text())
    profiles = json.loads((ROOT/'PROFILE_RESULTS.json').read_text())['profiles']
    rows = []
    widths = Widths()
    for r in scan['selected']:
        p = next(p for p in profiles if p['Q'] == r['Q'] and p['lambda_phi'] == r['lambda_phi'])
        obj = Object(p)
        frac = r['witness_population']['object_halo_fraction']
        _, _, signals = detector(obj, r['v_phi_GeV'], frac)
        strength, curves = original.curve_grid(signals, ys, toys, 1200)
        ul = original.limits(strength, curves, sum(signals))
        kappas = 1e-9*np.sqrt(ul)
        pct = np.quantile(kappas, [.1586552539, .5, .8413447461])
        upper = min(r['kappa_LZ_same_mass_count_ESTIMATE'], r['kappa_LZ_shape_matched_ESTIMATE'] or np.inf)
        lower = max(pct[1], r['kappa_for_lifetime_1s_central'])
        joint_kappa = float(np.sqrt(lower*upper)) if lower < upper else None
        joint_pop = None if joint_kappa is None else population(r['Q'], obj.mass(r['v_phi_GeV']), r['v_phi_GeV'], r['lambda_phi'], joint_kappa, widths.unit(r['m_s_GeV']))
        row = dict(Q=r['Q'], lambda_phi=r['lambda_phi'], v_phi_GeV=r['v_phi_GeV'],
                   kappa_projected90_percentiles=pct.tolist(),
                   conservative_of_two_LZ_screen_ESTIMATES=upper,
                   prompt_and_median_projection_screen_overlap=bool(lower < upper),
                   overlap_lower_kappa=float(lower), overlap_upper_kappa=float(upper),
                   witness_kappa=r['witness_kappa'], witness_above_median_projected_limit=bool(r['witness_kappa'] >= pct[1]),
                   expected_signals_at_median_projection=float(sum(signals)*(pct[1]/1e-9)**2))
        row['joint_screen_kappa'] = joint_kappa
        row['joint_screen_population'] = joint_pop
        row['joint_screen_expected_signals'] = None if joint_kappa is None else float(sum(signals)*(joint_kappa/1e-9)**2)
        rows.append(row)
        print(row, flush=True)
    out = dict(scope='Original 1000 public background-only projected likelihood toys, historical 20 tonne-year study; not a discovery forecast or contemporary observation. LZ mapping is approximate.',
               toys=len(toys), source_sha256=expected, rows=rows)
    (ROOT/'PROJECTION_RESULTS.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cache', type=Path)
    main(parser.parse_args().cache)
