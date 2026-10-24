#!/usr/bin/env python3
"""Add rigorous speed-of-light prefix bounds; compact the source tables."""
import json
import numpy as np
from scipy.integrate import quad
from streaming import ROOT,C,hubble

COLUMNS=['a','N','q_over_m','post_escape_to_observation_Mpc','prefix_path_mean_Mpc']


def load_run(path):
    r=json.loads(path.read_text())
    if r.get('escape_source_columns'):
        r['escape_source']=[dict(zip(r['escape_source_columns'],row)) for row in r['escape_source']]
        r.pop('escape_source_columns')
    return r


def light_prefix(af,ae):
    if ae<=af:return 0.
    return quad(lambda x:C/(np.exp(x)*hubble(np.exp(x))),np.log(af),np.log(ae),epsabs=1e-13,epsrel=2e-9)[0]


def main():
    rows=[]
    for path in sorted(ROOT.glob('RUN_*.json')):
        r=load_run(path);source=r['escape_source']
        weights=np.array([s['N'] for s in source]);ae=np.array([s['a'] for s in source])
        light=np.array([light_prefix(r['a_form'],a) for a in ae])
        for snap in r['snapshots']:
            keep=ae<=snap['a'];w=weights[keep]
            prefix=float(np.dot(w,light[keep])/max(sum(w),1e-100))
            snap['converted_mean_lightcone_prefix_bound_Mpc']=prefix
            snap['converted_total_traveled_path_mean_lower_Mpc']=snap['converted_post_escape_mean_Mpc']
            snap['converted_total_traveled_path_mean_upper_Mpc']=snap['converted_post_escape_mean_Mpc']+prefix
            snap['active_grid_path_estimates_certified']=False
        r['redshift_at_cumulative_escaped_fraction']=r.pop('last_inelastic_escape_redshift_quantiles',r.get('redshift_at_cumulative_escaped_fraction',{}))
        r['path_caution']='Active-grid path moments failed a free-redshift control and are diagnostics only. Physical conclusions use analytic post-escape horizons and a speed-of-light prefix bound. These path lengths do not fix a spatial transfer function.'
        # A compact numeric source table avoids repeated JSON field names.
        r['escape_source_columns']=COLUMNS
        r['escape_source']=[[s[c] for c in COLUMNS] for s in source]
        path.write_text(json.dumps(r,separators=(',',':'),allow_nan=False)+'\n')
        rows.append({k:v for k,v in r.items() if k not in ['escape_source','escape_source_columns','first_conversion_history']})
    old=json.loads((ROOT/'STREAMING_RESULTS.json').read_text())
    old['runs']=rows
    (ROOT/'STREAMING_RESULTS.json').write_text(json.dumps(old,indent=2,allow_nan=False)+'\n')
    control=json.loads((ROOT/'NO_REACTION_CONTROL.json').read_text())
    s=control['snapshots'][-1]
    error=s['all_X_total_path_mean_Mpc']/s['original_X_collisionless_Mpc']-1
    (ROOT/'REJECTED_PATH_DIAGNOSTIC.json').write_text(json.dumps(dict(
        control='Zero reactions must reproduce the independent collisionless horizon for both initial particles.',
        active_grid_path_relative_error=error,passed=False,
        consequence='Do not use active-grid path moments as a certified free-streaming result. Released conclusions use the separate comoving-momentum receiver, analytic quadrature and the causal prefix bound.'),indent=2)+'\n')
    print('Analyzed and compacted',len(rows),'production runs; rejected grid-path error',error)


if __name__=='__main__':main()
