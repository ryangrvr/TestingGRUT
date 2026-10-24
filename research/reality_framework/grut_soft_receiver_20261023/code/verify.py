#!/usr/bin/env python3
import argparse,hashlib,json
import numpy as np
from boundary import ROOT


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--with-traces',action='store_true');args=parser.parse_args()
    boundary=json.loads((ROOT/'BOUNDARY_RESULTS.json').read_text())
    audit=json.loads((ROOT/'TRACE_AUDIT_RESULTS.json').read_text())
    controls=boundary['controls'];checks=[]
    def check(name,ok,evidence):checks.append(dict(name=name,passed=bool(ok),evidence=evidence))
    exact=max(c['complex_amplitude_error'] for c in controls if c['kind']!='independent_grid_and_tolerance')
    check('Independent zero-damping and constant-slab analytic solutions',exact<1e-8,exact)
    grid=[c for c in controls if c['kind']=='independent_grid_and_tolerance']
    finest=max(c['complex_amplitude_errors'][-1] for c in grid)
    order=min(min(c['observed_orders']) for c in grid)
    check('Independent boundary-value grids converge to the layer ODE',finest<1e-4 and order>1.8,dict(max_finest_complex_amplitude_error=finest,minimum_observed_order=order))
    flux=max(c.get('flux_balance_absolute_error',0) for c in controls)
    tolerance=max(c.get('reflected_energy_tolerance_difference',0) for c in controls)
    passive=all(0<=s['reflected_energy']<=1+1e-10 for row in boundary['cases'] for s in row['samples'])
    check('Positive flux loss, passivity and ODE tolerance control',flux<1e-8 and tolerance<1e-8 and passive,dict(max_flux_error=flux,max_tolerance_shift=tolerance))
    parseval=max(v[key] for row in audit['rows'] for v in row['variants'] for key in ['relative_energy_parseval_error','relative_charge_parseval_error'])
    check('Measured signed currents agree with time/Fourier Parseval accounting',parseval<2e-12,parseval)
    raw=dict(performed=False,trace_count=0,two_probe_fields_checked=False)
    if args.with_traces:
        parent=ROOT.parent/'grut_core_conversion_20261023'
        for row in audit['rows']:
            trace=np.load(parent/'data'/f"{row['name']}_trace.npz")['trace']
            assert hashlib.sha256(trace.tobytes()).hexdigest()==row['trace_sha256']
        matched=audit['matched_two_probe_control'];columns=[0,1,2,3,4,7,8,9,10,11,12,13,14]
        for name in matched['names']:
            final=np.loadtxt(parent/'data'/f'{name}_final.csv',delimiter=',',skiprows=1)
            history=np.loadtxt(parent/'data'/f'{name}_history.csv',delimiter=',',skiprows=1)
            assert hashlib.sha256(final.tobytes()).hexdigest()==matched['final_numeric_sha256']
            assert hashlib.sha256(history[:,columns].tobytes()).hexdigest()==matched['physical_history_numeric_sha256']
        raw=dict(performed=True,trace_count=len(audit['rows']),two_probe_fields_checked=True)
    receiver_gate=all(row['passes_one_percent_gate_on_sampled_band'] for row in boundary['cases'] if row['name']=='original')
    out=dict(scope='Author numerical controls. The receiver/source gate is separate and fails.',
        numerical_controls=checks,passed=sum(c['passed'] for c in checks),total=len(checks),raw_validation=raw,
        original_receiver_one_percent_gate_passed=receiver_gate,complete_asymptotic_particle_source_certified=False,
        warnings=['Continuum linear boundary benchmark is not a full discrete nonlinear formation calculation.',
            'Conditional pair estimates are not a certified antiparticle abundance.',
            'The matched-layer formula is an analytic design reference; no time-domain PML is implemented.'])
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for c in checks:print('PASS' if c['passed'] else 'FAIL',c['name'],c['evidence'])
    print('RAW_VALIDATION',raw,'ORIGINAL_RECEIVER_GATE',receiver_gate)
    assert all(c['passed'] for c in checks)


if __name__=='__main__':main()
