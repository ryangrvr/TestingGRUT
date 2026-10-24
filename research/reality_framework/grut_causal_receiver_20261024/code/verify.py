#!/usr/bin/env python3
import argparse,hashlib,json
import numpy as np
from run import ROOT


def read(name):return json.loads((ROOT/name).read_text())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--with-raw',action='store_true');args=parser.parse_args()
    layer=read('DISCRETE_LAYER_RESULTS.json');operator=read('OPERATOR_CONTROLS.json');audit=read('SOURCE_AUDIT_RESULTS.json');interface=read('I025_RESULTS.json')
    checks=[]
    def check(name,condition,evidence):checks.append(dict(name=name,passed=bool(condition),evidence=evidence))
    maxR=max(r['reflected_energy'] for r in layer['rows']);zero=max(r['zero_layer_complex_error'] for r in layer['zero_layer_controls'])
    check('Actual discrete massive-wave layer meets sampled reflection target with exact lattice leads',maxR<.01 and zero<1e-8,dict(max_reflected_energy=maxR,zero_layer_amplitude_error=zero))
    local=operator['local']
    check('Actual executable has the Hamiltonian gradient and energy/charge boundary-work identities',local['gradient_energy_identity_relative_error']<1e-11 and max(v for k,v in local.items() if k!='gradient_energy_identity_relative_error')<1e-7,local)
    spectra=operator['spectral_controls'];growth=max(r['maximum_real_eigenvalue'] for r in spectra);amplification=max(r['maximum_RK4_eigen_amplification'] for r in spectra)
    check('Representative complete finite-grid linear spectra have no growing RK4 modes',growth<=1e-9 and amplification<=1+1e-10,dict(max_real_eigenvalue=growth,max_step_amplification=amplification))
    packets=[read(f'{name}_RESULTS.json') for name in ['free_X','free_X_time','free_s','free_soft']]
    packet_error=max(r['free_comparison']['exact_time_physical_energy_error_over_initial'] for r in packets)
    boundary_error=max(r['free_comparison']['same_RK4_bulk_physical_energy_error_over_initial']+r['free_comparison']['reference_fast_spectral_tail_energy_fraction'] for r in packets)
    check('Time-domain packets match independent free propagation, including the 60000-unit soft case',packet_error<1e-3 and boundary_error<1e-6,dict(max_full_physical_energy_error_fraction=packet_error,max_same_RK4_error_plus_reference_tail=boundary_error))
    refined=[read(f'{name}_RESULTS.json') for name in ['P025_time','P0125','P00625','P025_remote_time']]
    ee=max(r['maximum_energy_plus_receiver_relative_error'] for r in refined);qe=max(r['maximum_charge_plus_receiver_relative_error'] for r in refined)
    check('Refined formation runs close physical-field plus signed boundary ledgers',ee<1e-5 and qe<1e-6,dict(max_energy_error=ee,max_charge_error=qe))
    grid=[read(f'{name}_RESULTS.json') for name in ['P025_time','P0125','P00625']]
    values=[r['core60_late_mean_fraction'] for r in grid];spread=max(values)/min(values)-1
    check('Three spatial grids agree on the declared finite-time retained charge to one percent',spread<.01,dict(h=[r['parameters']['h'] for r in grid],retained_fractions=values,relative_spread=spread,Richardson_certified=False))
    radius=audit['boundary_radius_control']
    check('Moving the layer from radius 360 to 720 preserves the controlled finite-time interior',abs(radius['core_charge_relative_shift'])<.001 and radius['raw_complex_field_L2_squared_relative_difference_inside_300']<1e-5,radius)
    cut=interface['maximum_cut_budget_relative_error'];identities=interface['field_dynamics_bitwise_identical_to_P025_time']
    check('Exact extraction cut closes energy/charge without altering field dynamics',max(cut.values())<1e-5 and all(identities.values()),dict(cut_errors=cut,bitwise_identities=identities))
    parseval=max(v[k] for row in audit['rows'] for v in row['variants'] for k in ['relative_energy_parseval_error','relative_charge_parseval_error'])
    modes=read('BOUND_MODE_DIAGNOSTICS.json');mode_error=max(r['eigenvector_residual'] for r in modes['rows'])
    check('Classical-current Fourier accounting and independently applied bound-mode stencil are consistent',parseval<1e-12 and mode_error<1e-8,dict(Parseval_error=parseval,bound_mode_stencil_residual=mode_error))
    raw=dict(performed=False,arrays=0)
    if args.with_raw:
        count=0
        for path in sorted(ROOT.glob('*_RESULTS.json')):
            r=json.loads(path.read_text())
            if 'raw_numeric_sha256' not in r:continue
            for key,expected in r['raw_numeric_sha256'].items():
                array=np.load(ROOT/'data'/f"{r['name']}_{key}.npz")[key]
                assert hashlib.sha256(array.tobytes()).hexdigest()==expected,(r['name'],key)
                count+=1
        raw=dict(performed=True,arrays=count)
    coarse=[]
    for name in ['P025','P025_remote']:
        r=read(f'{name}_RESULTS.json');coarse.append(dict(name=name,energy_error=r['maximum_energy_plus_receiver_relative_error'],passes_refined_energy_target=r['maximum_energy_plus_receiver_relative_error']<1e-5))
    # Keep the coarse failures and unresolved physical source outside the pass count.
    out=dict(scope='Author numerical controls of a causal exterior and finite-time prepared core; not independent replication, an asymptotic particle source or a completed cosmology.',
        checks=checks,passed=sum(c['passed'] for c in checks),total=len(checks),raw_validation=raw,
        coarse_budget_controls=coarse,long_soft_packet_ledger_energy_drift=packets[-1]['maximum_energy_plus_receiver_relative_error'],
        complete_asymptotic_particle_source_certified=False,nonlinear_vacuum_matching_theorem=False,Richardson_certified=False)
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for c in checks:print('PASS' if c['passed'] else 'FAIL',c['name'],c['evidence'],flush=True)
    print('RAW',raw,'COARSE_BUDGET_CONTROLS',coarse,flush=True)
    assert all(c['passed'] for c in checks)


if __name__=='__main__':main()
