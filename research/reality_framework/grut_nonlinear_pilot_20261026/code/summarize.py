"""Summarize measured candidate statistics and retain unsolved physical outputs."""
from pilot import *

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];analysis=json.loads((ROOT/'ANALYSIS.json').read_text())['cases'];rows=[]
    for tag,case in analysis.items():
        rr=case['rows'];last=rr[-1]
        values=np.interp([160.,180.,200.,220.,240.],[r['t'] for r in rr],[r['catalogs'][1]['Q_fraction_in_masks'] for r in rr]);catalog=last['catalogs'][1]
        rows.append(dict(tag=tag,N=last['N'],L=last['L'],light_crossings=last['t']/last['L'],
            charge_drift=runs[tag]['charge_relative_drift'],energy_drift=runs[tag]['energy_relative_drift'],
            first_candidate_snapshot_time=next((r['t'] for r in rr if r['catalogs'][1]['candidate_count']),None),
            final_count=catalog['candidate_count'],final_Qmask_over_netQ=catalog['Q_fraction_in_masks'],
            late_Qmask_over_netQ_mean=float(values.mean()),late_Qmask_over_netQ_std=float(values.std()),
            late_sampling='Interpolated stored snapshots at times160,180,200,220,240; scatter is not an independent ensemble error.',
            final_candidate_Qmask=[c['Q_mask'] for c in catalog['candidates']],
            final_candidate_Emask_over_Qmask=[c['E_mask_over_Qm'] for c in catalog['candidates']],
            final_threshold_Qmask_over_netQ={str(c['threshold']):c['Q_fraction_in_masks'] for c in last['catalogs']},
            final_phi_negative_fraction=last['phi_negative_fraction'],
            final_boundary_shell_075=last['spectra']['boundary_shell_energy_fractions']['0.75'],
            final_interaction_residual=last['spectra']['vacuum_reference']['absolute_interaction_residual_over_ref']))
    rows.sort(key=lambda r:r['tag'])
    checks=json.loads((ROOT/'ENGINE_CHECKS.json').read_text())['groups'];checks.update(json.loads((ROOT/'ANALYSIS_CHECKS.json').read_text())['groups'])
    q=max(r['charge_relative_drift'] for r in runs.values());e=max(r['energy_relative_drift'] for r in runs.values())
    checks['nonlinear_conservation']=dict(passed=q<1e-9 and e<1e-3,charge=q,energy=e,
        scope='Discrete cold Hamiltonian conservation, not convergence of a quantum plasma or object abundance.')
    recon=[];parseval=[]
    for tag,case in analysis.items():
        last=case['rows'][-1];source=runs[tag]['final'];recon.append(abs(last['rho']/source['rho']-1))
        parseval.extend(r['spectra']['charge_Parseval_error'] for r in case['rows'])
    checks['snapshot_and_spectral_consistency']=dict(passed=max(recon)<1e-10 and max(parseval)<1e-8,
        energy_reconstruction=max(recon),charge_spectrum_Parseval=max(parseval))
    save('VERIFICATION_RESULTS.json',dict(groups=checks,independent_reviewer_ran=False,
        physical_population_converged=False,physical_f_sym=None))
    save('RESULTS.json',dict(rows=rows,initial_energy_m4=C.EX+C.V0,initial_charge_m3=C.Q0,
        source_parent='23356538581c64e49c255e6cc895373d5d297c96',
        nonlinear_3D_fields_evolved=True,object_charge_imposed=False,cold_seed_histories_only=True,
        candidate_definition='Gaussian-smoothed current>4mean and smoothed phi²<.25v²; masks connected on periodic lattice. Other thresholds2,8 retained.',
        bound_object_population_certified=False,physical_initial_particle_spectrum=None,
        physical_f_sym=None,physical_kinetic_threshold=None,
        scope='Finite cold classical boxes after several light crossings. Candidate-mask charge is not final retained charge or the symmetric-matter fraction. No continuum abundance or isolated core mass has been obtained.'))
    print('SUMMARY',[(r['tag'],r['final_count'],r['late_Qmask_over_netQ_mean'],r['final_interaction_residual']) for r in rows],flush=True)
    if not all(g['passed'] for g in checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
