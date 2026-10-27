"""Retain measured output and null physical conclusions."""
from cw_common import *

def main():
    final=read('FINAL_AUDIT.json')['cases'];history=read('HISTORY_AUDIT.json')['cases'];rows=[]
    for tag,r in final.items():
        picked=[c for c in r['candidates'] if 'spheres' in c]
        rows.append(dict(tag=tag,N=r['N'],L=r['L'],t=r['t'],candidate_count=len(r['candidates']),
            scalar_Hamiltonian_fraction=r['whole_box']['E_phi']/r['whole_box']['E_total'],
            largest_compact_masks=[dict(id=c['candidate']['id'],Q=c['mask']['Q'],
                E_over_Q=c['mask']['E_over_Q'],charged_E_over_Q=c['mask']['E_X_over_Q'],
                omega_fit=c['mask']['omega_fit'],kinetic_excess_fraction=c['mask']['kinetic_excess_fraction'],
                stationary_charge_force_residual=c['mask']['stationary_X_force_residual_over_m2'],
                positive_current_phase_coherence=c['mask']['phase_coherence_weighted_positive_current'],
                neutral_sign_shell=c['neutral_sign_shell']) for c in picked]))
    main=next(r for r in rows if r['tag']=='mixed_N96')
    save('RESULTS.json',dict(accepted_review='870a2ab51b74c3ab0a3d33fa8997ece376b80fb1',
        source_pilot='815955656099aeb684d909bf21c940ca9add8e87',
        review_static_elimination_wording_corrected=True,new_dynamical_histories=0,
        rows=rows,main_compact_masks=main['largest_compact_masks'],
        exact_local_flux_and_exchange_implemented=True,
        actual_core_wave_decomposition_derived=False,physical_bound_population=None,
        physical_particle_spectrum=None,physical_f_sym=None,physical_kinetic_threshold=None,
        scope='Existing cold fields, instantaneous necessary diagnostics and an exact ledger. No core mass, asymptotic particle count or cosmological abundance inferred.'))
    table=['| Cold control | Scalar energy / total | Excess charged kinetic fraction, selected masks |',
           '|---|---:|---|']
    for r in rows:
        table.append('| '+r['tag']+' | '+f"{100*r['scalar_Hamiltonian_fraction']:.2f}%"+' | '+', '.join(f"{100*c['kinetic_excess_fraction']:.1f}%" for c in r['largest_compact_masks'])+' |')
    (ROOT/'RESULT_TABLE.md').write_text('\n'.join(table)+'\n')
    print('SUMMARY',[(r['tag'],r['scalar_Hamiltonian_fraction'],[c['kinetic_excess_fraction'] for c in r['largest_compact_masks']]) for r in rows])

if __name__=='__main__':main()
