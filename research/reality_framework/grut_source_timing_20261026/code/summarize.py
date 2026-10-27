"""Keep conditional chronology separate from an actual relic history."""
from timing_common import *

def main():
    response=read('RESPONSE_RESULTS.json');clocks=read('SOURCE_CLOCKS.json');scan=read('QUARTIC_SCAN.json');rows=[]
    for row in response['rows']:
        if not row['tag'].startswith('quiet'):continue
        clock=next(c for c in clocks['histories'] if c['tag']==row['tag'])
        rows.append(dict(tag=row['tag'],beta=row['beta'],Gamma_over_H=row['Gamma_over_H'],late_C=row['late_C'],
            supplied_phase_crossings=row['crossings'],homogeneous_NX_final_change=row['homogeneous_NX_final_change_from_norm'],
            decay_number_medians={k:v['number_quantiles'][2] for k,v in clock['sources'].items()},
            source_fractions_by_003rad_unity={k:next(c for c in v['at_response_crossings'] if c['phase_seed']==.03 and c['kind']=='amplitude_unity') for k,v in clock['sources'].items()},
            wavelength_diagnostics=clock['wavelength_diagnostics']))
    save('RESULTS.json',dict(parent='423d37b78950d30552179d4de0c1c465016bca6d',
        accepted_reviews=['a5cd6a1e2f406fa9e6a4c88da092fed5d3a519fc','bca4927214695ef7c668ebb1d20f26a9cba91ca7'],
        previously_pending_reviews_completed=True,charged_mode_histories=42,background_controls=4,
        end_tau=800,rows=rows,auxiliary_quartic_sampled_bands=scan['sampled_bands'],
        physical_f_sym=None,physical_kinetic_threshold=None,physical_nonlinear_breakup_time=None,
        physical_initial_covariance=None,physical_bound_population=None,
        production_source_normalizations_recomputed=False,
        scope='Conditional collisionless finite-time responses and pinned-width decay clocks. Demonstrates loss of the homogeneous linear description before prescribed heating for tested seeds, not nonlinear fragmentation products or final abundance.'))
    lines=['| Quiet history | .03-rad amplitude reaches1 by sampled tau | Scale factor at that sample | 1e-6-rad amplitude reaches1 by sampled tau |',
           '|---|---:|---:|---:|']
    for r in rows:
        c=r['supplied_phase_crossings'];a=c[0]['amplitude_unity'];b=c[1]['amplitude_unity']
        lines.append(f"| ({r['beta']:g},{r['Gamma_over_H']:g}) | {a['tau_bracket'][1]:g} | {a['a_bracket'][1]:.3f} | {b['tau_bracket'][1]:g} |")
    (ROOT/'RESULT_TABLE.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
