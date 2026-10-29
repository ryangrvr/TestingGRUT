"""Audit clocks, threshold sensitivities, and parent/new modeled energy ledgers."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_m1_expansion_clock_20261028'
OLD=ROOT.parent/'grut_m1_local_cloud_20261029'
NC=1.2141071422639381e-15
EREF=187479.9981052459
def read(p):return json.loads(p.read_text())

def main():
    parent=read(PARENT/'profiled_h0025_96_to10.json');new=read(ROOT/'late_kin_h0025_96.json')
    oldsum=read(PARENT/'SUMMARY.json');p=parent['snapshots'][-1];n=new['snapshots'][-1]
    oldrow=next(x for x in oldsum['runs'] if x.get('tag')=='profiled_h0025_96_to10') if 'runs' in oldsum else None
    if oldrow is None:
        # Find the retained entry without relying on a particular top-level key.
        oldrow=next(x for rows in oldsum.values() if isinstance(rows,list) for x in rows
                    if isinstance(x,dict) and x.get('tag')=='profiled_h0025_96_to10')
    rnew=new['held_unassigned_recomputed_mX'];rold=oldrow['unassigned_initial_energy_mX']
    decompositions={}
    for label,s,r,full in [('parent',p,rold,oldrow['total_excess_including_held_unassigned_percent']),
                           ('through_20000',n,rnew,new['full_modeled_final_excess_percent_with_remainder'])]:
        decompositions[label]={'Qcore':s['Qcore'],'Nplus':s['Nplus'],'Nminus_solver_residual':s['Nminus'],
             'cold_excess_percent':100*s['stationary_cold_excess_fraction'],
             'gas_kinetic_percent':100*s['gas_kinetic_mX']/EREF,'held_remainder_percent':100*r/EREF,
             'full_percent':full,'held_remainder_mX':r,'work_redshift_mX':s['work_redshift_mX']}
    threshold=1+2*.8817204919524214;late=[]
    for width in (200,400,800):
        z=np.load(ROOT/f'late_h0025_t{width}.npz');e=z['energy_centers'];counts=z['plus_positive'].sum(axis=0)
        d=read(ROOT/f'late_h0025_t{width}.json')
        late.append({'window':width,'positive_projection_above_reference_reverse_threshold_approx':float(counts[e>threshold].sum()),
                     'plus':d['species']['plus']['positive_projection_number'],
                     'minus':d['species']['minus']['positive_projection_number'],
                     'neutral':d['species']['neutral']['positive_projection_number']})
    cr=read(OLD/'cross_e280_q64.json')['result']['rows'][1];rate=cr['instantaneous_final_rate_per_mX_t']
    early=cr['during_source_events'];T1=.01/rate+10000-early/rate
    previous={str(round(x['tau'],3)):x['Nminus'] for x in parent['snapshots'] if x['tau']<=1000.00001}
    out={'ncore_over_mX3':NC,'unit_optical_time_mX_t':1/NC,
         'tau_at_mX_t_1e7':NC*1e7,'tau_at_mX_t_1e10':NC*1e10,
         'hypothetical_time_to_0p01_distinct_core_events':T1,
         'hypothetical_events_at_mX_t_1e10':early+(1e10-10000)*rate,
         'clock_scope':'tau is the model density clock, not a universal reaction optical depth. The 1e10 cross-core extrapolation freezes all direct spectra and neglects interactions, new radiation, daughters and expansion.',
         'parent_Nminus_at_sampled_tau':previous,
         'reference_exchange_threshold_mX':threshold,'late_source_time_window_controls':late,
         'mass_reference_mX':EREF,'energy_decompositions_at_scale_factor_10':decompositions,
         'scope':'Positive projected late source, conditional parent expansion transport; no formation, complete angular matching, core backreaction or all-future source.'}
    (ROOT/'SCENARIO_AUDIT.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
