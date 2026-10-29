"""Summarize pressure, wave-return and modal controls from generated outputs."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
W=ROOT/'.work'

def table(tag,kind):return np.genfromtxt(W/(tag+'_'+kind+'.csv'),delimiter=',',names=True)
def window(data,r,a,b):
    x=data[(data['tau']>=a-1e-9)&(data['tau']<=b+1e-9)]; t=x['tau']; assert abs(t[0]-a)<1e-8 and abs(t[-1]-b)<1e-8
    avg=lambda k:float(np.trapezoid(x[k],t)/(b-a))
    tr=avg(f'trace{r}_center');surface=avg(f'surface{r}');rate=float((x[f'virial{r}'][-1]-x[f'virial{r}'][0])/(b-a))
    return {'trace_mean_mX':tr,'pressure_integral_mean_mX':tr/3,'trace_face_mean_mX':avg(f'trace{r}_face'),
            'surface_mean_mX':surface,'virial_change_rate_mX':rate,'virial_identity_residual_mX':tr-surface-rate,
            'instantaneous_trace_min_mX':float(x[f'trace{r}_center'].min()),'instantaneous_trace_max_mX':float(x[f'trace{r}_center'].max()),
            'fraction_samples_negative':float(np.mean(x[f'trace{r}_center']<0)),'samples':len(x)}

def main():
    out={'scope':'Flat-space finite-region diagnostics. Energy and integrated trace are scaled by 1/g^2=100; pressure is trace/3. The surface term and virial change are retained. Do not identify an unweighted finite-time mean with integrated FRW work.', 'early':{},'late':{},'stationary':{},'replay':{},'feedback':{},'axis':{}}
    arrays={}
    for tag,h in [('sampled01',.1),('sampled005',.05),('sampled0025',.025)]:
        s=table(tag,'stress');out['early'][str(h)]={str(r):window(s,r,100,1000) for r in (20,60)}
        arrays[tag+'_stress']=s
    for tag,h in [('stat01',.1),('stat005',.05),('stat0025',.025)]:
        s=table(tag,'stress');out['stationary'][str(h)]={'t0_trace60_center_mX':float(s['trace60_center'][0]),'t0_trace60_face_mX':float(s['trace60_face'][0]),'mean_trace60_0_200_mX':window(s,60,0,200)['trace_mean_mX']}
    for tag,h in [('pulse01',.1),('pulse005',.05),('pulse0025',.025)]:
        s=table(tag,'stress');assert abs(s['tau'][-1]-20000)<1e-9,(tag,s['tau'][-1]); out['late'][str(h)]={str(r):window(s,r,10000,20000) for r in (20,60)}
        hist=table(tag,'history');out['replay'][str(h)]={'end_mXt':float(hist['tau'][-1]),'max_energy_plus_absorber_drift_mX':float(100*np.max(abs(hist['E']+hist['E_abs']-hist['E'][0]))),
            'max_charge_plus_absorber_drift':float(100*np.max(abs(hist['Q']+hist['Q_abs']-hist['Q'][0]))),'E20_end_mX':float(100*hist['E20'][-1]),'Q20_end':float(100*hist['Q20'][-1])}
        if h in (.05,.025):arrays[tag+'_stress']=s;arrays[tag+'_history']=hist
    late=table('late_fine0025','stress');late_hist=table('late_fine0025','history');prior_hist=table('pulse0025','history')
    assert late['tau'][-1]==1000
    out['late_restart']={'absolute_start_mXt':20000,'absolute_end_mXt':21000,
        'initial_core_E20_difference_mX':float(100*(late_hist['E20'][0]-prior_hist['E20'][-1])),
        'initial_core_Q20_difference':float(100*(late_hist['Q20'][0]-prior_hist['Q20'][-1])),
        'windows':{str(r):{'20000_21000_dt025':window(late,r,0,1000),'20000_21000_dt1':window(late[::4],r,0,1000),
                           '20100_21000_dt025':window(late,r,100,1000),'20100_21000_dt1':window(late[::4],r,100,1000)} for r in (20,60)}}
    arrays['late_fine0025_stress']=late;arrays['late_fine0025_history']=late_hist
    p=np.load(ROOT.parent/'grut_m1_expansion_clock_20261028'/'ringdown_h0025_history.npz')['history'];q=table('pulse0025','history');n=len(p)
    assert np.array_equal(p['tau'],q['tau'][:n]);out['replay']['parent_max_abs_difference_unscaled']={k:float(np.max(abs(p[k]-q[k][:n]))) for k in p.dtype.names}
    old=json.loads((ROOT.parent/'grut_m1_core_matching_20261029'/'late_h0025_analysis.json').read_text())
    bytime={int(row['mX_t_absolute']):row for row in old['core_checkpoints']}
    for at in (10000,15000,20000):
        item=q[np.argmin(abs(q['tau']-at))];o=bytime[at]
        out['replay'][f'parent_checkpoint_{at}_differences']={'E20_mX':float(100*item['E20']-o['E20_mX']),'Q20':float(100*item['Q20']-o['Q20'])}
    baseline=table('open_2000','history');nphysical=1.2141071422639381e-15
    for tag,r in [('reflect_80',80),('reflect_120',120),('reflect160_retry',160)]:
        d=table(tag,'history');assert d['tau'][-1]==2000;assert np.array_equal(d['tau'],baseline['tau'])
        at=lambda t:np.flatnonzero(d['tau']==t)[0]
        item={'R':r,'cell_density_over_physical_ncore':float(3/(4*np.pi*r**3)/nphysical),'snapshots':{},'max_total_energy_drift_fraction':float(np.max(abs(d['E']-d['E'][0]))/d['E'][0]),'max_total_charge_drift_fraction':float(np.max(abs(d['Q']-d['Q'][0]))/d['Q'][0])}
        for t in (100,400,1000,2000):
            i=at(t);item['snapshots'][str(t)]={'delta_E20_mX':float(100*(d['E20'][i]-baseline['E20'][i])),'delta_Q20':float(100*(d['Q20'][i]-baseline['Q20'][i]))}
        out['feedback'][tag]=item;arrays[tag+'_history']=d
    arrays['open_2000_history']=baseline
    radial=table('radial_mono300','history');axis=table('axis_mono300','history');assert np.array_equal(radial['tau'],axis['tau'])
    out['axis']['pure_monopole_max_differences_unscaled']={k:float(np.max(abs(axis[k]-radial[k]))) for k in ('E','Q','E20','Q20')}
    out['axis']['budget_controls']={}
    for tag in ('quad02','quad01','quad005','quad_l6','quad_nq16','quad01_t1000','quad005_t1000','quad_l6_t1000'):
        h=table(tag,'history');out['axis']['budget_controls'][tag]={'end':float(h['tau'][-1]),
            'max_relative_energy_plus_absorber_drift':float(np.max(abs(h['E']+h['Eabs']-h['E'][0]))/h['E'][0]),
            'max_relative_charge_plus_absorber_drift':float(np.max(abs(h['Q']+h['Qabs']-h['Q'][0]))/h['Q'][0])}
        arrays[tag+'_history']=h
    for path in sorted(ROOT.glob('axis_spec*.json')):
        z=json.loads(path.read_text());m=z['modes'];ns=sum(x['raw_neutral_energy'] for x in m); cs=sum(x['raw_charged_energy'] for x in m)
        out['axis'][path.stem]={'h':z['h'],'lmax':z['lmax'],'window':z['window_width'],'end':z['time_end'],'neutral_total_mX':ns,'charged_total_mX':cs,
            'neutral_quadrupole_fraction':next(x['raw_neutral_energy'] for x in m if x['ell']==2)/ns,
            'neutral_quadrupole_energy_mX':next(x['raw_neutral_energy'] for x in m if x['ell']==2),
            'neutral_quadrupole_positive_projection':next(x['species']['neutral']['positive_projection_number'] for x in m if x['ell']==2),
            'max_parseval_error':max(x['parseval_max_abs_error'] for x in m)}
    for path in sorted(ROOT.glob('axis_direction*.json')):
        z=json.loads(path.read_text());out['axis'][path.stem]={'synthetic_error':z['synthetic_hankel_split_max_absolute_error'],'k_min':z['k_min'],'gradient_proxy':z['number_weighted_angular_gradient_proxy']}
    audit=json.loads((ROOT.parent/'grut_m1_core_matching_20261029'/'SCENARIO_AUDIT.json').read_text());a=audit['energy_decompositions_at_scale_factor_10']['parent'];b=audit['energy_decompositions_at_scale_factor_10']['through_20000'];eref=audit['mass_reference_mX']
    out['ledger']={'held_drop_mX':a['held_remainder_mX']-b['held_remainder_mX'],'cold_increase_mX':(b['cold_excess_percent']-a['cold_excess_percent'])*eref/100,
                  'gas_kinetic_increase_mX':(b['gas_kinetic_percent']-a['gas_kinetic_percent'])*eref/100,'work_increase_mX':b['work_redshift_mX']-a['work_redshift_mX'],
                  'parent_full_percent':a['full_percent'],'through_20000_full_percent':b['full_percent']}
    np.savez_compressed(ROOT/'field_diagnostics.npz',**arrays)
    (ROOT/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
