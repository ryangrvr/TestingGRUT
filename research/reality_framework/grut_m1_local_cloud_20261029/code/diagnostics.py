"""Time-source, matching-radius and first-passage summary."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid
from local_rates import pulse_series

ROOT=Path(__file__).resolve().parents[1]

def main():
    windows=[]
    for width in (200,400,800):
        z=np.load(ROOT/f'source_t{width}.npz');t,e,p,m=pulse_series(z,5.);de=e[1]-e[0]
        row={'width_mX_inverse':width}
        for s,flux in [('plus',p.sum(axis=1)*de),('minus',m.sum(axis=1))]:
            cum=cumulative_trapezoid(flux,x=t,initial=0);total=cum[-1]
            row[s]={'number':float(total),'t10_t50_t90':np.interp([.1,.5,.9],cum/total,t).tolist(),
                    'fraction_after_1000':float(1-np.interp(1000,t,cum)/total)}
        windows.append(row)
    h=[]
    for r in (40,60,80):
        d=json.loads((ROOT/f'source_h005_r{r}.json').read_text());x=json.loads((ROOT/f'continuous_h005_r{r}.json').read_text())['result']
        h.append({'radius':r,'actual_radius':d['probe_actual_mX_inverse'],'source_counts':{s:d['species'][s]['positive_projection_number'] for s in ('plus','minus','neutral')},
                  'raw_charge':d['raw_flux']['charge'],'raw_total_energy_mX':d['raw_flux']['charged_energy_mX']+d['raw_flux']['neutral_energy_mX'],
                  'collinear_primary_NWA_events_outside_surface':x['expected_unattenuated_events'],'trace_sha256':d['trace_sha256']})
    hist=[]
    for r in (40,80):
        p=ROOT/'.work'/f'radial_r{r}_history.csv'
        a=np.genfromtxt(p,names=True,delimiter=',');hist.append(a)
    errors={k:float(np.max(abs(hist[0][k]-hist[1][k]))) for k in ('E','Q','E_abs','Q_abs','E20','Q20')}
    out={'time_window_controls':windows,'matching_radius_controls':h,'same_dynamics_history_max_absolute_differences':errors,
         'scope':'Windowed positive projections and direct-beam NWA sensitivity. Radius changes include different finite-time matching and phase-space approximations; no interacting cloud or inferred angular distribution.'}
    (ROOT/'SOURCE_DIAGNOSTICS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
