"""Restart and same-lattice stationary-core diagnostics for late field runs."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
from scipy.optimize import curve_fit

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_m1_expansion_clock_20261028'
sys.path.insert(0,str(PARENT/'code'))
from source_ringdown import stationary_lattice

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--tag',choices=['h005','h0025'],required=True)
    p.add_argument('--old-history',type=Path);p.add_argument('--old-final',type=Path);a=p.parse_args()
    fine=a.tag=='h0025';h=.025 if fine else .05
    oldhist=a.old_history or (PARENT/'.work'/f'ringdown_{a.tag}_history.csv')
    oldfinal=a.old_final or (PARENT/'.work'/f'ringdown_{a.tag}_final.csv')
    prefix=ROOT/'.work'/('late_'+a.tag)
    before=np.genfromtxt(oldhist,names=True,delimiter=',');after=np.genfromtxt(str(prefix)+'_history.csv',names=True,delimiter=',')
    assert abs(before['tau'][-1]-10000)<1e-9 and abs(after['tau'][-1]-10000)<1e-9
    restart={k:float(after[k][0]-before[k][-1]) for k in ('E','Q','E12','Q12','E20','Q20','E30','Q30')}
    earlier=[]
    for dt in (4000,6000,8000):
        r=before[int(np.argmin(abs(before['tau']-dt)))];q=float(r['Q20']/.01);e=float(r['E20']/.01)
        stationary=stationary_lattice(q,h,200.)
        earlier.append({'mX_t_absolute':float(r['tau']),'Q20':q,'E20_mX':e,
                        'stationary_inside20_mX':stationary['energy_inside20_mX'],
                        'excess_inside20_over_same_grid_stationary_mX':e-stationary['energy_inside20_mX'],
                        'stationary_equation_max_residual':stationary['max_stationary_equation_residual']})
    rows=[]
    for dt in (0,1000,5000,7500,10000):
        idx=int(np.argmin(abs(after['tau']-dt)));r=after[idx]
        q=float(r['Q20']/.01);e=float(r['E20']/.01);stationary=stationary_lattice(q,h,200.)
        rows.append({'mX_t_absolute':10000+float(r['tau']),'Q20':q,'E20_mX':e,
                     'stationary_inside20_mX':stationary['energy_inside20_mX'],
                     'excess_inside20_over_same_grid_stationary_mX':e-stationary['energy_inside20_mX'],
                     'stationary_equation_max_residual':stationary['max_stationary_equation_residual'],
                     'additional_probe_energy_mX':float(r['Eflux_probe']/.01),
                     'additional_probe_charge':float(r['Qflux_probe']/.01)})
    E=(after['E']+after['E_abs'])/.01;Q=(after['Q']+after['Q_abs'])/.01
    full=earlier+rows;t=np.array([x['mX_t_absolute'] for x in full]);y=np.array([x['excess_inside20_over_same_grid_stationary_mX'] for x in full]);u=t-t[0]
    def no_floor(x,A,tau):return A*np.exp(-x/tau)
    def with_floor(x,A,tau,c):return c+A*np.exp(-x/tau)
    ab,_=curve_fit(no_floor,u,y,p0=(y[0],20000.),bounds=([0.,1.],[1e7,1e7]),maxfev=20000)
    abc,_=curve_fit(with_floor,u,y,p0=(y[0]*.8,15000.,y[-1]*.3),bounds=([0.,1.,0.],[1e7,1e7,1e7]),maxfev=20000)
    fits={'zero_floor':{'amplitude_at_4000_mX':float(ab[0]),'decay_time_mX_inverse':float(ab[1]),
                        'max_abs_residual_mX':float(np.max(abs(y-no_floor(u,*ab))))},
          'free_nonnegative_floor':{'amplitude_at_4000_mX':float(abc[0]),'decay_time_mX_inverse':float(abc[1]),
                                    'floor_mX':float(abc[2]),'max_abs_residual_mX':float(np.max(abs(y-with_floor(u,*abc))))},
          'scope':'Empirical fits to eight sparse inside-20 energy checkpoints, not a complex quasinormal-mode measurement or a certified long-time asymptote.'}
    source=[]
    for width in (200,400,800):
        d=json.loads((ROOT/f'late_{a.tag}_t{width}.json').read_text());item={'width':width,'raw_flux':d['raw_flux']}
        for species in ('plus','minus','neutral'):
            s=d['species'][species];item[species]={k:s[k] for k in ('signed_number','positive_projection_number','negative_bin_weight','signed_energy_mX','split_k_gt_0.1')}
        source.append(item)
    out={'radial_step':h,'restart_time_mX_t':10000,'end_time_mX_t':20000,
         'input_final_sha256':sha(oldfinal),'new_trace_sha256':sha(Path(str(prefix)+'_trace.csv')),
         'restart_differences_in_unscaled_engine_units':restart,
         'checks':{'segment_energy_plus_absorber_max_abs_drift_mX':float(np.max(abs(E-E[0]))),
                   'segment_charge_plus_absorber_max_abs_drift':float(np.max(abs(Q-Q[0]))),
                   'first_trace_time':float(after['tau'][0]),'last_trace_time':float(after['tau'][-1])},
         'earlier_core_checkpoints':earlier,'core_checkpoints':rows,'phenomenological_fits':fits,'late_window_sources':source,
         'scope':'Isolated radial field run restarted from the prior final state; distant absorber, no interacting gas, imposed h decay or expansion. Core20 excess and earlier unassigned ledger residual are not automatically the same energy.'}
    (ROOT/f'late_{a.tag}_analysis.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
