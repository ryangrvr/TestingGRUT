"""Illustrative re-injection sensitivity within the inherited kinetic model.

Inject 20% of the recorded, projected 10000--20000 spectrum after the a=10
checkpoint, once at a=10 and once at a later chosen scale factor. Match the ledger by subtracting
its *instantaneous* core-plus-gas energy cost from an explicitly held residual.
This is a controlled counterexample, not a prediction of future core emission.
"""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
RF=ROOT.parent
sys.path.insert(0,str(RF/'grut_m1_core_matching_20261029'/'code'))
import late_kinetics
from hybrid import Hybrid
import expansion

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--late',type=float,default=100);ap.add_argument('--final',type=float,default=1000);ap.add_argument('--out',default='timing_sensitivity.json');args=ap.parse_args()
    if not 10<args.late<args.final:raise ValueError('require 10 < late < final')
    t0=time.monotonic(); source_dir=RF/'grut_m1_core_matching_20261029'
    z=np.load(source_dir/'late_kin_h0025_64.npz');e=z['E'];n=len(e);y0=z['Y'][:,-1]
    if not np.allclose(e,expansion.grid(0,8.,'translated',64),rtol=0,atol=0):raise ValueError('grid mismatch')
    source,audit=late_kinetics.append_late(e,np.zeros(2*n),'h0025');extra=.2*source
    previous=json.loads((source_dir/'late_kin_h0025_64.json').read_text())
    held=previous['held_unassigned_recomputed_mX'];epsilon=expansion.clock()['H0_mX2_over_ncore0']
    model=Hybrid(e,'full',8,epsilon,.00025,True)
    def run(a0,a1,y,points):
        sol=solve_ivp(model.rhs,(np.log(a0),np.log(a1)),y,method='BDF',jac=model.jac,rtol=1e-6,atol=1e-9,t_eval=np.log(points))
        if not sol.success:raise RuntimeError(sol.message)
        ledger=[]
        for x,yy in zip(sol.t,sol.y.T):
            o=model.observables_hybrid(x,yy);ledger.append(o['total_stationary_plus_gas_plus_work_mX']-o['interpolated_generator_energy_defect_mX'])
        return sol,model.observables_hybrid(sol.t[-1],sol.y[:,-1]),float(np.max(np.abs(np.array(ledger)-ledger[0])))
    print('built',round(time.monotonic()-t0,1),'nodes',len(e),'held',held,flush=True)
    base,basefinal,base_drift=run(10,args.final,y0,[10,args.late,args.final]);print('baseline',round(time.monotonic()-t0,1),basefinal['Qcore'],flush=True)
    rows=[]
    for j,ainj in [(0,10),(1,args.late)]:
        before=base.y[:,j].copy();after=before.copy();after[:2*n]+=extra
        qb=model.terms(before)[0];qa=model.terms(after)[0]
        source_gas=float(e@(extra[:n]+extra[n:]));cost=float(model.mass(qa)-model.mass(qb)+source_gas)
        if not (0<cost<held):raise ValueError(f'cost {cost} not covered by held residual {held}')
        sol,final,drift=run(ainj,args.final,after,[ainj,args.final]);print('branch',ainj,round(time.monotonic()-t0,1),final['Qcore'],flush=True)
        cold=lambda v:float(model.mass(v['Qcore']))+v['Nplus']+v['Nminus']
        dc=cold(final)-cold(basefinal);dw=final['work_redshift_mX']-basefinal['work_redshift_mX']
        dk=final['gas_kinetic_mX']-basefinal['gas_kinetic_mX'];dd=final['interpolated_generator_energy_defect_mX']-basefinal['interpolated_generator_energy_defect_mX']
        rows.append({'a_inject':ainj,'delta_Q_at_injection':float(qa-qb),'injected_gas_energy_mX':source_gas,'held_residual_cost_mX':cost,
                     'final':final,'delta_cold_stationary_plus_rest_mX':dc,'delta_cold_plus_held_mX':dc-cost,
                     'delta_work_mX':dw,'delta_gas_kinetic_mX':dk,'delta_generator_defect_mX':dd,
                     'ledger_identity_residual_mX':dc-cost+dw+dk-dd,'saved_point_ledger_drift_mX':drift,
                     'minimum_saved_population':float(sol.y[:2*n].min()),'solver':{'nfev':sol.nfev,'njev':sol.njev,'nlu':sol.nlu}})
    out={'scope':'Illustrative instantaneous source re-injected after a=10 with the same projected finite-time spectrum scaled to 20%; no claim that a real core emits this, that held residual is available excitation, or that constant-g radiation expansion is valid throughout the interval. Frozen prescribed source; homogeneous mixing and quasistatic profile approximations inherited. High-a numerical diffusion and finite energy support persist.',
         'parent_source':'late_kin_h0025_64','source_fraction':.2,'held_residual_mX':held,'a_final':args.final,
         'late_source_audit':audit,'baseline_final':basefinal,'baseline_ledger_drift_mX':base_drift,'branches':rows,'seconds':time.monotonic()-t0}
    (ROOT/args.out).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
