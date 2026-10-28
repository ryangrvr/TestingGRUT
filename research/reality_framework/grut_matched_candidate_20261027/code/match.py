"""Recompute a connected conditional candidate using the existing reviewed engine.

No old file or old output is mutated. H^(3/2) is checked by actual ODE runs.
"""
from pathlib import Path
import os,sys,json,math,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RF=ROOT.parent
WORK=Path(os.environ.get('GRUT_MATCHED_WORKDIR','/workspace/scratch/0c6d6cf183c6/matched_candidate_runs'))
WORK.mkdir(parents=True,exist_ok=True)
os.environ['GRUT_REHEATING_WORKDIR']=str(WORK)
sys.path.insert(0,str(RF/'grut_reheating_preparation_20261026/code'))
import reheating as old


def save(name,x):
    def cast(v):
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError(type(v).__name__)
    (ROOT/name).write_text(json.dumps(x,indent=2,allow_nan=False,default=cast)+'\n')


def run(exe,tag,beta,gamma,scale=1.,delta=0.,tol=2e-11,end=800.,samples=16000):
    original_h=old.H0
    old.H0=original_h*scale
    try:
        result,arr=old.run(exe,tag,beta,gamma,eta2=(old.H0/old.MRED)**2,
                           delta=delta,end=end,tol=tol,samples=samples)
    finally:
        old.H0=original_h
    result['uniform_scale']=scale
    result['absolute_yield_magnitude']=abs(result['yield_over_YQ'])*old.INP['benchmark']['YQ']
    result['physical_H_GeV']=original_h*scale
    # Keep actual source-transfer and late homogeneous drift distinct.
    for t in (20.,50.,80.,200.,800.):
        i=int(np.argmin(abs(arr[:,0]-t)))
        nx=arr[i,6];nl=arr[i,7]
        result.setdefault('charge_clock',[]).append({'tau':float(arr[i,0]),'NX':float(nx),'Nlight':float(nl)})
    mask=arr[:,0]>=640
    result['radiation_C']=float(np.mean(arr[mask,1]**2*arr[mask,2]))
    result['radiation_C_range_fraction']=float(np.ptp(arr[mask,1]**2*arr[mask,2])/result['radiation_C'])
    return result,arr


def main():
    exe=old.build();rows=[];arrays={}
    inp=old.INP;benchmark=inp['benchmark'];light=inp['light_sector']
    D=benchmark['DM_per_visible_entropy_GeV'];epq=light['mX_GeV']*benchmark['original_E_over_Qm']
    required=D/epq
    for beta,gamma in ((3.,.5),(5.,1.)):
        tag=f'b{beta:g}_g{gamma:g}'
        base,arr=run(exe,tag+'_baseline',beta,gamma)
        rows.append(base)
        factor=(required/base['absolute_yield_magnitude'])**(2/3)
        new,na=run(exe,tag+'_matched',beta,gamma,factor)
        rows.append(new)
        new['target_yield_magnitude']=required
        new['target_mass_ratio']=new['absolute_yield_magnitude']*epq/D
        print('MATCH',tag,'scale',factor,'mass_ratio',new['target_mass_ratio'],flush=True)
        # Declared compact export: early complete state and selected later clocks.
        for suffix,data in (('baseline',arr),('matched',na)):
            keep=(data[:,0]<=30.000001) | (np.arange(len(data))%100==0)
            arrays[tag+'_'+suffix]=data[keep]
    central=next(x for x in rows if x['tag']=='b3_g0.5_matched')
    for scale in (.5,2.):
        r,_=run(exe,'scale_'+str(scale),3.,.5,scale,samples=4000)
        rows.append(r)
    for delta in (-.001,.001):
        r,_=run(exe,'phase_'+str(delta),3.,.5,central['uniform_scale'],delta=delta,samples=8000)
        rows.append(r)
    tight,_=run(exe,'matched_tighter',3.,.5,central['uniform_scale'],tol=2e-12,samples=8000)
    rows.append(tight)
    np.savez_compressed(ROOT/'BACKGROUND_CONTROLS.npz',**arrays)
    neg=next(x for x in rows if x['tag']=='phase_-0.001')['absolute_yield_magnitude']
    pos=next(x for x in rows if x['tag']=='phase_0.001')['absolute_yield_magnitude']
    slope=(math.log(pos)-math.log(neg))/.002
    # Frozen standard uncorrelated-CDI benchmark; no current combined fit.
    beta_iso,As=.038,2.1e-9
    allowed=beta_iso/(1-beta_iso)*As
    save('MATCHED_RESULTS.json',{'scope':'homogeneous physical-reheating estimator and cold final-population matching hypothesis',
        'target_YQ_for_all_reference_objects_no_pairs':required,
        'target_DM_per_visible_entropy_GeV':D,'object_energy_per_charge_GeV':epq,
        'rows':rows,'columns':old.COLS+[f'state_{i}' for i in range(19)],
        'phase_response':{'definition':'derivative of ln absolute homogeneous charge yield versus supplied Sigma preparation phase',
                          'dlnY_dphase':slope,'phase_step_rad':.001,
                          'uncorrelated_CDI_beta_benchmark':beta_iso,'curvature_power_benchmark':As,
                          'allowed_charge_isocurvature_power':allowed,
                          'required_phase_sqrt_power_if_linear_transfer':math.sqrt(allowed)/abs(slope) if slope else None,
                          'phase_holder_constructed':False},
        'actual_final_object_fraction':None,'actual_final_symmetric_fraction':None,
        'cosmic_state_derived':False,'vacuum_density_derived':False})


if __name__=='__main__':main()
