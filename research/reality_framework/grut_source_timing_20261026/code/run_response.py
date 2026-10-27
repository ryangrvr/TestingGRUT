"""Collisionless, conditional response past the earlier tau80 endpoint."""
from timing_common import *
from concurrent.futures import ThreadPoolExecutor,as_completed

def summary(tag,beta,decay,e,ks,times,bg,vec,control):
    amp,dq,nx=responses(beta,bg,vec);a=np.exp(bg[:,15]);window=(times>=64)&(times<=80)
    norm=float(np.trapezoid(nx[window],times[window])/(times[window][-1]-times[window][0]));charge=dq/abs(norm)
    brackets=[]
    for seed in [.03,1e-6,1e-9]:
        brackets.append(dict(phase_seed_radians=seed,
            amplitude_10percent=crossings(times,a,seed*amp,ks,.1),
            amplitude_unity=crossings(times,a,seed*amp,ks,1.),
            local_charge_unity=crossings(times,a,seed*charge,ks,1.)))
    C=a*a*bg[:,16];last=times>=.9*times[-1]
    row=dict(tag=tag,beta=beta,Gamma_over_H=decay,e=e,ks=ks,controls=control,
        late_NX_normalization_through80=norm,final_homogeneous_NX=float(nx[-1]),
        homogeneous_NX_final_change_from_norm=float(nx[-1]/norm-1),
        final_a=float(a[-1]),late_C=float(C[-1]),late_C_relative_range=float(np.ptp(C[last])/abs(C[-1])),
        matched_radiation_tau_offset=float(times[-1]-a[-1]**2/(2*C[-1])),crossings=brackets,
        physical_breakup_time=None,physical_initial_covariance=None,physical_f_sym=None,
        scope='Linear responses are extrapolated mathematical solutions; each supplied seed loses small-perturbation control at its own crossing. No nonlinear continuation or actual decay/collision backreaction.')
    selected=np.array([int(np.argmin(abs(ks-k))) for k in [.02,2.13,5.8]])
    np.savez_compressed(ROOT/(tag+'.npz'),times=times,ks=ks,bg=bg,amplitude_per_radian=amp,
                        charge_fraction_per_radian=charge,selected_vectors=vec[:,selected],selected_vector_ks=ks[selected])
    return row

def main():
    exe=build();early=np.load(MODE/'MODE_PROFILES.npz')['times']
    times=np.unique(np.r_[early,np.arange(80.,800.00001,.5)])
    jobs=[('quiet3',3.,.5,1e-7,KS),('quiet5',5.,1.,1e-7,KS),
          ('gauge3',3.,.5,1e-6,np.array([.02,2.13,5.8])),('gauge5',5.,1.,1e-6,np.array([.02,2.13,5.8]))]
    rows=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        tasks={pool.submit(integrate,exe,tag,beta,decay,ks,times,e): (tag,beta,decay,e,ks) for tag,beta,decay,e,ks in jobs}
        for future in as_completed(tasks):
            tag,beta,decay,e,ks=tasks[future];bg,vec,control=future.result()
            row=summary(tag,beta,decay,e,ks,times,bg,vec,control);rows.append(row)
            save('RESPONSE_RESULTS.json',dict(rows=sorted(rows,key=lambda r:r['tag']),
                physical_f_sym=None,physical_breakup_time=None,source_normalization_recomputed=False))
            print('RESPONSE',tag,'seconds',control['seconds'],'C',row['late_C'],'crossings',row['crossings'],flush=True)

if __name__=='__main__':main()
