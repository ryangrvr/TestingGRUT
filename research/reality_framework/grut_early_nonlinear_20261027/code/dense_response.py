"""Repair the sparse-sampling diagnostic, keeping old sealed files unchanged."""
from en_common import *

def main():
    exe=T.build();times=np.arange(0,240.00001,.025);ks=np.r_[.02,.1,.4,.8,1.2,1.6,2.,2.13,2.4,2.8,3.2,5.8]
    rows=[]
    for b,gam in [(3.,.5),(5.,1.)]:
        tag=f'dense_b{b:g}';bg,vec,ctrl=T.integrate(exe,tag,b,gam,ks,times)
        amp,dq,nx=T.responses(b,bg,vec);a=np.exp(bg[:,15])
        norms=[];qref=[]
        for v in bg:
            _,_,_,q,p=M.matrices(v,b);norms.append(apo(q[:2],p[:2]));qref.append(q)
        full=np.linalg.norm(vec[:,:,:2],axis=2)/np.array(norms)[:,None]
        late=(times>=64)&(times<=80);norm=np.trapezoid(nx[late],times[late])/(times[late][-1]-times[late][0])
        row=dict(beta=b,Gamma_over_H=gam,controls=ctrl,
            preparation_phase_003_amplitude_warning=T.crossings(times,a,.03*amp,ks,.1),
            preparation_phase_003_pericentre_unity=T.crossings(times,a,.03*amp,ks,1.),
            preparation_phase_003_orbit_vector_unity=T.crossings(times,a,.03*full,ks,1.),
            max_amplitude_through80=float(.03*np.max(abs(amp[times<=80]))),
            mean_current_dense64_80=float(norm),scope='Supplied Sigma preparation-phase derivative; single cosine modes. Pericentre and orbit-scale thresholds differ; no nonlinear evolution inferred.')
        rows.append(row);selected=np.array([int(np.argmin(abs(ks-k))) for k in [.02,2.,2.13,2.4]])
        np.savez_compressed(ROOT/(tag+'.npz'),times=times,selected_ks=ks[selected],
            selected_amplitude_per_radian=amp[:,selected],selected_orbit_vector_per_radian=full[:,selected],
            selected_charge_response_per_radian=dq[:,selected]/abs(norm),
            maximum_amplitude_per_radian=np.max(abs(amp),axis=1),maximum_orbit_vector_per_radian=np.max(full,axis=1),a=a)
        print('DENSE',b,row['preparation_phase_003_amplitude_warning'],row['max_amplitude_through80'],flush=True)
    v,_=M.prepared(3.);q=M.realvec(v[2:8:2]+1j*v[3:8:2]);seed=T.phase_seed(3.)
    save('DENSE_CORRECTIONS.json',dict(rows=rows,seed_definition=dict(
        dtheta_X_per_radian=(q[0]*seed[1]-q[1]*seed[0])/np.sum(q[:2]**2),
        dtheta_Sigma_per_radian=(q[4]*seed[5]-q[5]*seed[4])/np.sum(q[4:6]**2)),
        zero_k_note='No exponential auxiliary growth at exact k0; secular Jordan drift can remain. A full-system relative X/Sigma phase is not a pure gauge rotation.',
        physical_breakup_time=None))

if __name__=='__main__':main()
