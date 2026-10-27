"""Same-author independent identities and necessary-condition controls."""
from cw_common import *

def main():
    rng=np.random.default_rng(217061);N=7;L=5.
    q=rng.normal(size=(3,N,N,N));q[2]+=8
    p=rng.normal(size=q.shape);f=ledger(q,p,L);mask=rng.random((N,N,N))>.4
    rr=region(q,p,f,mask,L);groups={}
    error=max(rr['kinetic_identity_relative'],rr['rotation_fit_identity_relative'])
    om=rr['omega_fit'];dv=(L/N)**3
    costs=[]
    for offset in [-1,-.1,.1,1]:
        trial=om+offset
        cost=.5*np.sum(((p[0]+trial*q[1])**2+(p[1]-trial*q[0])**2)[mask])*dv
        costs.append(cost-(rr['kinetic_radial']+rr['kinetic_angular_dispersion']))
    qnode=q.copy();qnode[:2,0]=0
    fn=ledger(qnode,p,L);rn=region(qnode,p,fn,mask,L)
    node_error=max(rn['kinetic_identity_relative'],rn['rotation_fit_identity_relative'])
    groups['kinetic_identity_and_unique_rotation_minimum']=dict(passed=max(error,node_error)<1e-13 and min(costs)>0 and rn['kinetic_at_exact_nodes']>0,
        relative_error=error,exact_node_relative_error=node_error,trial_excesses=costs)

    expected_X=f['dEX_flux']+f['work_to_X'];expected_phi=f['dEphi_flux']-f['work_to_X']
    charge_direct=q[0]*f['force'][1]-q[1]*f['force'][0]
    charge_error=float(np.max(abs(charge_direct-f['dQ']))/max(1,np.max(abs(charge_direct))))
    derivatives=[]
    for h in [1e-5,2e-6,4e-7]:
        plus=ledger(q+h*p,p+h*f['force'],L);minus=ledger(q-h*p,p-h*f['force'],L)
        errors=[]
        for key,want in [('eX',expected_X),('ephi',expected_phi)]:
            got=(plus[key]-minus[key])/(2*h)
            errors.append(float(np.max(abs(got-want))/max(1,np.max(abs(want)))))
        derivatives.append(dict(step=h,X=errors[0],phi=errors[1]))
    err=min(max(r['X'],r['phi']) for r in derivatives)
    groups['local_energy_charge_equations']=dict(passed=err<1e-8 and charge_error<1e-13,
        energy_difference_checks=derivatives,charge_relative=charge_error)

    # Separate boundary traversal rather than the divergence assembly in ledger.
    boundary=np.zeros(3);d=L/N
    for i in np.ndindex(N,N,N):
        for ax in range(3):
            for step,coef in [(1,4/3/d**2),(2,-1/12/d**2)]:
                k=list(i);k[ax]=(k[ax]+step)%N;k=tuple(k)
                weight=float(mask[i])-float(mask[k])
                a=q[(slice(None),)+i];b=q[(slice(None),)+k]
                u=p[(slice(None),)+i];v=p[(slice(None),)+k]
                boundary[0]+=.5*coef*np.dot((b-a)[:2],(u+v)[:2])*weight
                boundary[1]+=.5*coef*(b-a)[2]*(u+v)[2]*weight
                boundary[2]+=coef*(a[0]*b[1]-a[1]*b[0])*weight
    want=np.array([f['dEX_flux'][mask].sum(),f['dEphi_flux'][mask].sum(),f['dQ'][mask].sum()])
    be=float(np.max(abs(boundary-want))/max(1,np.max(abs(want))))
    global_derivatives=[float(np.mean(f[k])) for k in ['dEX_flux','dEphi_flux','dQ']]
    groups['independent_boundary_bonds']=dict(passed=be<1e-13 and max(abs(x) for x in global_derivatives)<1e-11,
        boundary_relative=be,global_flux_divergence_means=global_derivatives)

    # A stationary homogeneous charged control: passing these conditions is NOT binding.
    r=3.;z=np.sqrt(C.NU**2-C.G**2/C.LP*r*r);omega=np.sqrt(C.G**2*z*z+C.LX*r*r)
    qs=np.zeros_like(q);ps=np.zeros_like(q);qs[0]=r;qs[2]=z;ps[1]=omega*r
    fs=ledger(qs,ps,L);cr=region(qs,ps,fs,np.ones_like(mask),L)
    ce=max(cr['kinetic_excess_fraction'],cr['stationary_X_force_residual_over_m2'],cr['phi_velocity_RMS_over_mv'],cr['phi_acceleration_RMS_over_m2v'])
    groups['stationary_rotation_control_not_bound_object']=dict(passed=ce<1e-13,maximum_residual=ce,
        localized=False,bound_object_certified=False)

    old=A.local_energy(q,p,L);total=f['eX']+f['ephi']
    equal=float(abs(np.mean(old)-np.mean(total))/max(1,abs(np.mean(old))))
    fft=float(abs(P.energy(q,p,L)-np.mean(total))/max(1,abs(np.mean(total))))
    groups['symmetric_forward_FFT_global_energy']=dict(passed=max(equal,fft)<1e-13,
        forward_relative=equal,FFT_relative=fft,
        local_allocation_difference_nonzero=bool(np.max(abs(old-total))>0))

    qc=q.copy();pc=p.copy();qc[1]*=-1;pc[1]*=-1
    fc=ledger(qc,pc,L);rc=region(qc,pc,fc,mask,L)
    qz=q.copy();pz=p.copy();qz[2]*=-1;pz[2]*=-1
    fz=ledger(qz,pz,L);rz=region(qz,pz,fz,mask,L)
    checks=[abs(rc['Q']+rr['Q']),abs(rc['omega_fit']+rr['omega_fit']),
            abs(rc['outward_charge_flux']+rr['outward_charge_flux'])]
    for k in ['E_total','kinetic_radial','kinetic_angular_dispersion','outward_energy_flux_X','reversible_work_phi_to_X']:
        checks += [abs(rc[k]-rr[k]),abs(rz[k]-rr[k])]
    symmetry=float(max(checks)/max(1,rr['E_total']))
    groups['conjugation_and_neutral_parity']=dict(passed=symmetry<1e-13,scaled_error=symmetry)
    save('VERIFICATION_RESULTS.json',dict(groups=groups,independent_reviewer_ran=False,physical_f_sym=None))
    print((ROOT/'VERIFICATION_RESULTS.json').read_text())
    assert all(g['passed'] for g in groups.values())

if __name__=='__main__':main()
