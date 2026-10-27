"""Export measured classical field covariances, never an assumed particle gas."""
from en_common import *

def main():
    runs=read('RUNS.json')['cases'];out={};checks=[]
    for tag,r in runs.items():
        N,L,s,raw=snapshot(r['snapshots'][-1]);q=np.array(raw[:3]);p=np.array(raw[3:]);k2=symbol(N,L)
        qh=rfftn(q,axes=(-3,-2,-1),norm='forward');ph=rfftn(p,axes=(-3,-2,-1),norm='forward')
        v=np.array([qh[0],qh[1],ph[0],ph[1],qh[2],ph[2]])
        f=np.fft.fftfreq(N)*2*np.pi*N/L;last=np.fft.rfftfreq(N)*2*np.pi*N/L
        physical=np.sqrt(f[:,None,None]**2+f[None,:,None]**2+last[None,None,:]**2)
        edges=np.linspace(0,physical.max()*1.00001,51);bins=np.clip(np.digitize(physical.ravel(),edges)-1,0,49)
        weights=np.full(N//2+1,2.);weights[[0,-1]]=1.;weights=np.broadcast_to(weights[None,None,:],k2.shape).copy();weights[0,0,0]=0.
        cov=np.zeros((50,6,6));covk=cov.copy()
        for i in range(6):
            for j in range(i,6):
                val=weights*np.real(v[i].conj()*v[j]);c=np.bincount(bins,weights=val.ravel(),minlength=50)
                ck=np.bincount(bins,weights=(val*k2).ravel(),minlength=50)
                cov[:,i,j]=cov[:,j,i]=c;covk[:,i,j]=covk[:,j,i]=ck
        meanq=q.mean(axis=(1,2,3));meanp=p.mean(axis=(1,2,3));Qmean=meanq[0]*meanp[1]-meanq[1]*meanp[0]
        Q=float(np.mean(q[0]*p[1]-q[1]*p[0]));Qexport=Qmean+cov[:,0,3].sum()-cov[:,1,2].sum()
        Ekin=.5*(cov[:,2,2]+cov[:,3,3]+cov[:,5,5]).sum()+.5*np.sum(meanp**2)
        Egrad=.5*(covk[:,0,0]+covk[:,1,1]+covk[:,4,4]).sum()
        r2=q[0]**2+q[1]**2;U=np.mean(LX/4*r2*r2+G2/2*r2*q[2]**2+LP/4*q[2]**4)
        E=energy(q,p,L);qe=abs(Qexport/Q-1);ee=abs((Ekin+Egrad+U)/E-1)
        checks.append(dict(tag=tag,charge_relative=qe,energy_relative=ee))
        out[tag+'_covariance']=cov;out[tag+'_k2_covariance']=covk;out[tag+'_k_edges']=edges
        out[tag+'_mean_q']=meanq;out[tag+'_mean_p']=meanp;out[tag+'_potential_density']=U;out[tag+'_L']=L
    np.savez_compressed(ROOT/'SPECTRAL_EXPORT.npz',**out)
    save('EXPORT_CHECKS.json',dict(passed=max(max(x['charge_relative'],x['energy_relative']) for x in checks)<1e-12,rows=checks,
        component_order=['chi_X_real','chi_X_imag','pi_X_real','pi_X_imag','chi_phi','pi_phi'],
        convention='Forward FFT; half-space weights account for both signs. Binned real symmetric classical equal-time covariance, with k0 removed and mean stored separately. k2_covariance uses the exact FD4 symbol. Potential density stays separate.',
        physical_distribution=None,scope='Measured field input for a future spectral/kinetic match. Contains no quantum commutator, renormalization, pole measurement or particle assignment.'))

if __name__=='__main__':main()
