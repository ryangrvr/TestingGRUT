"""Mode-resolved finite-window signed flux for the axisymmetric open run."""
from pathlib import Path
import argparse
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
MH=np.sqrt(6.)

def analyze(tag,h,width,lmax,bins=280):
    data=np.genfromtxt(ROOT/'.work'/(tag+'_trace.csv'),names=True,delimiter=',');t=data['tau'];dt=float(t[1]-t[0]);N=len(t)
    n=round(width/dt);assert n%4==0;hop=n//2;centers=np.arange(0,N,hop,dtype=int)
    if centers[-1]!=N-1:raise ValueError('endpoint must be a window center')
    off=np.arange(-n//2,n//2);win=np.cos(np.pi*off/n);cover=np.zeros(N)
    for c in centers:
        idx=c+off;ok=(idx>=0)&(idx<N);cover[idx[ok]]+=win[ok]**2
    if np.max(abs(cover-1))>1e-12:raise ValueError('window partition failed')
    nf=4*n;freq=2*np.pi*np.fft.fftfreq(nf,dt);period=nf*dt
    radius=(int(np.floor(60/h))+.5)*h;edges=np.linspace(1,8,bins+1);mid=(edges[1:]+edges[:-1])/2
    results=[];arrays={'energy_edges':edges,'energy_centers':mid,'time_centers':t[centers]}
    for ell in range(0,lmax+1,2):
        area=4*np.pi*radius**2/(.01*(2*ell+1));species={s:np.zeros((len(centers),bins)) for s in ('plus','minus','neutral')}
        allflux=np.zeros(3)
        for j,c in enumerate(centers):
            idx=c+off;ok=(idx>=0)&(idx<N);window=np.zeros(n);window[ok]=win[ok]
            def field(name):
                x=np.zeros(n);x[ok]=data[name+'_l'+str(ell)][idx[ok]]*window[ok];return dt*np.fft.fft(x,n=nf)
            Z=field('x')+1j*field('y');DZ=field('dx')+1j*field('dy');P=field('px')+1j*field('py')
            H=field('s_minus1');DH=field('ds');PH=field('ps')
            Q=area/period*np.imag(np.conj(Z)*DZ);CE=-area/period*np.real(np.conj(P)*DZ);HE=-area/period*np.real(np.conj(PH)*DH)
            allflux+=np.array([Q.sum(),CE.sum(),HE.sum()])
            for name,mask in [('plus',freq<-1),('minus',freq>1),('neutral',abs(freq)>MH)]:
                e=abs(freq[mask]);weights=Q[mask] if name=='plus' else -Q[mask] if name=='minus' else HE[mask]/e
                species[name][j]=np.histogram(e,bins=edges,weights=weights)[0]
        Z=data['x_l'+str(ell)]+1j*data['y_l'+str(ell)];DZ=data['dx_l'+str(ell)]+1j*data['dy_l'+str(ell)];P=data['px_l'+str(ell)]+1j*data['py_l'+str(ell)]
        raw=np.array([dt*area*np.imag(np.conj(Z)*DZ).sum(),-dt*area*np.real(np.conj(P)*DZ).sum(),
                      -dt*area*(data['ps_l'+str(ell)]*data['ds_l'+str(ell)]).sum()])
        rr={'ell':ell,'raw_charge':float(raw[0]),'raw_charged_energy':float(raw[1]),'raw_neutral_energy':float(raw[2]),
            'parseval_max_abs_error':float(np.max(abs(allflux-raw))),'species':{}}
        for name,a in species.items():
            b=np.maximum(a,0);arrays[f'{name}_l{ell}_signed']=a;arrays[f'{name}_l{ell}_positive']=b
            rr['species'][name]={'signed_number':float(a.sum()),'positive_projection_number':float(b.sum()),
                                 'negative_bin_weight':float(np.maximum(-a,0).sum()),'signed_energy_mX':float((a*mid).sum())}
        results.append(rr)
    return {'tag':tag,'h':h,'lmax':lmax,'window_width':width,'probe_radius':radius,'time_end':float(t[-1]),
            'window_partition_max_error':float(np.max(abs(cover-1))),'modes':results,
            'scope':'Finite-time angular-mode classical flux with positive spectral projection. Angular modes are coherent waves, not a derived two-particle angle distribution; near-threshold and window artifacts remain.'},arrays,data

def main():
    p=argparse.ArgumentParser();p.add_argument('--tag',required=True);p.add_argument('--h',type=float,required=True)
    p.add_argument('--lmax',type=int,default=4);p.add_argument('--width',type=float,default=200.);p.add_argument('--out',required=True);a=p.parse_args()
    d,arrays,trace=analyze(a.tag,a.h,a.width,a.lmax);arrays['trace']=trace
    np.savez_compressed(ROOT/(a.out+'.npz'),**arrays)
    (ROOT/(a.out+'.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');print(json.dumps(d,indent=2))

if __name__=='__main__':main()
