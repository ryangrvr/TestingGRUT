"""Partition-of-unity time-frequency analysis of the stored radial probe.

This is an operational semiclassical source diagnostic. Frequency-local flux
is signed and window-dependent; positive-bin projection is only a sensitivity
prescription for a later kinetic estimate, not a derived particle density.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
G2=.01
MH=np.sqrt(6.)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--trace',required=True)
    p.add_argument('--width',type=float,default=400.)
    p.add_argument('--bins',type=int,default=280)
    p.add_argument('--emax',type=float,default=8.)
    p.add_argument('--h',type=float,default=.025)
    p.add_argument('--probe',type=float,default=60.)
    p.add_argument('--tag',required=True)
    args=p.parse_args()
    tracepath=Path(args.trace);data=np.genfromtxt(tracepath,names=True,delimiter=',')
    t=data['tau'];dt=float(t[1]-t[0]);N=len(t);T=t[-1]
    assert np.max(abs(np.diff(t)-dt))<1e-8
    n=int(round(args.width/dt));assert n%4==0
    hop=n//2;center=np.arange(0,N,hop,dtype=int)
    if center[-1]!=N-1:raise ValueError('Trace endpoint must be a window center')
    off=np.arange(-n//2,n//2);win=np.cos(np.pi*off/n)
    cover=np.zeros(N)
    for c in center:
        idx=c+off;ok=(idx>=0)&(idx<N);cover[idx[ok]]+=win[ok]**2
    covererr=float(np.max(abs(cover-1)))
    assert covererr<1e-12
    nf=4*n;freq=2*np.pi*np.fft.fftfreq(nf,dt);period=nf*dt
    radius=(int(np.floor(args.probe/args.h))+.5)*args.h
    area=4*np.pi*radius**2/G2
    fac=4*np.pi/(G2*period)
    edges=np.linspace(1.,args.emax,args.bins+1)
    mid=(edges[:-1]+edges[1:])/2
    arrays={s:np.zeros((len(center),args.bins)) for s in ('plus','minus','neutral')}
    # Outward and inward split is only a diagnostic at k>0.1; near mass shell
    # threshold the decomposition is ill-conditioned, so no clipping is used.
    split={s:{'out':0.,'in':0.,'net_from_components':0.,'net_direct':0.} for s in arrays}
    signed_outside={s:{'number':0.,'energy_mX':0.} for s in arrays}
    all_signed={'charge':0.,'charged_energy_mX':0.,'neutral_energy_mX':0.}
    norms=[]
    masks={'plus':freq<-1,'minus':freq>1,'neutral':abs(freq)>MH}
    for j,c in enumerate(center):
        idx=c+off;ok=(idx>=0)&(idx<N);window=np.zeros(n);window[ok]=win[ok]
        norms.append(float(np.sum(window**2)*dt))
        def field(name):
            x=np.zeros(n);x[ok]=data[name][idx[ok]]*window[ok]
            return dt*np.fft.fft(x,n=nf)
        Z=field('x')+1j*field('y');DZ=field('dx')+1j*field('dy')
        PZ=field('px')+1j*field('py');H=field('s_minus1');DH=field('ds');PH=field('ps')
        q=area/period*np.imag(np.conj(Z)*DZ)
        ec=-area/period*np.real(np.conj(PZ)*DZ)
        eh=-area/period*np.real(np.conj(PH)*DH)
        all_signed['charge']+=float(q.sum());all_signed['charged_energy_mX']+=float(ec.sum());all_signed['neutral_energy_mX']+=float(eh.sum())
        for label in ('plus','minus','neutral'):
            m=masks[label];ee=abs(freq[m]);Ns=(q[m] if label=='plus' else -q[m] if label=='minus' else eh[m]/ee)
            arrays[label][j]=np.histogram(ee,bins=edges,weights=Ns)[0]
            outside=(ee<edges[0])|(ee>edges[-1]);signed_outside[label]['number']+=float(Ns[outside].sum());signed_outside[label]['energy_mX']+=float(Ns[outside]@ee[outside])
            field_m=H[m] if label=='neutral' else Z[m];dfield_m=DH[m] if label=='neutral' else DZ[m]
            mass=MH if label=='neutral' else 1.;k=np.sqrt(np.maximum(ee**2-mass**2,0));usable=k>.1
            direction=-np.sign(freq[m]);chi=radius*field_m;dchi=field_m+radius*dfield_m
            aa=.5*(chi[usable]+dchi[usable]/(1j*direction[usable]*k[usable]));bb=chi[usable]-aa
            no=fac*k[usable]*abs(aa)**2;ni=fac*k[usable]*abs(bb)**2
            split[label]['out']+=float(no.sum());split[label]['in']+=float(ni.sum())
            split[label]['net_from_components']+=float((no-ni).sum())
            if label=='neutral':
                # The flux diagnostic uses physical temporal derivatives; a
                # finite time window can spoil its equality to omega*N.
                split[label]['net_direct']+=float(Ns[usable].sum())
            else:split[label]['net_direct']+=float(Ns[usable].sum())
    # Conservation of integrated current before frequency masking follows
    # solely from window partition and Parseval; it is not source positivity.
    Z=data['x']+1j*data['y'];DZ=data['dx']+1j*data['dy'];PZ=data['px']+1j*data['py']
    raw={'charge':float(dt*area*np.imag(np.conj(Z)*DZ).sum()),
         'charged_energy_mX':float(-dt*area*np.real(np.conj(PZ)*DZ).sum()),
         'neutral_energy_mX':float(-dt*area*(data['ps']*data['ds']).sum())}
    res={};positive={}
    for s,a in arrays.items():
        positive[s]=np.maximum(a,0)
        res[s]={'signed_number':float(a.sum()),'signed_energy_mX':float((a*mid[None,:]).sum()),
                'negative_bin_weight':float(np.maximum(-a,0).sum()),
                'positive_projection_number':float(positive[s].sum()),
                'outside_support':signed_outside[s], 'split_k_gt_0.1':split[s]}
    arrays_out={'energy_edges':edges,'energy_centers':mid,'window_center_mX_t':t[center],
                'window_norm_mX_t':np.array(norms),'window_width_mX_t':np.array(args.width),
                'probe_radius_mX_inverse':np.array(radius)}
    for s in arrays:
        arrays_out[s+'_signed']=arrays[s];arrays_out[s+'_positive']=positive[s]
    np.savez_compressed(ROOT/(args.tag+'.npz'),**arrays_out)
    result={'configuration':vars(args),'trace_sha256':hashlib.sha256(tracepath.read_bytes()).hexdigest(),
            'samples':N,'sample_dt_mX_inverse':dt,'probe_actual_mX_inverse':radius,
            'window_partition_max_error':covererr,'raw_flux':raw,'partitioned_all_frequency_flux':all_signed,
            'parseval_error':{k:all_signed[k]-raw[k] for k in raw},'species':res,
            'scope':'Time-frequency windowed radial field flux; positive projection is a model sensitivity, not a certified time-local particle flux. Off-shell, finite-window interference and missing earlier/later sources are retained as limitations.'}
    (ROOT/(args.tag+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
