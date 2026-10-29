"""Reproduce one inherited history and recover all three outgoing spectra."""
from setup import *
import subprocess

def main():
    src = DYN/'code/radial.cpp'
    exe = WORK/('radial_'+sha(src)[:16])
    if not exe.exists():
        subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off',str(src),'-o',str(exe)],check=True)
    case=json.loads((DYN/'EXCITATION_RESULTS.json').read_text())['cases']['pulse10']
    par=case['params']; h,dt,L,ab,amp=[par[k] for k in ('h','dt','L','ab','amp')]
    r=(np.arange(round(L/h))+.5)*h;f,s=profile(r);shape=np.exp(-(r/12)**2)
    f*=1+amp*shape;s+=amp*(1-s)*shape
    vol=4*np.pi/3*((r+h/2)**3-(r-h/2)**3)
    py=-OMEGA*f;py*=2000/np.sum(-f*py*vol)
    initial=WORK/'initial.csv'
    np.savetxt(initial,np.c_[r,f,0*r,s,0*r,py,0*r],delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
    prefix=WORK/'pulse10'
    subprocess.run([str(exe),str(initial),str(prefix),str(h),str(dt),str(L),'1000','60',str(ab),'.2'],check=True,stdout=subprocess.DEVNULL)
    hist=np.genfromtxt(str(prefix)+'_history.csv',names=True,delimiter=',')
    with np.load(DYN/'EXCITATION_DATA.npz') as z: ref=z['pulse10_history']
    error=max(float(np.max(abs(hist[n]-ref[n])))/max(1,float(np.max(abs(ref[n])))) for n in hist.dtype.names)
    if error>1e-10: raise RuntimeError('Inherited history did not reproduce')
    d=np.genfromtxt(str(prefix)+'_trace.csv',names=True,delimiter=',');sample_dt=d['tau'][1]-d['tau'][0];d=d[:-1]
    nf=4*len(d);period=nf*sample_dt;omega=2*np.pi*np.fft.fftfreq(nf,sample_dt)
    transform=lambda x:sample_dt*np.fft.fft(x,n=nf)
    z=transform(d['x']+1j*d['y']);dz=transform(d['dx']+1j*d['dy'])
    area=4*np.pi*((int(np.floor(60/h))+.5)*h)**2/G2
    q=area/period*np.imag(np.conj(z)*dz)
    eh=-area/period*np.real(np.conj(transform(d['ps']))*transform(d['ds']))
    arrays={};stats={}
    for label,mask,weight in [('plus',omega < -1,q),('minus',omega > 1,-q),('neutral',abs(omega)>MH,eh/np.maximum(abs(omega),1e-300))]:
        e=abs(omega[mask]);w=weight[mask];order=np.argsort(e);e=e[order];w=w[order]
        # Finite-window signed artifacts are diagnosed; positive transport counts
        # are rescaled to the measured signed total, not to a desired outcome.
        target=w.sum();positive=np.maximum(w,0);positive*=target/positive.sum()
        arrays[label+'_E']=e;arrays[label+'_N']=positive
        stats[label]={'number':float(target),'energy':float(e@positive),
                      'negative_weight_fraction':float(np.maximum(-w,0).sum()/max(target,1e-300))}
    np.savez_compressed(ROOT/'SOURCE.npz',**arrays)
    out={'history_scaled_error':error,'params':par,'spectra':stats,
         'initial_core_charge':QREF-stats['plus']['number']+stats['minus']['number'],
         'trace_sha256':sha(str(prefix)+'_trace.csv'),'history_sha256':sha(str(prefix)+'_history.csv'),
         'source_kind':'supplied pulse, not cosmological formation'}
    write_json('SOURCE_CHECKS.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__': main()
