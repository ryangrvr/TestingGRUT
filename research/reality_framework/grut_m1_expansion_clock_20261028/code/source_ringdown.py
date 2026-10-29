"""Extend the supplied radial pulse; measure late emission without imposing decay."""
from support import *
from scipy import sparse
from scipy.sparse.linalg import spsolve
import subprocess, argparse, time

def geometry(h,L):
    r=(np.arange(round(L/h))+.5)*h
    vol=4*np.pi/3*((r+h/2)**3-(r-h/2)**3)
    face=4*np.pi*(r+h/2)**2/h;face[-1]*=2
    left=np.r_[0,face[:-1]]
    lap=sparse.diags([face[:-1]/vol[:-1],-(face+left)/vol,left[1:]/vol[1:]],[1,0,-1],format='csr')
    return r,vol,face,lap

def stationary_lattice(q,h,L):
    """Stationary FV equations matching radial.cpp, including its outer face."""
    r,v,face,lap=geometry(h,L);n=len(r);f,s=fixed.profile(r);om=fixed.OMEGA
    def residual(f,s,om):
        rs=-lap@s+(3*(s*s-1)+f*f)*s-.6*(s-1)**2
        rs[-1]-=face[-1]/v[-1]  # Dirichlet s=1 at the outer face
        return np.r_[-lap@f+(s*s+.1*f*f-om*om)*f,rs,
                     (om*(v@f**2)/(.01*q)-1)]
    for step in range(20):
        res=residual(f,s,om)
        if max(abs(res))<2e-10:break
        ff=-lap+sparse.diags(s*s+.3*f*f-om*om)
        ss=-lap+sparse.diags(9*s*s-3+f*f-1.2*(s-1))
        fs=sparse.diags(2*f*s)
        jac=sparse.bmat([[ff,fs,sparse.csr_matrix((-2*om*f)[:,None])],
          [fs,ss,sparse.csr_matrix((n,1))],
          [sparse.csr_matrix((2*om*v*f/(.01*q))[None,:]),sparse.csr_matrix((1,n)),
           sparse.csr_matrix([[v@f**2/(.01*q)]])]],format='csr')
        delta=spsolve(jac,-res);alpha=1.
        for _ in range(15):
            nf=f+alpha*delta[:n];ns=s+alpha*delta[n:2*n];no=om+alpha*delta[-1]
            if np.linalg.norm(residual(nf,ns,no))<np.linalg.norm(res):break
            alpha*=.5
        else:raise RuntimeError('Stationary lattice Newton failed')
        f,s,om=nf,ns,no
    error=float(max(abs(residual(f,s,om))));assert error<2e-9,error
    grad=.5*face*(np.diff(np.r_[f,0.])**2+np.diff(np.r_[s,1.])**2)
    pot=.75*(s*s-1)**2-.2*(s-1)**3+.5*s*s*f*f+.025*f**4
    cell=(v*(.5*om*om*f*f+pot)+grad)/.01
    return {'Q':q,'Omega':om,'energy_mX':float(cell.sum()),
            'energy_inside20_mX':float(cell[r<20].sum()),'max_stationary_equation_residual':error}

def spectrum(trace,begin,end,h):
    data=trace[(trace['tau']>=begin)&(trace['tau']<end)]
    dt=float(data['tau'][1]-data['tau'][0]);nf=4*len(data);period=nf*dt
    om=2*np.pi*np.fft.fftfreq(nf,dt)
    ft=lambda x:dt*np.fft.fft(x,n=nf)
    z=ft(data['x']+1j*data['y']);dz=ft(data['dx']+1j*data['dy'])
    pz=ft(data['px']+1j*data['py'])
    area=4*np.pi*((int(np.floor(60/h))+.5)*h)**2/.01
    charge=area/period*np.imag(np.conj(z)*dz)
    eflux=-area/period*np.real(np.conj(pz)*dz)
    hflux=-area/period*np.real(np.conj(ft(data['ps']))*ft(data['ds']))
    out={};arrays={}
    for label,mask,number,energy in [('plus',om<-1,charge,eflux),
            ('minus',om>1,-charge,eflux),('neutral',abs(om)>fixed.MH,hflux/np.maximum(abs(om),1e-300),hflux)]:
        e=abs(om[mask]);w=number[mask];en=energy[mask]
        order=np.argsort(e)
        out[label]={'signed_number':float(w.sum()),'negative_number_weight':float(np.maximum(-w,0).sum()),
           'frequency_times_number_energy_mX':float(e@w),'direct_energy_flux_mX':float(en.sum())}
        # Retain signed bins; these are diagnostics, not positivity-certified kinetic sources.
        arrays[label+'_E']=e[order];arrays[label+'_N_signed']=w[order];arrays[label+'_Eflux']=en[order]
    out['all_frequency_direct_energy_flux_mX']=float(eflux.sum()+hflux.sum())
    out['all_frequency_charge_flux']=float(charge.sum())
    return out,arrays

def compact_spectrum(arrays):
    """Diagnostic bins only; preserve signed number/energy and negative weights.

    Fine Fourier arrays remain reproducible from the retained code and local
    work traces. These compact arrays are not a certified kinetic source.
    """
    edges=np.linspace(0,4*np.pi,1001);out={'edges':edges}
    for label in ('plus','minus','neutral'):
        e=arrays[label+'_E'];n=arrays[label+'_N_signed'];ef=arrays[label+'_Eflux']
        out[label+'_N_signed']=np.histogram(e,bins=edges,weights=n)[0]
        out[label+'_N_negative_part']=np.histogram(e,bins=edges,weights=np.maximum(-n,0))[0]
        out[label+'_EN_signed']=np.histogram(e,bins=edges,weights=e*n)[0]
        out[label+'_Eflux']=np.histogram(e,bins=edges,weights=ef)[0]
    return out

def run(h,dt,tag):
    started=time.monotonic();L=200.;end=10000.;src=fixed.DYN/'code/radial.cpp'
    exe=WORK/('radial_'+sha(src)[:16])
    if not exe.exists():subprocess.run(['g++','-O3','-std=c++17','-ffp-contract=off',str(src),'-o',str(exe)],check=True)
    r,vol,face,lap=geometry(h,L);f,s=fixed.profile(r);shape=np.exp(-(r/12)**2)
    f*=1+.1*shape;s+=.1*(1-s)*shape
    py=-fixed.OMEGA*f;py*=2000/np.sum(-f*py*vol)
    initial=WORK/(tag+'_initial.csv')
    np.savetxt(initial,np.c_[r,f,0*r,s,0*r,py,0*r],delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
    prefix=WORK/tag
    subprocess.run([str(exe),str(initial),str(prefix),str(h),str(dt),str(L),str(end),'60','100','.2'],check=True)
    hist=np.genfromtxt(str(prefix)+'_history.csv',delimiter=',',names=True)
    trace=np.genfromtxt(str(prefix)+'_trace.csv',delimiter=',',names=True)
    startrow=hist[np.argmin(abs(hist['tau']-1000))];last=hist[-1]
    original_error=None
    if h==.2 and dt==.025:
        with np.load(fixed.DYN/'EXCITATION_DATA.npz') as z:ref=z['pulse10_history']
        original_error=max(float(np.max(abs(hist[k][:len(ref)]-ref[k])))/max(1,float(np.max(abs(ref[k])))) for k in hist.dtype.names)
        assert original_error<1e-10,original_error
    checkpoints=[]
    for t in (1000,3000,10000):
        row=hist[np.argmin(abs(hist['tau']-t))]
        stationary=stationary_lattice(float(row['Q20']/.01),h,L)
        checkpoints.append({'mX_t':t,'core20_charge':float(row['Q20']/.01),
          'core20_energy_mX':float(row['E20']/.01),
          'lattice_stationary_at_core20_charge':stationary,
          'core20_excess_over_same_grid_stationary_mX':float(row['E20']/.01-stationary['energy_inside20_mX'])})
    spectra={}
    for begin,stop in [(0,1000),(1000,10000),(3000,10000)]:
        key=f'{begin}_{stop}';sp,arrays=spectrum(trace,begin,stop,h);spectra[key]=sp
        np.savez_compressed(ROOT/(tag+'_spectrum_'+key+'.npz'),**compact_spectrum(arrays))
    np.savez_compressed(ROOT/(tag+'_history.npz'),history=hist)
    energy=hist['E']+hist['E_abs'];charge=hist['Q']+hist['Q_abs']
    out={'configuration':{'h':h,'dt':dt,'L':L,'tmax_mX':end,'absorber_start':100,'probe':60,'amplitude':.1},
      'source_code_sha256':sha(src),'history_sha256':sha(str(prefix)+'_history.csv'),
      'original_history_scaled_error':original_error,
      'initial_pulse_energy_mX':float(hist['E'][0]/.01),
      'initial_stationary_lattice':stationary_lattice(200000,h,L),
      'late_1000_to_10000':{'probe_energy_mX':float((last['Eflux_probe']-startrow['Eflux_probe'])/.01),
        'probe_charge':float((last['Qflux_probe']-startrow['Qflux_probe'])/.01),
        'absorber_energy_mX':float((last['E_abs']-startrow['E_abs'])/.01),
        'absorber_charge':float((last['Q_abs']-startrow['Q_abs'])/.01)},
      'checks':{'max_relative_global_energy_drift':float(np.max(abs(energy/energy[0]-1))),
        'final_global_energy_error_mX':float((energy[-1]-energy[0])/.01),
        'late_global_energy_change_mX':float((energy[-1]-energy[np.argmin(abs(hist['tau']-1000))])/.01),
        'max_relative_global_charge_drift':float(np.max(abs(charge/charge[0]-1)))},
      'checkpoints':checkpoints,'signed_spectrum_diagnostics':spectra,'seconds':time.monotonic()-started,
      'interpretation':'Classical isolated radial evolution only; no quantum or portal decay, imposed heat injection, gas collisions or expansion.'}
    save(tag+'.json',out);print(json.dumps(out,default=lambda x:x.tolist() if hasattr(x,'tolist') else x),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--h',type=float,default=.2);p.add_argument('--dt',type=float,default=.025)
    p.add_argument('--tag',default='ringdown_h02');a=p.parse_args();run(a.h,a.dt,a.tag)
