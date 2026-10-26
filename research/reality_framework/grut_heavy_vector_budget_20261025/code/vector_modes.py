"""Finite-band primordial longitudinal modes on the supplied Higgs norm.

Reuses the author's own b448f25 canonical Proca method. No GMR abundance
normalization, reviewer mode solver, or full scalar/Goldstone evolution.
"""
import json, math
from pathlib import Path
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from background import ROOT, WORK, HH, write

_profiles={}


def load_profile(name='actual'):
    if name not in _profiles:
        file=ROOT/'HIGGS_PROFILE.npz'
        if name=='long':
            raw=np.load(WORK/'long.npz')['samples']
            aa,ff=raw[:,0],raw[:,14]
        else:
            data=np.load(file);aa,ff=data['a'],data['F_over_H']
        # Adjacent output targets can differ by less than one log(a) ulp
        # at the endpoint. They are the same sample for interpolation.
        logarithms=np.log(aa)
        assert np.all(np.diff(logarithms)>=0)
        unique,index=np.unique(logarithms,return_index=True)
        aa,ff=aa[index],ff[index]
        _profiles[name]=(aa,ff,PchipInterpolator(unique,ff*ff),float((ff[-1]**2-(6*HH)**2)*aa[-1]**2))
    return _profiles[name]


def mode(task):
    e,profile,x,end_factor,rtol,start_x=task
    aa,ff,interpolator,tail=load_profile('long' if profile=='long' else 'actual')
    mv=6*e*HH;k=x*math.sqrt(mv)
    mi=e*float(ff[0]) if profile!='constant' else mv
    a0=k/start_x
    phase=np.exp(1j*start_x)
    initial_mass2=(a0*mi)**2
    z2=initial_mass2/(k*k+initial_mass2);zz=math.sqrt(z2)
    omega=math.sqrt(k*k+initial_mass2)
    u=phase/math.sqrt(2*omega)
    zp_over_z=a0*k*k/(k*k+initial_mass2)
    initial=np.array([u/zz,zz*(-1j*omega-zp_over_z)*u])

    def mass2(a,infl):
        if infl:return mi*mi
        if profile=='constant':return mv*mv
        if profile=='smooth':return mv*mv+(mi*mi-mv*mv)/(a*a)
        if a<=aa[-1]:return e*e*float(interpolator(math.log(a)))
        return mv*mv+e*e*tail/(a*a)

    def rhs(l,y,infl):
        a=math.exp(l);am2=a*a*mass2(a,infl)
        inertia=am2/(k*k+am2)
        conversion=1/a if infl else a
        return np.array([conversion*y[1]/inertia,-conversion*am2*y[0]])

    atol=np.array([abs(initial[0])*1e-13,max(abs(initial[1])*1e-13,1e-30)])
    inf=solve_ivp(lambda l,y:rhs(l,y,True),(math.log(a0),0),initial,
                  method='DOP853',rtol=rtol,atol=atol)
    assert inf.success
    end=end_factor/math.sqrt(mv)
    # Explicitly stop at the end of the tabulated actual release before
    # taking the asymptotic tail; no solver steps across that junction.
    intervals=[0.,math.log(aa[-1]),math.log(end)] if profile in ['actual','long'] else [0.,math.log(end)]
    state=inf.y[:,-1];nfev=inf.nfev
    for left,right in zip(intervals[:-1],intervals[1:]):
        rad=solve_ivp(lambda l,y:rhs(l,y,False),(left,right),state,
                      method='DOP853',rtol=rtol,atol=atol)
        assert rad.success
        state=rad.y[:,-1];nfev+=rad.nfev
    B,pi=state;am2=end*end*mass2(end,False);inertia=am2/(k*k+am2)
    energy=.5*(abs(pi)**2/inertia+am2*abs(B)**2)
    occupation=energy/math.sqrt(k*k+am2)-.5
    wronskian=np.imag(np.conj(B)*pi)
    return dict(e=e,profile=profile,k_over_kstar=x,k_over_H=k,occupation=float(occupation),
                dN_dlnk=float(k**3*occupation/(2*math.pi**2)),end_factor=end_factor,
                rtol=rtol,start_x=start_x,nfev=nfev,
                absolute_Wronskian_error=float(abs(wronskian+.5)),
                scaled_Wronskian_error=float(abs(wronskian+.5)/max(abs(B)*abs(pi),1.)),
                m_inf_over_H=mi,m_vac_over_H=mv)


def integrate(rows):
    rows=sorted(rows,key=lambda r:r['k_over_H'])
    yy=np.array([r['dN_dlnk'] for r in rows]);xx=np.log([r['k_over_H'] for r in rows])
    return dict(N_A_band=float(np.trapezoid(yy,xx)),
                positive_only_integral=float(np.trapezoid(np.maximum(yy,0),xx)),
                negative_integral=float(np.trapezoid(np.minimum(yy,0),xx)),nodes=len(rows),
                k_min=rows[0]['k_over_H'],k_max=rows[-1]['k_over_H'])


def main():
    base=[(1e-7,'actual',float(x),20.,2e-10,1000.) for x in np.geomspace(.2,5.,97)]
    smooth=[(1e-7,'smooth',float(x),20.,2e-10,1000.) for x in np.geomspace(.2,5.,49)]
    changed_e=[(1e-6,'actual',float(x),20.,2e-10,1000.) for x in np.geomspace(.2,5.,49)]
    constant=[(1e-7,'constant',float(x),20.,2e-10,1000.) for x in np.geomspace(.2,5.,25)]
    controls=[]
    for x in [.3,1.,3.]:
        controls.extend([(1e-7,'actual',x,20.,2e-10,1000.),
                         (1e-7,'actual',x,30.,2e-10,1000.),
                         (1e-7,'actual',x,20.,2e-11,2000.),
                         (1e-7,'long',x,20.,2e-10,1000.)])
    tasks=base+smooth+changed_e+constant+controls
    with Pool(4) as pool:
        out=[]
        for idx,result in enumerate(pool.imap(mode,tasks),1):
            out.append(result)
            if idx%20==0:print('MODE_PROGRESS',idx,len(tasks),flush=True)
    n0=len(base);n1=n0+len(smooth);n2=n1+len(changed_e);n3=n2+len(constant)
    spectra=dict(actual_fine=integrate(out[:n0]),actual_coarse=integrate(out[:n0:2]),
                 smooth=integrate(out[n0:n1]),e_sensitivity=integrate(out[n1:n2]),
                 constant_vacuum=integrate(out[n2:n3]))
    result=dict(spectra=spectra,modes=out[:n3],controls=out[n3:],
                scope='One longitudinal Proca polarization; standard BD inflation and supplied mass profile. Finite positive band only; no coupled physical-phase/scalar production or thermal relic evolution.',
                normalization='d^3k/(2pi)^3 -> k^3 dlnk/(2pi^2), one longitudinal polarization; occupation=Hamiltonian/omega-1/2.')
    write(ROOT/'VECTOR_RESULTS.json',result)
    print('VECTOR_BAND',json.dumps(spectra),flush=True)


if __name__=='__main__':main()
