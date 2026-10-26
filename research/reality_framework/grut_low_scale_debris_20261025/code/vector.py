"""Light longitudinal Proca source with explicit changing-mass approximation.

This tests the withdrawn constant-vacuum-mass substitution. It does not
replace the coupled release-induced scalar/vector perturbation problem.
"""
import numpy as np
from scipy.integrate import solve_ivp
from model import minimum,yield_factor,YQ
from modes import vacuum
from run import ROOT,write


def mode(k,mv,mi,profile,end_factor=30,rtol=2e-10,start_x=1000):
    a0=k/start_x
    def m2(a,infl):
        if profile=='constant':return mv*mv
        return mi*mi if infl else mv*mv+(mi*mi-mv*mv)/(a*a)
    mstart=np.sqrt(m2(a0,True));z2=a0*a0*mstart*mstart/(k*k+a0*a0*mstart*mstart)
    z=np.sqrt(z2);v=np.exp(1j*start_x)/np.sqrt(2*k)
    zprime_over_z=a0*k*k/(k*k+a0*a0*mstart*mstart)
    initial=np.array([v/z,z*(-1j*k-zprime_over_z)*v])
    def rhs(l,y,infl):
        a=np.exp(l);mass=m2(a,infl);zz=a*a*mass/(k*k+a*a*mass)
        dt=1/a if infl else a
        return np.array([dt*y[1]/zz,-dt*a*a*mass*y[0]])
    atol=np.array([abs(initial[0])*1e-13,max(abs(initial[1])*1e-13,1e-30)])
    inf=solve_ivp(lambda l,y:rhs(l,y,True),[np.log(a0),0],initial,method='DOP853',rtol=rtol,atol=atol)
    astar=1/np.sqrt(mv);end=end_factor*astar
    rad=solve_ivp(lambda l,y:rhs(l,y,False),[0,np.log(end)],inf.y[:,-1],method='DOP853',rtol=rtol,atol=atol)
    assert inf.success and rad.success
    B,P=rad.y[:,-1];mass=m2(end,False);zz=end*end*mass/(k*k+end*end*mass)
    energy=.5*(abs(P)**2/zz+end*end*mass*abs(B)**2)
    occupation=energy/np.sqrt(k*k+end*end*mass)-.5
    scaled_w=abs(np.imag(np.conj(B)*P)+.5)/max(abs(B)*abs(P),1.)
    return dict(occupation=float(occupation),scaled_Wronskian_error=float(scaled_w),
                nfev=inf.nfev+rad.nfev,end_factor=end_factor,rtol=rtol,start_x=start_x)


def main():
    import json
    d=json.loads((ROOT/'DECAY_RESULTS.json').read_text());p=vacuum();H=d['rows'][0]['H_GeV']
    mi=minimum(H,10.)['vector_mass_over_H'];mv=p['mA']/H;kstar=np.sqrt(mv)
    rows=[]
    for factor in [.1,.3,1.,3.]:
        a=mode(factor*kstar,mv,mi,'constant');b=mode(factor*kstar,mv,mi,'changing')
        rows.append(dict(k_over_kstar=factor,constant=a,changing=b,number_ratio=b['occupation']/a['occupation']))
    tight=mode(kstar,mv,mi,'changing',rtol=2e-11,start_x=2000)
    shorter=mode(kstar,mv,mi,'changing',end_factor=20)
    ratios=[r['number_ratio'] for r in rows];estimates=[]
    for r in d['rows']:
        # Stated GMR constant-mass scaling, then multiply by the calculated
        # transfer ratio. This abundance normalization remains approximate.
        rho_s=4.4e-10*np.sqrt(p['mA']/6e-15)*(H/1e14)**2
        N0=rho_s/(p['mA']*yield_factor(H))
        R_unit=p['mA']/H*np.sqrt(np.pi*H/(2*p['GammaAXX']))/r['EX']
        estimates.append(dict(candidate=r['candidate'],constant_mass_naive_pairs_per_net=N0/r['NX'],
            corrected_pairs_per_net_range=[N0*min(ratios)/r['NX'],N0*max(ratios)/r['NX']],
            corrected_R_range=[N0*min(ratios)*R_unit,N0*max(ratios)*R_unit]))
    write(ROOT/'VECTOR_RESULTS.json',dict(H_GeV=H,mvac_over_H=mv,minflation_over_H=mi,
        frozen_angle_mass_ratio_squared=(mv/mi)**2,modes=rows,estimates=estimates,
        tighter_start_and_tolerance=tight,shorter_late_endpoint=shorter,
        scope='Original Proca mode solver for m(a)^2=mvac^2+(minf^2-mvac^2)/a^2 after release and stationary inflation before it. BD initial state; approximate GMR abundance normalization. Not the full oscillating-Higgs or release-induced physical-phase source.',
        withdrawn='The 1060-pair figure in 2793dab used a constant vacuum mass through inflation and was withdrawn by f1aa7df. It is not used as evidence against either low-scale release.'))
    print('VECTOR_MASS_TRANSFER',[(x['k_over_kstar'],x['number_ratio']) for x in rows],flush=True)
    print('VECTOR_PAIRS',estimates,flush=True)


if __name__=='__main__':main()
