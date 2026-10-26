"""Disjoint decay ledger, including distributed and boosted vector decays.

All tree widths are independently derived in modes.vacuum. Particle
numbers are projections of the undamped prepared classical releases;
this is a postprocessed ledger, not a full interacting kinetic evolution.
"""
import json,math
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_laguerre
from scipy.optimize import brentq
from modes import vacuum
from run import ROOT,write

YQ=4.8294009076988175e-9;DM=4.4e-10


def vector_average(H,gamma,mass,q,a0):
    """A emitted with comoving momentum q at a0. Integrate its proper-time
    survival in RD; returned energy a*E(a)=sqrt(q^2+a^2*m^2), GeV.
    Mean over exp(-optical depth) uses Gauss-Laguerre in optical depth.
    """
    b=q/mass;scale=np.sqrt(H/gamma)
    def primitive(a):
        if b==0:return .5*a*a
        return .5*(a*np.sqrt(a*a+b*b)-b*b*np.arcsinh(a/b))
    p0=primitive(a0)
    def optical(a):return gamma/H*(primitive(a)-p0)
    # The nonrelativistic late limit is checked numerically for every source.
    xs,ws=roots_laguerre(64);leading=mass*np.sqrt(2*H/gamma)
    energy=leading*math.gamma(1.5);residual=0.
    for x,w in zip(xs,ws):
        upper=max(a0+scale*np.sqrt(2*x+1),b,scale)
        while optical(upper)<x:upper*=2
        aa=brentq(lambda a:optical(a)-x,a0,upper,xtol=1e-3,rtol=1e-12)
        energy+=w*(np.sqrt(q*q+aa*aa*mass*mass)-leading*np.sqrt(x))
        residual=max(residual,abs(optical(aa)-x)/max(x,1))
    # Exact energy-work balance for one particle: d(a E)/da=m^2 a/(a E).
    # Average decay energy equals birth conformal energy + expansion work.
    birth=np.sqrt(q*q+a0*a0*mass*mass)
    def work(l):
        a=np.exp(l);depth=optical(a)
        return 0. if depth>740 else np.exp(-depth)*mass*mass*a*a/np.sqrt(q*q+a*a*mass*mass)
    upper=max(a0+50*scale,20*b)
    ledger=birth+quad(work,np.log(a0),np.log(upper),epsabs=0,epsrel=1e-9,limit=300)[0]
    return dict(mean_conformal_energy_GeV=float(energy),birth_conformal_energy_GeV=float(birth),
                independent_expansion_work_energy_GeV=float(ledger),
                quadrature_energy_work_relative_error=float(abs(energy-ledger)/ledger),
                optical_inverse_max_relative_error=float(residual),q_GeV=float(q),a_birth=float(a0))


def main():
    modes=json.loads((ROOT/'MODE_RESULTS.json').read_text());p=vacuum();out=[]
    # At these couplings A is deeply nonrelativistic long before decay.
    # Test two very different emission epochs and both physical birth momenta.
    H=float(modes['ensemble_means']['10']['H_GeV']);checks=[]
    for a0 in [10.,1000.]:
        for mom in [p['pA_from_Sigma'],p['pA_from_YH']]:
            checks.append(vector_average(H,p['GammaAXX'],p['mA'],a0*mom,a0))
    for index in [10,17]:
        r=modes['ensemble_means'][str(index)];H=r['H_GeV'];nx=abs(r['NX_window']);ex=r['EX_window']
        nh,nl,ns=r['NYH'],r['NYL'],r['NSIG'];ny=nh+nl
        mean_a=lambda G:np.sqrt(np.pi*H/(2*G))
        rY=ny*p['mYL']/H*mean_a(p['GammaYL3'])/ex
        # Includes no unknown vector production. Heavy-cascade bracket is
        # [0, total Y number]; late NH is a central estimate, not an early
        # decay solution in the co-rotating mixing era.
        avec=p['mA']/H*mean_a(p['GammaAXX'])/ex
        rS=2*ns*avec;rH=nh*avec;rHmax=ny*avec
        pair_limit=.01*DM/(2*p['mX']*YQ)
        pairY=1.5*ny/nx;pairS=2*ns/nx;pairH=nh/nx
        rows=[]
        for x in modes['endpoint_rows']:
            if x['candidate']!=index:continue
            n=x['NYH']+x['NYL'];ee=x['EX_window'];xx=abs(x['NX_window'])
            ry=n*p['mYL']/H*mean_a(p['GammaYL3'])/ee
            unit=p['mA']/H*mean_a(p['GammaAXX'])/ee
            rows.append(dict(member=x['member'],R_Y=ry,R_Sigma=2*x['NSIG']*unit,
                             R_cascade=x['NYH']*unit,R_lower=ry+2*x['NSIG']*unit,
                             R_central=ry+(2*x['NSIG']+x['NYH'])*unit,
                             pairs_central=(1.5*n+2*x['NSIG']+x['NYH'])/xx))
        out.append(dict(candidate=index,H_GeV=H,cH=10,delta=r['delta'],epsilon_net=r['eps_net'],
            ymult_raw=r['YX_over_YQ_mean'],NX=nx,EX=ex,NYH=nh,NYL=nl,NSIG=ns,
            R_Y=rY,R_Sigma=rS,R_heavy_cascade_central=rH,R_heavy_cascade_max=rHmax,
            R_lower=rY+rS,R_central=rY+rS+rH,R_upper=rY+rS+rHmax,
            pairs_Y=pairY,pairs_Sigma=pairS,pairs_heavy_cascade=pairH,
            pairs_total_central=pairY+pairS+pairH,one_percent_pair_allowance=pair_limit,
            pairs_total_over_allowance=(pairY+pairS+pairH)/pair_limit,
            a_mean_decay=dict(YL=mean_a(p['GammaYL3']),Sigma=mean_a(p['GammaSigmaAA']),
                             YH_cascade=mean_a(p['GammaYHLA']),A=mean_a(p['GammaAXX'])),
            ensemble_rows=rows))
    write(ROOT/'DECAY_RESULTS.json',dict(vacuum=p,rows=out,vector_time_dilation_controls=checks,
        conserved_number_assumption='Late projected numbers conserved to decay. Every Y counted once after the cascade; vector energy is not added a second time to its source. No backreaction/thermalization history solved.',
        cascade_bracket='R_lower omits all heavy-cascade vectors; R_upper assigns at most one cascade vector to every Y quantum. The late NYH is a central diagnostic only.',
        radiation_deposit='All deposited radiation is compared in a^4 rho/H_i^4 units; Y uses the light daughter mass after cascade. Sterile vectors redshift with momentum before their own delayed decay.'))
    for r in out:print('DEBRIS',r['candidate'],'Y',r['R_Y'],'Sigma',r['R_Sigma'],'total range',r['R_lower'],r['R_central'],r['R_upper'],flush=True)


if __name__=='__main__':main()
