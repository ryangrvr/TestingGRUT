"""Conservative high-energy contact-cascade control in an assumed soft gas.

Massless unbroken-phase MB bath, central CM scattering angles only,
soft-bath number/energy feedback. This tests hard-source relaxation, not
condensate breakup, the thermal phase transition or a full collision kernel.
"""
import json,math,sys
import numpy as np
from scipy.special import gammainc,gamma
import ks_common as K
from moment_controls import save

Cxx=(4*K.LX)**2/(16*np.pi);Cxs=(2*K.GX**2)**2/(16*np.pi);Css=(6*K.LP)**2/(16*np.pi)

def soft_conversion(N,U,a,dt):
    total=N.sum();Q=N[0]-N[1];T=U/(3*total)
    root=np.sqrt(4*total*total-3*Q*Q);plus=(-total+root)/3;minus=(-total-root)/3
    ratio=(N[2]-plus)/(N[2]-minus)
    k=(Cxs/2)*a/(8*T*T)
    ratio*=np.exp(-1.5*k*(plus-minus)*dt)
    ns=(plus-ratio*minus)/(1-ratio)
    return np.array([(total-ns+Q)/2,(total-ns-Q)/2,ns])

def source_step(src,a,b):
    par=src['parent'];N=par['N_A0_over_H3'];g=par['gamma'];c=par['p_star_GeV']/K.H_I
    x,y=g*a*a/2,g*b*b/2
    dn=N*(np.exp(-x)-np.exp(-y))
    energy=N*c*np.sqrt(2/g)*gamma(1.5)*(gammainc(1.5,y)-gammainc(1.5,x))
    return dn,energy

def deposit(hard,soft,energy,species,number,q,grid,cut):
    if number==0:return energy
    if q<=cut:
        soft[species]+=number;return energy+number*q
    j=np.searchsorted(grid,q)
    if j==0 or j>=len(grid):raise ValueError('Momentum grid does not cover deposition')
    f=(q-grid[j-1])/(grid[j]-grid[j-1]);hard[species,j-1]+=number*(1-f);hard[species,j]+=number*f
    return energy

def run(key,n=161,dx=.04,phi_fraction=0.):
    src=K.SRC['vector_sources'][key];NQ=K.N_X
    mediator_key='e1e-6' if key.startswith('e1e-6') else 'e1e-7'
    mediator_pairs=K.SRC['mediator_debris']['totals'][mediator_key]['pairs_per_net_charge_total']
    total=6.171997264164328*NQ+2*mediator_pairs*NQ
    soft=np.array([(total*(1-phi_fraction)+NQ)/2,(total*(1-phi_fraction)-NQ)/2,total*phi_fraction])
    med=K.R_MED[mediator_key]
    U0=K.E_X*(1+med);U=U0
    predicted_Tref=(U0+src['per_sign_totals']['conformal_energy_both_signs'])/(3*total)*K.H_I/K.A_REF
    Tc=K.RATES.T_c*K.MX if hasattr(K,'RATES') else K.Rates().T_c*K.MX
    end=min(K.A_REF/4,K.A_REF*predicted_Tref/(2*Tc))
    start=src['parent']['a_unit']*.001
    grid=np.geomspace(.03,5e6,n);hard=np.zeros((3,n));a=start
    expectedN=total;expectedE=U0;maxE=maxN=maxQ=0.;mins=1.;steps=0;traj=[]
    channels={}
    for i in range(3):
        for j in range(3):
            ch=[]
            if i<2 and j<2:
                ch.append((Cxx/(2 if i==j else 1),(i,j)))
                if i!=j:ch.append((Cxs/2,(2,2)))
            elif i==2 and j==2:ch=[(Css/2,(2,2)),(Cxs,(0,1))]
            else:ch=[(Cxs,(i,j))]
            channels[i,j]=ch
    while a<end:
        b=min(a*np.exp(dx),end);step=np.log(b/a);mid=np.sqrt(a*b)
        N=soft.copy();H=hard.copy();En=U
        dn,de=source_step(src,a,b)
        if dn>0 and de>0:
            q=de/dn
            En=deposit(H,N,En,0,dn,q,grid,0.)
            En=deposit(H,N,En,1,dn,q,grid,0.)
        N=soft_conversion(N,En,mid,step/2)
        T=En/(3*N.sum());cut=10*T;recoil=2*T
        for j,q in enumerate(grid):
            if q<=cut:
                N+=H[:,j];En+=q*H[:,j].sum();H[:,j]=0
        incoming=H.copy();new=np.zeros_like(H);delta=np.zeros(3);dE=0.
        for j in range(n-1,-1,-1):
            q=grid[j]
            if q<=cut:
                # Interpolation of an outgoing energy just above the cut
                # can place weight in the adjacent lower soft bin.
                delta+=incoming[:,j];dE+=q*incoming[:,j].sum();continue
            for species in range(3):
                rate=[]
                for target in range(3):
                    for coeff,out in channels[species,target]:
                        # Exact leading hard/MB rate; retain only half of
                        # isotropic CM phase space for central scatterings.
                        G=.5*coeff*mid*N[target]/(4*q*T)
                        rate.append((G,target,out))
                totalrate=sum(r[0] for r in rate)
                number=incoming[species,j]/(1+step*totalrate);new[species,j]=number
                for G,target,out in rate:
                    events=step*G*number;delta[target]-=events;dE-=events*recoil
                    for z in (.3,.4,.5,.6,.7):
                        ev=events/5
                        for child,qout in zip(out,[z*(q+recoil),(1-z)*(q+recoil)]):
                            if qout<=cut:delta[child]+=ev;dE+=ev*qout
                            else:
                                pos=np.searchsorted(grid,qout)
                                if pos>j:raise RuntimeError('Cascade is not triangular')
                                f=(qout-grid[pos-1])/(grid[pos]-grid[pos-1])
                                incoming[child,pos-1]+=ev*(1-f);incoming[child,pos]+=ev*f
        N+=delta;En+=dE
        if N.min()<-1e-9*total or En<=0:
            raise RuntimeError('Bath positivity failed; decrease step')
        N=soft_conversion(N,En,mid,step/2)
        soft,hard,U=N,new,En;a=b;steps+=1
        expectedN+=2*dn;expectedE+=2*de
        actualE=U+float((hard*grid).sum());actualN=soft.sum()+hard.sum();actualQ=soft[0]-soft[1]+hard[0].sum()-hard[1].sum()
        maxE=max(maxE,abs(actualE/expectedE-1));maxN=max(maxN,abs(actualN/expectedN-1));maxQ=max(maxQ,abs(actualQ/NQ-1));mins=min(mins,float(soft.min()))
        if steps%5==0 or a==end:
            Tvis=U/(3*soft.sum())*K.H_I/a
            traj.append(dict(a=a,soft_T_GeV=Tvis,hard_energy_fraction=float((hard*grid).sum()/actualE),
                hard_number_fraction=float(hard.sum()/actualN),scalar_fraction=float(soft[2]/soft.sum())))
    Tconf=U/(3*soft.sum());kinetic=[]
    for species in (0,1):
        gg=sum(coeff*soft[target] for target in range(3) for coeff,out in channels[species,target])
        kinetic.append(.5*a*gg/(8*Tconf*Tconf))
    remain=src['parent']['N_A0_over_H3']*np.exp(-src['parent']['gamma']*a*a/2)
    return dict(source=key,n=n,dx=dx,initial_phi_fraction=phi_fraction,start_a=start,end_a=end,steps=steps,
        final=traj[-1],remaining_parent_fraction=remain/src['parent']['N_A0_over_H3'],soft_elastic_rate_over_H=kinetic,
        checks=dict(max_energy=maxE,max_number=maxN,max_charge=maxQ,min_soft_population=mins),trajectory=traj,
        scope='Assumed released soft gas with R2 number; massless contact/central-angle cascade plus soft 2->2 conversion. No initial coherent R1 evolution, number-changing relaxation, phase transition or object formation.')

if __name__=='__main__':
    rows=[]
    for key,n,dx in [('e1e-7_band',161,.04),('e1e-7_band',321,.02),('e1e-7_broad',161,.04),('e1e-6_band',161,.04)]:
        r=run(key,n,dx);rows.append(r);save('RELAXATION_CONTROLS.json',dict(rows=rows));print('RELAX',key,n,r['final'],r['checks'],flush=True)
