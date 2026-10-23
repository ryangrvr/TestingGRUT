#!/usr/bin/env python3
"""Outgoing linear-field diagnostics and explicitly conditional decay folding.

Near-threshold coherent flux is retained as unresolved when the local
mass correction is too large. No missing energy is forced into particles.
"""
import json
from pathlib import Path
import numpy as np
from scipy.signal.windows import tukey
from run_core import load_trace

ROOT=Path(__file__).resolve().parents[1]
G,V=.1,.924089637375362
MX=G*V;MS=np.sqrt(.06)*V;MD=np.sqrt(.0101)*V


def field_spectrum(trace,run,field,alpha=.05):
    t=trace[:,0];dt=t[1]-t[0];n=len(t);T=n*dt
    window=tukey(n,alpha)
    # Recorded probe index is rounded to the nearest cell centre in C++.
    c=run['parameters'];r=(round(c['probe']/c['h']-.5)+.5)*c['h']
    if field=='X':
        value=trace[:,1]+1j*trace[:,2];der=trace[:,4]+1j*trace[:,5]
        momentum=trace[:,7]+1j*trace[:,8];mass=1.
        perturb=abs((1+trace[:,3])**2-1+.1*abs(value)**2)
    else:
        value=trace[:,3];der=trace[:,6];momentum=trace[:,9];mass=np.sqrt(6)
        perturb=abs(3*((1+trace[:,3])**2-1)+trace[:,1]**2+trace[:,2]**2)
    omega=-2*np.pi*np.fft.fftfreq(n,dt)
    F=dt*np.fft.fft(value*window);D=dt*np.fft.fft(der*window);P=dt*np.fft.fft(momentum*window)
    area=4*np.pi*r*r
    enet=-area*np.real(P.conj()*D)/T
    qnet=area*np.imag(F.conj()*D)/T
    k2=omega*omega-mass*mass
    valid=k2>1e-4
    k=np.sign(omega[valid])*np.sqrt(k2[valid])
    outgoing=.5*(F[valid]+(D[valid]+F[valid]/r)/(1j*k))
    incoming=.5*(F[valid]-(D[valid]+F[valid]/r)/(1j*k))
    eout=area*abs(omega[valid]*k)*abs(outgoing)**2/T
    ein=area*abs(omega[valid]*k)*abs(incoming)**2/T
    number=eout/abs(omega[valid])
    # This cut controls local effective mass corrections, not all possible
    # particle/interaction or finite-radius approximation errors.
    controlled=k2[valid]>10*np.max(perturb)
    energy_time=-area*dt*np.sum(window**2*np.real(momentum.conj()*der))
    q_time=area*dt*np.sum(window**2*np.imag(value.conj()*der))
    diagnostics=dict(field=field,window_alpha=alpha,probe_radius_gv_units=r,
        max_local_mass_squared_correction=float(max(perturb)),
        controlled_k_squared_lower=float(10*max(perturb)),
        total_windowed_net_energy=float(sum(enet)),time_domain_windowed_net_energy=float(energy_time),
        total_windowed_net_charge=float(sum(qnet)),time_domain_windowed_net_charge=float(q_time),
        net_propagating_energy_from_measured_momenta=float(sum(enet[valid])),
        net_energy_from_out_minus_in_linear_relation=float(sum(eout-ein)),
        outgoing_linear_energy=float(sum(eout)),incoming_linear_energy=float(sum(ein)),
        controlled_outgoing_linear_energy=float(sum(eout[controlled])),
        controlled_incoming_linear_energy=float(sum(ein[controlled])),
        controlled_fraction_of_outgoing_linear_energy=float(sum(eout[controlled])/sum(eout)),
        outgoing_positive_frequency_Nbar=float(sum(number[omega[valid]>0])),
        outgoing_negative_frequency_Nbar=float(sum(number[omega[valid]<0])),
        controlled_positive_frequency_Nbar=float(sum(number[controlled & (omega[valid]>0)])),
        controlled_negative_frequency_Nbar=float(sum(number[controlled & (omega[valid]<0)])),
        endpoint_field_amplitude=float(abs(value[-1])),
        first_field_amplitude=float(abs(value[0])))
    return diagnostics,dict(omega=omega[valid],eout=eout,ein=ein,Nbar=number,controlled=controlled)


def fold(sx,ss,point):
    """Only the declared controlled/outgoing band; physical units."""
    old=ROOT.parent/'grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json'
    pop=json.loads(old.read_text())['rows'][1]['joint_screen_population']
    bx,bd,bv=[pop[k] for k in ('hidden_X_branch','hidden_D_branch','visible_branch')]
    edges=np.r_[0.,np.geomspace(1e-5,1000,241)]
    nx=np.zeros(len(edges)-1);nanti=nx.copy();nd=nx.copy();energy_in=0.;energy_out=0.
    mask=sx['controlled'];freq=sx['omega'][mask]
    number=sx['Nbar'][mask]/G**2
    momentum=np.sqrt(freq*freq-1)*MX
    nx+=np.histogram(momentum[freq>0]/MX,edges,weights=number[freq>0])[0]
    nanti+=np.histogram(momentum[freq<0]/MX,edges,weights=number[freq<0])[0]
    direct_E=float(sum(sx['eout'][mask])*V/G)
    direct_Q=float(sum(number*np.sign(freq)))
    # Real-s positive and negative FFT frequencies together comprise its
    # one physical real-field energy; no extra factor of two is inserted.
    mask=ss['controlled'];Es=abs(ss['omega'][mask])*MX
    ns=ss['Nbar'][mask]/G**2
    c,w=np.polynomial.legendre.leggauss(32);w=w/2
    for mass,branch,label in [(MX,bx,'X'),(MD,bd,'D')]:
        if not np.any(ns):continue
        gam=Es/MS;beta=np.sqrt(np.maximum(1-1/gam**2,0))
        pstar=np.sqrt(MS*MS/4-mass*mass)
        # Equal marginals for both daughters by isotropy. For X there is
        # one of each charge; D has two identical real particles.
        ed=gam[:,None]*(MS/2+beta[:,None]*pstar*c[None,:])
        pd=np.sqrt(np.maximum(ed*ed-mass*mass,0))
        weight=branch*ns[:,None]*w[None,:]
        hist=np.histogram(pd.ravel()/MX,edges,weights=weight.ravel())[0]
        if label=='X':nx+=hist;nanti+=hist
        else:nd+=2*hist
        energy_out+=float(2*np.sum(weight*ed))
    scalar_energy=float(np.dot(ns,Es));visible_energy=bv*scalar_energy
    energy_out+=visible_energy
    energy_in=scalar_energy
    result=dict(scope='Controlled outgoing linear band only, folded through free vacuum decays; unresolved coherent/soft flux is omitted and disclosed, not assigned particles.',
        N_X=float(sum(nx)),N_Xstar=float(sum(nanti)),N_D=float(sum(nd)),
        direct_free_X_net_charge=direct_Q,returned_net_charge=float(sum(nx)-sum(nanti)),
        direct_charged_field_energy_GeV=direct_E,scalar_energy_before_decay_GeV=scalar_energy,
        visible_decay_energy_GeV=visible_energy,
        scalar_to_daughters_energy_relative_residual=energy_out/max(energy_in,1e-100)-1,
        stable_rest_energy_GeV=float(MX*(sum(nx)+sum(nanti))+MD*sum(nd)),
        emitted_energy_at_birth_GeV=direct_E+scalar_energy-visible_energy,
        X_pair_branch=bx,D_pair_branch=bd,visible_branch=bv)
    result['kinetic_energy_at_birth_GeV']=result['emitted_energy_at_birth_GeV']-result['stable_rest_energy_GeV']
    table=np.column_stack([edges[:-1],edges[1:],nx,nanti,nd])
    return result,table


def main():
    runs=json.loads((ROOT/'CORE_RESULTS.json').read_text())['runs']
    point=json.loads((ROOT/'PARAMETER_POINT.json').read_text())
    results=[]
    for run in runs:
        trace=load_trace(run)
        dx,sx=field_spectrum(trace,run,'X');ds,ss=field_spectrum(trace,run,'s')
        folded,table=fold(sx,ss,point)
        name=run['name']
        np.savetxt(ROOT/'data'/f'{name}_controlled_particle_spectrum.csv',table,delimiter=',',
                   header='p_over_mX_low,p_over_mX_high,N_X,N_Xstar,N_D',comments='')
        windows=[]
        for a in [.025,.1]:
            ex,_=field_spectrum(trace,run,'X',a);es,_=field_spectrum(trace,run,'s',a)
            windows.append(dict(alpha=a,X_net=ex['total_windowed_net_energy'],
                                scalar_net=es['total_windowed_net_energy']))
        results.append(dict(run=name,X=dx,scalar=ds,conditional_particle_fold=folded,
                            window_controls=windows))
        print(name,'X/s Eflux',dx['total_windowed_net_energy'],ds['total_windowed_net_energy'],
              'controlled X/s fractions',dx['controlled_fraction_of_outgoing_linear_energy'],
              ds['controlled_fraction_of_outgoing_linear_energy'],flush=True)
    (ROOT/'SPECTRAL_RESULTS.json').write_text(json.dumps(dict(
        scope='Finite-time/radius outgoing-field diagnostics and partial free-decay particle spectra, not a complete asymptotic population.',
        controlled_band_definition='k^2 > 10 times the largest local effective-mass-squared correction, plus k^2>1e-4.',
        rows=results),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
