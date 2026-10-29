"""Conditional M1 free flight from the saved kinetic nodes, not a full transfer function.

Radiation: photons/e± in equilibrium with instantaneous neutrino decoupling at
2 MeV and conserved electromagnetic entropy; after that T_nu ∝ 1/a. Neff=3.044,
standard Planck-like H0 and matter density. Actual M1 gas collisions after the
checkpoint and its perturbation correlations are not included.
"""
from pathlib import Path
import json, math
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator
from scipy.special import expit
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
T0K=2.7255; KB_EV=8.617333262145e-5; T0=T0K*KB_EV/1e6 # MeV
TI=2.0; ME=.51099895; H0=67.4; h=H0/100; OM=.315; NEFF=3.044; C=299792.458
OG=2.4728e-5/h**2; ONU=OG*7/8*(4/11)**(4/3)*NEFF; OR=OG+ONU; OL=1-OM-OR

def e_dof(T):
    y=ME/T
    def w(x):return expit(-np.sqrt(x*x+y*y))
    ir=quad(lambda x:x*x*np.sqrt(x*x+y*y)*w(x),0,80,epsabs=1e-10,epsrel=2e-9)[0]
    ip=quad(lambda x:x**4/np.sqrt(x*x+y*y)*w(x),0,80,epsabs=1e-10,epsrel=2e-9)[0]
    rho=4/(2*np.pi**2)*ir;press=4/(6*np.pi**2)*ip
    return 30*rho/np.pi**2,45*(rho+press)/(2*np.pi**2)

def background(N=4500,a_end=None):
    ts=np.geomspace(.005,TI,500);v=np.array([e_dof(t) for t in ts]);rg=PchipInterpolator(np.log(ts),v[:,0]);sg=PchipInterpolator(np.log(ts),v[:,1]);
    ai=T0/TI*(2/(2+float(sg(np.log(TI)))))**(1/3)
    ae=OR/OM;aa=np.geomspace(ai,ae if a_end is None else a_end,N);ag=T0/ts*(2/(2+v[:,1]))**(1/3)
    invert=PchipInterpolator(np.log(ag[::-1]),np.log(ts[::-1]));Tc=np.where(aa<=ag[0],np.exp(invert(np.minimum(np.log(aa),np.log(ag[0])))),T0/aa)
    tc=np.maximum(Tc,.005);gre=np.where(Tc>.005,rg(np.log(tc)),0.)
    heat=aa*Tc/T0
    E=np.sqrt(OM/aa**3+(OG*heat**4*(1+gre/2)+ONU)/aa**4+OL)
    return aa,E,ai,ae,{'initial_T_MeV':TI,'T0_K':T0K,'H0_km_s_Mpc':H0,'Omega_m':OM,'Omega_gamma':OG,'Omega_radiation_today':OR,'N_eff':NEFF,
       'a_initial':ai,'a_equality':ae,'z_equality':1/ae-1,'electron_g_rho_initial':float(rg(np.log(TI))),'electron_g_s_initial':float(sg(np.log(TI))),
       'definition':'Standard background; instantaneous neutrino decoupling at 2 MeV, equilibrium e± and photons, entropy transfer, and constant Neff. Neglects residual neutrino heating details and interactions of M1 species.'}

def quantile(values,weights,q):
    order=np.argsort(values);v=values[order];w=weights[order];return np.interp(q,(np.cumsum(w)-.5*w)/w.sum(),v)

def run_one(bins,aa,E,ai):
    z=np.load(RF/'grut_m1_core_matching_20261029'/f'late_kin_h0025_{bins}.npz');energy=z['E'];n=len(energy);y=z['Y'][:,-1];raw=y[:n]+y[n:2*n];weights=np.maximum(raw,0);p=np.sqrt(energy*energy-1);v=p/energy
    source=json.loads((RF/'grut_m1_core_matching_20261029'/f'late_kin_h0025_{bins}.json').read_text());last=source['snapshots'][-1];mref=187479.9981052459
    cold=mref*(1+last['stationary_cold_excess_fraction']);mass_fraction=weights.sum()/(cold+source['held_unassigned_recomputed_mX'])
    # log-a quadrature, exact relativistic speed; path is an upper bound only
    # within the declared collisionless background after the checkpoint.
    Q=p[:,None]*ai/aa;speed=Q/np.sqrt(1+Q*Q)
    dist=(C/H0)*np.trapezoid(speed/(aa*E),np.log(aa),axis=1)
    N=weights.sum();moments={'number':float(N),'charge_fraction':float(N/200000),'mass_fraction_of_cold_plus_held':float(mass_fraction),
        'rms_velocity_c':float(np.sqrt(weights@(v*v)/N)),'mean_velocity_c':float(weights@v/N),'kinetic_energy_mX':float(weights@(energy-1)),
        'fraction_floor_bin':float(weights[0]/N),'negative_clipped_weight':float(np.maximum(-raw,0).sum()),
        'mean_free_streaming_Mpc':float(weights@dist/N),'rms_free_streaming_Mpc':float(np.sqrt(weights@(dist*dist)/N)),
        'p10_Mpc':float(quantile(dist,weights,.1)),'p50_Mpc':float(quantile(dist,weights,.5)),'p90_Mpc':float(quantile(dist,weights,.9)),
        'max_Mpc':float(dist[weights>1e-8].max()),'cold_total_mX':float(cold),'held_unassigned_mX':float(source['held_unassigned_recomputed_mX'])}
    return moments,{'E':energy,'weights':weights,'p_at_2MeV':p,'v_at_2MeV':v,'lambda_Mpc':dist}

def save_psd(bins,arr,ai,logsigma=.18):
    # Smooth the discrete count measure in log momentum.  CLASS needs f0(q),
    # where dN ∝ q^2 f0(q) dq, hence f0 ∝ (dN/d ln q)/q^3.
    pscale=.05;qnodes=arr['p_at_2MeV']/pscale;weights=arr['weights'];q=np.geomspace(.001,25,2400)
    dif=(np.log(q[:,None])-np.log(qnodes[None,:]))/logsigma
    dlog=np.exp(-dif*dif/2)@weights/(np.sqrt(2*np.pi)*logsigma)
    f=dlog/q**3
    # Small, smooth, decreasing tails ensure CLASS extrapolation stays finite.
    f+=1e-20*np.exp(-q/3)
    norm=np.trapezoid(q*q*f,q);f/=norm
    path=ROOT/f'psd_{bins}.dat';np.savetxt(path,np.c_[q,f],fmt='%.14e')
    return {'p_scale_in_mX':pscale,'log_smoothing_sigma':logsigma,'number_integral_q2f':float(np.trapezoid(q*q*f,q)),
            'q_mean_from_smoothed_psd':float(np.trapezoid(q**3*f,q)),'q_mean_from_nodes':float(weights@qnodes/weights.sum()),
            'q2_mean_from_smoothed_psd':float(np.trapezoid(q**4*f,q)),'q2_mean_from_nodes':float(weights@(qnodes*qnodes)/weights.sum()),
            'warning':'Smoothed construction is an additional phase-space ansatz; moments shift by lognormal-kernel factors. Zero or negative numerical node weights clipped; no angular or perturbation correlations are derived.'}

def main():
    aa,E,ai,ae,bg=background();out={'background':bg,'populations':{},'psd':{},'scope':'Per-particle collisionless trajectories from the conditional a/a0=10 M1 kinetic output to equality. Spectrum is a numerical homogeneous-population model, not a derived relic phase space. Paths are not matter-power transfer functions.'};arrays={'a':aa,'E_over_H0':E}
    for b in (64,96):
        s,arr=run_one(b,aa,E,ai);out['populations'][str(b)]=s;out['psd'][str(b)]=save_psd(b,arr,ai)
        for k,v in arr.items():arrays[f'{k}_{b}']=v
        if b==96:
            # Sensitivity to a delayed start; this does not model collisions.
            weights=arr['weights'];mom=arr['p_at_2MeV'];q=mom[:,None]*ai/aa;v=q/np.sqrt(1+q*q)
            out['delayed_ballistic_mean_Mpc']={}
            for factor in (1,2,4,10):
                sel=aa>=factor*ai;dist=(C/H0)*np.trapezoid(v[:,sel]/(aa[sel]*E[sel]),np.log(aa[sel]),axis=1)
                out['delayed_ballistic_mean_Mpc'][str(factor)]=float(weights@dist/weights.sum())
            # Smooth analytic PSD with exactly the discrete <q²>; it is an
            # alternate assumed distribution for CLASS convergence controls.
            q=np.geomspace(.0001,25,2000);qnodes=mom/.05;q2=float(weights@(qnodes*qnodes)/weights.sum());q0=np.sqrt(2*q2/3)
            f=np.exp(-(q/q0)**2)+1e-20*np.exp(-q/3);f/=np.trapezoid(q*q*f,q)
            np.savetxt(ROOT/'psd_maxwell96.dat',np.c_[q,f],fmt='%.14e')
            out['psd']['maxwell96']={'q0':float(q0),'target_node_mean_q2':q2,'mean_q2_from_psd':float(np.trapezoid(q**4*f,q)),
                'mean_q_from_psd':float(np.trapezoid(q**3*f,q)),'scope':'Moment-matched Maxwell-like ansatz, not the measured M1 distribution.'}
            a_today,E_today,_,_,_=background(6500,a_end=1.);late,late_arr=run_one(b,a_today,E_today,ai)
            out['collisionless_through_today']={k:v for k,v in late.items() if 'free_streaming' in k or k.endswith('_Mpc')}
            arrays['lambda_to_today_Mpc_96']=late_arr['lambda_Mpc']
            # Toy: instantaneous independent, isotropic displacements from a
            # common initial density perturbation; not a gravitating transfer.
            fmass=s['mass_fraction_of_cold_plus_held'];paths=arr['lambda_Mpc'];out['toy_ballistic_at_equality']={}
            for kh in (1.,10.,100.,1000.,10000.):
                k=kh*h;j0=np.sinc(k*paths/np.pi);amp=1-fmass+fmass*float(weights@j0/weights.sum())
                out['toy_ballistic_at_equality'][str(kh)]={'delta_power_fraction':amp*amp-1,'bound_on_abs_fraction_small_perturbation':float(fmass*k*k*(weights@(paths*paths)/weights.sum())/3)}
    # Resolution control in log(a) for the weighted mean.
    a2,E2,_,_,_=background(9000);control,_=run_one(96,a2,E2,ai);out['checks']={'mean_path_4500_vs_9000_fraction':out['populations']['96']['mean_free_streaming_Mpc']/control['mean_free_streaming_Mpc']-1}
    np.savez_compressed(ROOT/'free_stream.npz',**arrays);(ROOT/'free_stream.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
