"""Independent phase-space normalization, clock and stationary mass checks."""
from setup import *
from scipy.integrate import quad,simpson
from scipy.special import kv
import importlib.util
from cold_completion import Seed
import core

def main():
    # Direct Breit-Wigner angular integral, transformed to arctan coordinates
    # to resolve the actual narrow width, versus the independently derived NWA.
    bw=[];aa=MH*WIDTH
    for ei,ej in [(1.4,1.4),(1.2,2.0),(1.5,3.0),(2.,4.)]:
        pi,pj=np.sqrt(ei*ei-1),np.sqrt(ej*ej-1)
        slo,shi=2+2*(ei*ej-pi*pj),2+2*(ei*ej+pi*pj)
        lower,upper=np.arctan((slo-6)/aa),np.arctan((shi-6)/aa)
        integrand=lambda z:np.sqrt(max(0,1-4/(6+aa*np.tan(z))))
        full=GH**4/(128*np.pi*ei*ej*pi*pj*aa)*quad(integrand,lower,upper,epsabs=1e-11,epsrel=1e-10)[0]
        nwa=np.pi*GH**2/(8*ei*ej*pi*pj) if slo<6<shi else 0.
        bw.append({'Ei':ei,'Ej':ej,'finite_width_sigma_v':full,'narrow_width_sigma_v':nwa,'relative_difference':full/nwa-1})
    thermal=[]
    for T in (.15,.3,.8):
        # Direct two-particle phase-space integral after using the on-shell
        # delta function. Compare with n_h^eq <Gamma*m_h/E_h> independently.
        integral=quad(lambda w:np.exp(-w/T)*np.sqrt(w*w-6),MH,np.inf,epsabs=1e-20,epsrel=1e-10)[0]
        inv=GH**2*np.sqrt(1/3)/(32*np.pi**3)*integral
        decay=WIDTH*6*T/(2*np.pi**2)*kv(1,MH/T)
        thermal.append({'T':T,'inverse_decay_reaction_density':inv,'decay_reaction_density':decay,
                        'relative_error':inv/decay-1})
    spec=importlib.util.spec_from_file_location('m1_deformation_check',M1/'code/deformation.py')
    deformation=importlib.util.module_from_spec(spec);spec.loader.exec_module(deformation)
    names=['translated_full_96','translated_full_64','aligned_full_64_long']
    cases=[];initial_source=json.loads((ROOT/'SOURCE_CHECKS.json').read_text())
    qi=initial_source['initial_core_charge']
    targets=[('initial_stationary_core',qi,None)]
    for name in names:
        p=ROOT/(name+'.json')
        if p.exists():
            z=json.loads(p.read_text());row=z['snapshots'][-1];targets.append((name,row['Qcore'],row))
    solutions={}
    for name,q,row in targets:
        sol=deformation.solve(q,.002,Seed(),length=160,tol=2e-9);solutions[name]=sol
        r=np.linspace(0,160,24001);f,fp,s,sp,j=sol.sol(r);om=float(sol.p[0])
        potential=B/4*(s*s-1)**2-C*(s-1)**3+.5*s*s*f*f+A/4*f**4
        energy=4*np.pi/G2*simpson(r*r*(.5*(fp*fp+sp*sp)+.5*om*om*f*f+potential),x=r)
        result={'case':name,'Q':q,'Omega':om,'stationary_energy_mX':energy,
                'BVP_residual':float(sol.rms_residuals.max()),
                'tangent_energy_error_mX':EREF+OMEGA*(q-QREF)-energy}
        if row is not None:
            result['cold_mass_excess_fraction']=(energy+row['Nplus']+row['Nminus'])/EREF-1
            result['tangent_cold_mass_excess_fraction']=row['tangent_cold_mass_excess_fraction']
        cases.append(result)
    hist=np.genfromtxt(WORK/'pulse10_history.csv',names=True,delimiter=',')
    gas=sum(v['energy'] for v in initial_source['spectra'].values())
    pulse=hist['E'][0]/G2
    retained=pulse-cases[0]['stationary_energy_mX']-gas
    source_ledger={'pulse_initial_lattice_energy_mX':pulse,'reference_stationary_energy_mX':EREF,
        'outgoing_on_shell_spectral_energy_mX':gas,'recorded_probe_energy_mX':hist['Eflux_probe'][-1]/G2,
        'postpulse_stationary_core_energy_mX':cases[0]['stationary_energy_mX'],
        'unprocessed_core_and_matching_energy_mX':retained,
        'unprocessed_energy_over_reference':retained/EREF,
        'caution':'This remainder is retained, not erased. Its release into the gas is uncomputed; it includes residual core excitation and finite-grid/spectral matching error.'}
    point=json.loads((M1/'CANDIDATE_M1.json').read_text());T=.02;g=gs=10.75
    entropy=2*np.pi**2/45*gs*T**3;H=np.sqrt(np.pi**2*g/90)*T*T/2.435e18
    ncore=point['final_population_matching_hypothesis']['YQ_target']*entropy/QREF
    ratio=H*M*M/ncore
    clock={'T_SM_GeV':T,'gstar':g,'gstar_entropy':gs,'ncore_GeV3':ncore,'ncore_over_mX3':ncore/M**3,
        'H_GeV':H,'H_mX2_over_ncore':ratio,'H_times_duration_tau1000':ratio*1000,
        'H_times_duration_tau10000':ratio*10000,'h_width_over_mX':WIDTH,
        'assumed_eta':1,'assumed_Q':QREF,'formation_epoch_derived':False}
    kernel=np.load(ROOT/'CORE_KERNELS.npz');eh=np.sqrt(1+kernel['k']**2)+OMEGA;opened=eh>MH
    kh=np.sqrt(eh[opened]**2-6);w=kernel['W'][opened]
    hrate=np.pi*(w[:,0,2]+w[:,1,2])/(kh*eh[opened])
    clock['max_sampled_h_core_conversion_over_vacuum_decay'] = float(np.max((ncore/M**3)*hrate/(WIDTH*MH/eh[opened])))
    # Finite-Q sensitivity is measured, not substituted into the frozen-profile run.
    response=[]
    if 'translated_full_96' in solutions:
        oldprofile,oldomega=core.profile,core.OMEGA
        sol=solutions['translated_full_96'];newprofile=lambda r:(sol.sol(np.minimum(r,160))[0],sol.sol(np.minimum(r,160))[2])
        for k in (.4,1.,1.5):
            base,*_=core.strengths(k)
            core.profile=newprofile;core.OMEGA=float(sol.p[0]);new,*_=core.strengths(k)
            core.profile=oldprofile;core.OMEGA=oldomega
            response.append({'incident_antiparticle_k':k,'exchange_sigma_reference':np.pi/k**2*base[0,1],
                'exchange_sigma_at_final_charge':np.pi/k**2*new[0,1],
                'exchange_relative_change':new[0,1]/base[0,1]-1,
                'neutral_sigma_reference':np.pi/k**2*base[2,1],
                'neutral_sigma_at_final_charge':np.pi/k**2*new[2,1]})
        z=np.load(ROOT/'translated_full_96.npz');e=z['E'];Y=z['Y'][:-1];n=len(e)
        edges=np.r_[1,(e[1:]+e[:-1])/2,e[-1]+(e[-1]-e[-2])/2]
        phase=np.diff(np.maximum(edges*edges-1,0)**1.5)/(6*np.pi**2)
        clock['maximum_resolved_MB_occupation']=float(np.max(Y/np.tile(phase,2)[:,None])*ncore/M**3)
    write_json('PHYSICS_CHECKS.json',{'finite_width_comparison':bw,'thermal_inverse_decay_normalization':thermal,
        'stationary_mass_checks':cases,'source_energy_ledger':source_ledger,'clock':clock,
        'finite_Q_kernel_sensitivity':response})
    print(json.dumps({'mass':cases,'source_ledger':source_ledger,'clock':clock,'finite_Q':response},indent=2))

if __name__=='__main__':main()
