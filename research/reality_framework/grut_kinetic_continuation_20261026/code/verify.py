"""Same-author numerical checks; none certify the missing physical closures."""
import json
import numpy as np
from scipy.integrate import quad
from bootstrap import ROOT
import ks_common as K
from stage2_moments import N2,Case,Solver,RATES
from thermo_state import TH,species,initial_equilibrium
from moment_controls import save
from relax_control import source_step

def phase_bound():
    tc=RATES.T_c;rows=[]
    for R in [0.,K.R_VEC['e1e-7_band']+K.R_MED['e1e-7'],K.R_VEC['e1e-7_broad']+K.R_MED['e1e-7'],K.R_VEC['e1e-6_band']+K.R_MED['e1e-6']]:
        e=K.rho_onset(R)/(K.MX*K.N_Q0)
        ptemplate=.006571055435898394+.0028047657247771878
        for label,p in [('template_inventory',ptemplate),('template_plus_latent_s',ptemplate+1.2787693395),('orbit_count',2.5859986321)]:
            lower=max(0.,(e-1-K.RS*p)/(3*(1+2*p)))
            rows.append(dict(R_total=R,label=label,inventory_per_net=p,minimum_thermal_T_over_m=lower,
                exceeds_broken_phase_temperature=bool(lower>tc),T_c_over_m=tc))
    save('INITIAL_STATE_DOMAIN.json',dict(rows=rows,
        inequality='rho/(m nQ) <= 1 + (m_s/m) p + 3(T/m)(1+2p), for the supplied ideal Bose gases and objects with E/Q<=m.',
        scope='Necessary check on forcing a release bracket into the vacuum-mass thermal closure. A coherent release has no temperature assigned by this inequality. Latent scalar number is supplied, not derived.'))

def main():
    groups={};stoich=[]
    for ch,d in N2.items():
        change=np.array(d['rev'])-np.array(d['fwd'])
        stoich.append(dict(ch=ch,charge=int(change@[1,-1,0]),inventory=int(change@[0,1,1]),declared=d['dP']))
    helper=Solver(Case('check',41.707077,rule='b',catalysis=False))
    station=[]
    for t,u in [(7.,.8611951308217225),(1.,.2),(.1,.8)]:
        st=dict(t=t,up=u,um=-u,us=0.,n_c=0.,n_cb=0.,n_plus=1.,n_minus=1.,n_s=1.)
        station.append(max(abs(g) for g,dp,dn in helper.rates(st).values()))
    for t in [.1,1.,7.]:
        st=dict(t=t,up=K.OMEGA_REF,um=-K.OMEGA_REF,us=0.,n_c=1.,n_cb=0.,n_plus=1.,n_minus=1.,n_s=1.)
        station.append(abs(helper.rates(st)['capture'][0]))
    groups['network_inventory_and_equilibrium']=dict(passed=all(r['charge']==0 and r['inventory']==r['declared'] for r in stoich) and max(station)<1e-20,
        stoichiometry=stoich,stationarity_absolute=max(station),scope='MB-fugacity detailed balance, not the full nonzero-mu Bose collision operator.')

    errs=[];derrs=[]
    for t,u in [(.1,.7),(1.,.5),(7.,.86),(7.,1-1e-6)]:
        a=species(t,1.,u);b=TH.quad_species(t,1.,u)
        errs.append(max(abs(a[k]/b[k]-1) for k in ['n','rho','P']))
        h=1e-5
        ap=species(t*np.exp(h),1.,u);am=species(t*np.exp(-h),1.,u)
        for f,d in [('n','dn_dt'),('rho','drho_dt')]:
            der=(ap[f]-am[f])/(2*h*t);derrs.append(abs(der/a[d]-1))
    initerr=[];nq=K.N_Q0/K.MX**3
    for R in [0.,K.R_VEC['e1e-7_band']+K.R_MED['e1e-7']]:
        for cap,ecap in [(1-1e-6,1.),(K.OMEGA_REF,K.E_OVER_QM_REF)]:
            st=initial_equilibrium(nq,K.rho_onset(R)/K.MX**4,cap,ecap)
            q=species(st['t'],1.,st['up'])['n']-species(st['t'],1.,st['um'])['n']+st['n_cond']
            initerr.append(abs(q/nq-1))
    groups['Bose_state_and_initial_matching']=dict(passed=max(errs)<1e-8 and max(derrs)<1e-6 and max(initerr)<1e-10,
        quadrature_relative=max(errs),susceptibility_finite_difference=max(derrs),initial_charge=max(initerr))

    source_err=[]
    for key in ['e1e-7_band','e1e-7_broad','e1e-6_band']:
        src=K.SRC['vector_sources'][key];par=src['parent'];g=par['gamma'];N=par['N_A0_over_H3'];c=par['p_star_GeV']/K.H_I
        for lo,hi in [(.01,.3),(.3,2.),(2.,6.)]:
            # y=sqrt(g/2)*a; compare time-domain quadrature to closed moments.
            a=lo*np.sqrt(2/g);b=hi*np.sqrt(2/g)
            exact=source_step(src,a,b)
            numeric=(N*quad(lambda y:2*y*np.exp(-y*y),lo,hi,epsabs=0,epsrel=1e-12)[0],
                     N*c*np.sqrt(2/g)*quad(lambda y:2*y*y*np.exp(-y*y),lo,hi,epsabs=0,epsrel=1e-12)[0])
            source_err.extend(abs(x/y-1) for x,y in zip(exact,numeric))
    rel=json.loads((ROOT/'RELAXATION_CONTROLS.json').read_text())['rows']
    assert len(rel)==4, 'Run all four hard controls before verification'
    conservation=max(r['checks'][k] for r in rel for k in ['max_energy','max_number','max_charge'])
    coarse,fine=rel[:2];xx=np.log([r['a'] for r in coarse['trajectory']]);xf=np.log([r['a'] for r in fine['trajectory']])
    yy=np.array([r['hard_energy_fraction'] for r in coarse['trajectory']]);yf=np.array([r['hard_energy_fraction'] for r in fine['trajectory']])
    diff=max(abs(yy-np.interp(xx,xf,yf)))
    groups['conditional_hard_source_control']=dict(passed=max(source_err)<1e-10 and conservation<1e-10 and diff<.01,
        source_quadrature=max(source_err),conservation=conservation,refined_hard_energy_fraction_absolute=diff,
        scope='Conservative massless contact cascade in a supplied soft MB gas; no proof of coherent-release thermalization.')

    prod=json.loads((ROOT/'MOMENT_CONTROLS.json').read_text())['cases']
    follow=json.loads((ROOT/'FOLLOWUP_CONTROLS.json').read_text())['cases']
    assert len(prod)==13 and len(follow)==11
    assert all('error' not in r for r in [*prod.values(),*follow.values()])
    e=max(r['checks']['max_energy_work_error'] for r in prod.values());q=max(r['checks']['max_charge_error'] for r in prod.values())
    groups['production_moment_conservation']=dict(passed=e<1e-7 and q<1e-6,energy_work=e,charge=q,case_count=len(prod),
        scope='Conservation of the supplied moment recipes; not validation of formation or kinetic closure.')
    val=lambda r:r['final']['f_at_stop']
    refinement={
        'gas_Radau_tolerance':abs(val(follow['gas_band_coarse'])/val(prod['e1e7_band_gas'])-1),
        'gas_BDF_vs_Radau':abs(val(follow['gas_band_BDF'])/val(prod['e1e7_band_gas'])-1),
        'instant_Radau_tolerance':abs(val(follow['instant_band_coarse'])/val(prod['e1e7_band_instant'])-1)}
    groups['same_author_method_and_tolerance']=dict(passed=max(refinement.values())<1e-4,relative_changes=refinement,
        scope='Distinct numerical methods by one author; no independent verifier or critic ran.')
    B=K.Background();err=[]
    for x in [16.,17.,18.,19.]:err.append(abs(B.H_of_x(x+.5)/B.H_of_x(x)-np.exp(-1)))
    groups['late_radiation_background']=dict(passed=max(err)<1e-12,absolute_ratio_error=max(err),
        Tvis_at_x18_GeV=B.T_of_x(18.),scope='Supplied radiation background; no silent 10-keV floor. Not a full late cosmology.')
    phase_bound()
    save('VERIFICATION_RESULTS.json',dict(groups=groups,full_kinetic_complete=False,physical_formation_kernel_derived=False,
        verdict='Numerical control checks only. Physical abundance and threshold remain null.'))
    print(json.dumps(groups,indent=2,default=lambda v:v.item()),flush=True)
    if not all(r['passed'] for r in groups.values()):raise SystemExit(1)

if __name__=='__main__':main()
