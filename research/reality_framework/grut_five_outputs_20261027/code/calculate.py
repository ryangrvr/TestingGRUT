"""Five-output determination: matching, an ideal shape, and explicit remaining freedoms."""
from pathlib import Path
import json, math
import numpy as np
from scipy.integrate import simpson
from scipy.special import erf, spherical_jn
from scipy.optimize import minimize_scalar, brentq

ROOT=Path(__file__).resolve().parents[1]
RF=ROOT.parent
HEAVY=RF/'grut_heavy_vector_budget_20261025/INPUTS.json'
PROFILE=RF/'grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03_refined.csv'


def save(name, obj):
    (ROOT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def inputs():
    return json.loads(HEAVY.read_text())


def state_and_couplings():
    old=RF/'grut_whole_parts_20261027/HISTORY_MODEL.npz'
    with np.load(old) as d:
        embedding=d['embedding'];constraint=d['constraint'];gauss=d['global_neutrality']
    e=np.eye(2,dtype=complex)
    plus=(e[0]+e[1])/np.sqrt(2)
    kron=lambda *xs: __import__('functools').reduce(np.kron,xs)
    a=kron(e[0],e[1],e[0],plus)
    b=kron(e[1],e[0],e[0],plus)
    rows=[]
    for p,phase in ((.25,0),(.75,0),(.5,0),(.5,np.pi)):
        initial=np.sqrt(p)*a+np.sqrt(1-p)*np.exp(1j*phase)*b
        psi=embedding@initial
        rows.append({'p':p,'relative_phase':phase,'later_record_probability_plus_charge':p,
                     'X_charge':6*(2*p-1),'X_entropy_bits':-p*np.log2(p)-(1-p)*np.log2(1-p),
                     'global_purity':float(abs(np.vdot(psi,psi))**2),
                     'constraint_error':float(np.max(abs(constraint@psi))),
                     'neutrality_error':float(np.max(abs(gauss@psi))),
                     'joint_coherence':2*np.sqrt(p*(1-p))*np.cos(phase)})
    save('STATE_FAMILY.json',{'state_selected_by_closure':False,'scope':'finite physical model; not a cosmic state',
                              'blank_initial_records':True,'states':rows,'cosmic_preparation_law':None})

    par=inputs()['normalized_parameters']
    mass,mu,h,k=par['M_Y_GeV'],par['mu_GeV'],par['h_GeV'],par['kappa']
    r=np.sqrt(2)*mu*h/mass**2
    # Delta V = -B|X|^6 + (A X^6 + h.c.) in a real unitary-gauge convention.
    B=k*k/(mass*mass*(1-r*r)); A=r*B/2
    rng=np.random.default_rng(20261027)
    matrix_err=0.
    # Dimensionless comparison avoids subtracting vastly disparate powers of GeV.
    M2=1.; mixing=r/2
    H=np.diag([M2+2*mixing,M2-2*mixing])
    for _ in range(40):
        J=rng.normal()+1j*rng.normal()
        # y=(u+iv), V=(u,v)H(u,v)^T + 2(Re J, Im J).(u,v).
        force=np.array([J.real,J.imag])
        uv=-np.linalg.solve(H,force)
        direct=uv@H@uv+2*force@uv
        effective=(-abs(J)**2+mixing*(J*J+(J.conjugate())**2).real)/(1-r*r)
        matrix_err=max(matrix_err,abs(direct-effective)/max(abs(direct),1.))
    light=inputs()['light_sector']
    save('COUPLING_RELATIONS.json',{'scope':'quadratic tree matching in unitary gauge; same relation existed in gauge-bridge construction',
         'definition':'DeltaV=-B*abs(X)^6+(A*X^6+h.c.)',
         'r':float(r),'B_GeV_minus2':float(B),'A_GeV_minus2':float(A),
         'two_absA_over_B':float(2*abs(A)/B),'independent_matrix_error':float(matrix_err),
         'mX_GeV':light['mX_GeV'],'mphi_GeV':float(np.sqrt(2*light['lambdaPhi'])*light['vPhi_GeV']),
         'mphi_over_mX':float(np.sqrt(2*light['lambdaPhi'])/light['gX']),
         'dimensionless_couplings_selected':False,'full_RG_or_inflation_matching':False})


def abundance_and_vacuum():
    d=inputs();light=d['light_sector'];bench=d['benchmark']
    m=light['mX_GeV'];q=light['Q_reference'];Y=bench['YQ'];D=bench['DM_per_visible_entropy_GeV']
    epq=m*bench['original_E_over_Qm'];gap=m-epq;req=D/Y
    # Cold final species: positive monocharge cores, free net X, symmetric free pairs.
    # D_calc/Y = eta*epq + (1-eta)*m + 2*m*p, p=Ybar/Y_Q.
    eta_no_pairs=(m-req)/gap
    eta_max_primary=(m-req+.01*req)/gap
    eta_max_objects=(m-req)/(gap-.01*epq)
    pmax_primary=.01*D/(2*m*Y)
    rows=[]
    for eta in (0.,eta_no_pairs,.35,1.):
        base=Y*(eta*epq+(1-eta)*m)
        p_needed=(D-base)/(2*m*Y)
        rows.append({'supplied_core_charge_fraction':float(eta),'object_number_yield':float(eta*Y/q),
                     'object_mass_fraction_of_target':float(eta*Y*epq/D),
                     'base_mass_over_target_without_pairs':float(base/D),
                     'required_pair_number_per_net_charge_for_target':float(p_needed),
                     'required_pair_mass_fraction_of_target':float(2*m*Y*p_needed/D)})
    save('ABUNDANCE_LEDGER.json',{'scope':'necessary cold monocharge mass consistency, not predicted formation or final abundance',
        'assumptions':['all conserved light charge accounted in positive-charge objects or free X',
                       'one supplied object charge Q=200000 and its ground energy',
                       'no surviving scalar/vector reservoir, antiobjects, kinetic heat, or extra DM species',
                       'fixed supplied primordial yield; negligible late leakage'],
        'YQ_input':Y,'target_DM_per_visible_entropy_GeV':D,'mX_GeV':m,
        'object_Q_input':q,'object_mass_GeV':q*epq,'object_energy_per_charge_GeV':epq,
        'required_total_energy_per_net_charge_GeV':req,
        'primary_pair_allowance_denominator_GeV_per_entropy':D,
        'primary_pair_ratio_max':pmax_primary,
        'core_charge_fraction_if_no_pairs':eta_no_pairs,
        'core_charge_fraction_interval_at_1pct_total_DM_pairs':[eta_no_pairs,eta_max_primary],
        'secondary_core_charge_fraction_interval_at_1pct_actual_object_mass_pairs':[eta_no_pairs,eta_max_objects],
        'all_charge_in_reference_objects_mass_over_target':Y*epq/D,
        'minimum_YQ_multiplier_for_all_reference_objects_no_pairs':D/(Y*epq),
        'minimum_YQ_multiplier_for_all_reference_objects_at_1pct_total_DM_pairs':.99*D/(Y*epq),
        'conditional_rows':rows,'actual_core_charge_fraction':None,'actual_pair_ratio':None,
        'final_object_number_density':None,'final_f_sym':None})
    ref=2.52e-47;mp=2.435e18
    examples=[]
    for scale in (-1.,0.,1.):
        rho=scale*ref;matter=3*ref
        H=np.sqrt((matter+rho)/(3*mp*mp))
        examples.append({'vacuum_offset_GeV4':rho,'matter_density_GeV4':matter,
                         'H_GeV':float(H),'Friedmann_residual_GeV4':float(3*mp*mp*H*H-matter-rho)})
    hbar=6.582119569e-25;cm_per_GeV_inv=1.973269804e-14
    gev_to_g=1.7826619216278976e-24
    mass_density=ref*gev_to_g/cm_per_GeV_inv**3
    save('VACUUM_FREEDOM.json',{'scope':'observational reference conversion and allowed parameter family, not a vacuum-energy prediction',
        'rho_reference_GeV4_input':ref,'fourth_root_meV':ref**.25*1e12,
        'mass_density_g_cm3':mass_density,'M_reduced_GeV':mp,
        'rho_over_M4':ref/mp**4,'Lambda_GeV2':ref/mp**2,
        'examples_with_adjusted_initial_geometry':examples,
        'derivation':'V -> V+c leaves all field-potential derivatives unchanged; Tmunu -> Tmunu-c*gmunu',
        'renormalized_vacuum_parameter_selected':False,'predicted_rho_vacuum':None,
        'q_equilibrium_reference_does_not_supply_gap_or_future_exit':True})


def profile():
    p=np.loadtxt(PROFILE,delimiter=',',skiprows=1)
    rho,F,Fp,S,Sp=p.T
    d=inputs();l=d['light_sector'];g=l['gX'];v=l['vPhi_GeV'];q=l['Q_reference']
    omega=g*g*q/(4*np.pi*simpson(rho*rho*F*F,x=rho))
    density=.5*(Fp*Fp+Sp*Sp)+.5*omega**2*F*F+l['lambdaPhi']/g**2/4*(S*S-1)**2+.5*S*S*F*F+l['lambdaX']/g**2/4*F**4
    energy=4*np.pi*simpson(rho*rho*density,x=rho)*v/g
    weight=rho*rho*(1-S*S);norm=simpson(weight,x=rho)
    def form(q_GeV):
        arg=np.atleast_1d(q_GeV)/(g*v)
        return simpson(weight[None,:]*np.sinc(arg[:,None]*rho[None,:]/np.pi),x=rho,axis=1)/norm
    return energy,form,float(omega)


def inverse_speed(vmin,v0=220.,ve=232.,vesc=544.):
    x=np.asarray(vmin)/v0;y=ve/v0;z=vesc/v0
    nz=erf(z)-2*z*np.exp(-z*z)/np.sqrt(np.pi)
    a=erf(x+y)-erf(x-y)-4*y*np.exp(-z*z)/np.sqrt(np.pi)
    b=erf(z)-erf(x-y)-2*(z+y-x)*np.exp(-z*z)/np.sqrt(np.pi)
    return np.maximum(0.,np.where(x<z-y,a,np.where(x<z+y,b,0.))/(2*nz*ve))


def recoil():
    mass,form,omega=profile()
    # Fixed ideal grid. It is a shape diagnostic, with no exposure or event prediction.
    edges=np.linspace(5,40,141);centers=(edges[1:]+edges[:-1])/2
    nodes,weights=np.polynomial.legendre.leggauss(8)
    energies=centers[:,None]+np.diff(edges)[:,None]*nodes/2
    targets=[('Xe131',131),('Ar40',40),('Ge74',74)]
    eflat=energies.ravel();width=np.diff(edges)
    objects=[];kinematics=[];grids=[]
    for name,a in targets:
        ma=a*.93149410242
        q=np.sqrt(2*ma*eflat*1e-6)
        rn=np.sqrt((1.23*a**(1/3)-.6)**2+7*np.pi**2*.52**2/3-5*.9**2)
        x=q*rn/.1973269804
        nuclear=(3*spherical_jn(1,x)/x*np.exp(-.5*(q*.9/.1973269804)**2))**2
        scalar=(125.2**2/(125.2**2+q*q))**2
        oq=form(q)**2
        kinematics.append((name,ma,nuclear*scalar,oq))

    def shape(entry,dm_mass,object_form,halo):
        name,ma,nuc,fo=entry
        reduced=dm_mass*ma/(dm_mass+ma)
        vmin=299792.458*np.sqrt(ma*eflat*1e-6/(2*reduced*reduced))
        r=nuc*inverse_speed(vmin,*halo)*(fo if object_form else 1.)
        bins=width/2*np.sum(r.reshape(energies.shape)*weights[None,:],axis=1)
        return bins/bins.sum() if bins.sum()>0 else np.zeros_like(bins)

    def divergence(p,q):
        mask=p>0
        return float(np.sum(p[mask]*np.log(p[mask]/np.maximum(q[mask],1e-300))))

    def fit(ps,indices,halo):
        def objective(logm):
            return float(np.mean([divergence(ps[i],shape(kinematics[i],np.exp(logm),False,halo)) for i in indices]))
        grid=np.linspace(np.log(3),np.log(1e6),180)
        vals=np.array([objective(t) for t in grid]);j=int(np.argmin(vals))
        result=minimize_scalar(objective,bounds=(grid[max(0,j-1)],grid[min(len(grid)-1,j+1)]),method='bounded',options={'xatol':1e-11})
        return float(np.exp(result.x)),float(result.fun)

    main=None;controls=[];saved={}
    for v0,vesc in [(v0,vesc) for v0 in (200.,220.,240.) for vesc in (500.,544.,580.)]:
        halo=(v0,232.,vesc)
        ps=[shape(entry,mass,True,halo) for entry in kinematics]
        joint,kl=fit(ps,list(range(3)),halo)
        fits=[fit(ps,[i],halo) for i in range(3)]
        rows=[]
        for i,(name,ma,nuc,fo) in enumerate(kinematics):
            qs=shape(kinematics[i],joint,False,halo)
            rows.append({'target':name,'separate_best_WIMP_mass_GeV':fits[i][0],
                         'separate_KL_per_event':fits[i][1],
                         'KL_at_joint_mass':divergence(ps[i],qs),
                         'TV_at_joint_mass':float(np.sum(abs(ps[i]-qs))/2)})
            if v0==220 and vesc==544:
                saved[name+'_object_probability']=ps[i];saved[name+'_point_probability_at_joint_mass']=qs
        result={'halo_km_s':list(halo),'common_best_point_mass_GeV':joint,
                'equal_target_weighted_KL_per_signal_event':kl,'targets':rows}
        controls.append(result)
        if v0==220 and vesc==544: main=result
    np.savez_compressed(ROOT/'RECOIL_SHAPES.npz',edges_keV=edges,**saved)
    qs=np.array([0.,.02,.05,.08,.1])
    # Secondary, post-primary nuisance challenge: fix the object's baseline halo,
    # allow the comparator its own coarse-grid halo, and refit its common mass.
    baseline=[shape(entry,mass,True,(220.,232.,544.)) for entry in kinematics]
    adversarial=[]
    for cv0 in (200.,220.,240.):
        for cve in (220.,232.,244.):
            for cesc in (500.,544.,580.):
                cm,ckl=fit(baseline,list(range(3)),(cv0,cve,cesc))
                adversarial.append({'comparator_halo_km_s':[cv0,cve,cesc],
                                    'best_common_mass_GeV':cm,'KL_per_event':ckl})
    best_halo=min(adversarial,key=lambda x:x['KL_per_event'])
    first_node=brentq(lambda q:float(form(np.array([q]))[0]),.06,.09,xtol=1e-14)
    node_energies=[{'target':name,'ER_keV':float(first_node**2/(2*ma)*1e6),
                    'inside_primary_window':bool(5<=first_node**2/(2*ma)*1e6<=40)}
                   for name,ma,nuc,fo in kinematics]
    save('OBSERVABLE_PREDICTION.json',{
        'scope':'conditional ideal normalized monocharge recoil shapes, not a viable cosmology or detector likelihood',
        'object_Q_supplied':200000,'profile_mass_GeV':mass,'profile_Omega':omega,
        'reference_mass_relative_difference':mass/(inputs()['benchmark']['original_E_over_Qm']*inputs()['light_sector']['mX_GeV']*200000)-1,
        'window_keV':[5,40],'ideal_monoisotopes':['Xe131','Ar40','Ge74'],
        'separate_target_normalizations_free':True,'primary':main,'shared_halo_controls':controls,
        'secondary_comparator_halo_scan':{'object_halo_fixed':[220.,232.,544.],
             'timing':'added as a nuisance challenge after viewing the primary result',
             'grid_points':27,'best':best_halo,'full_nuisance_likelihood':False},
        'source_form_factor_squared':[{'q_GeV':float(q),'F_squared':float(f*f)} for q,f in zip(qs,form(qs))],
        'first_source_form_factor_zero_GeV':float(first_node),
        'first_zero_target_energies':node_energies,
        'zero_scope':'ideal elastic single-Q static profile; mixtures, resolution, inelastic channels and backgrounds can fill it',
        'future_data_or_significance_calculated':False,'absolute_event_rate':None,
        'portal_thermal_history_certified':False,'new_GRUT_specific_law':False})
    print(json.dumps({'profile_mass_GeV':mass,'joint_best_mass_GeV':main['common_best_point_mass_GeV'],
                       'joint_KL':main['equal_target_weighted_KL_per_signal_event']}))


if __name__=='__main__':
    state_and_couplings()
    abundance_and_vacuum()
    recoil()
