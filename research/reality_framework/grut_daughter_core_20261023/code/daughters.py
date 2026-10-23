#!/usr/bin/env python3
"""Point-specific daughter kinematics, cold annihilation and free streaming."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid, simpson
from scipy.optimize import brentq, minimize_scalar

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_portal_debris_window_20261023'
sys.path.insert(0,str(PARENT/'code'))
from widths import Widths, hidden, HBAR
from profiles import solve

G,LAM,LX,HD=.1,.03,.001,.01
V=.924089637375362
KAPPA=9.006401966916295e-10
MX=G*V;MS=np.sqrt(2*LAM)*V;MD=np.sqrt(G*G+HD*HD)*V
MH,VH=125.,246.
CM3S=.3893793721e-27*29979245800.
MPL=2.435e18
T0=2.7255*8.617333262e-14
GS0=2+7/8*6*(4/11)
H0=67.4*1000/3.085677581491367e22
OG=(np.pi**2*T0**4/15)/(3*MPL**2*(H0*HBAR)**2)
OR=OG*(1+(7/8)*(4/11)**(4/3)*3.046)
OM=.315;OL=1-OM-OR
C_MPC_S=299792458./3.085677581491367e22
YQ=4.8294009076988175e-9
ATTR=G**4/LAM-LX
WIDTHS=Widths()


class CSVSeed:
    def __init__(self,row):
        a=np.loadtxt(PARENT/row['profile_file'],delimiter=',',skiprows=1)
        self.x=a[:,0];self.p=np.array([row['Omega']])
        j=cumulative_trapezoid(self.x**2*a[:,1]**2,self.x,initial=0.)
        self.y=np.vstack([a[:,1:].T,j])
    def sol(self,r):return np.array([np.interp(r,self.x,y) for y in self.y])


def hubble(a):return H0*np.sqrt(OR/a**4+OM/a**3+OL)


def stream(aform,u0,aend=1.,radiation_only=False):
    def integrand(loga):
        a=np.exp(loga);u=u0*aform/a
        h=H0*np.sqrt(OR)/a**2 if radiation_only else hubble(a)
        return C_MPC_S*u/np.sqrt(1+u*u)/(a*h)
    return quad(integrand,np.log(aform),np.log(aend),epsabs=1e-13,epsrel=1e-10,limit=150)[0]


def decay_grid(gamma_s, mass=MX, points=96):
    """Unit-normalized isotropic two-body boosted daughter distribution."""
    c,w=np.polynomial.legendre.leggauss(points)
    beta=np.sqrt(1-1/gamma_s**2);estar=MS/2
    pstar=np.sqrt(estar*estar-mass*mass)
    energy=gamma_s*(estar+beta*pstar*c)
    return energy,np.sqrt(np.maximum(energy*energy-mass*mass,0))/mass,w/2


def ann(W, photon=10.):
    widths=WIDTHS;theta=KAPPA*V*VH/(MH*MH-MS*MS)
    gx,gd=hidden(V,LAM)
    gv=theta*theta*widths.unit(MS)
    prop=theta*(1/(W*W-MS*MS+1j*MS*(gx+gd+gv))
                -1/(W*W-MH*MH+1j*MH*.0041))
    unit=widths.unit(W,photon_coefficient=photon)
    return float(2*(2*G*G*V)**2*abs(prop)**2*unit/W)


def linear_case():
    n=YQ*(2*np.pi**2*10.75*.02**3/45)
    f2=n/MX
    for _ in range(10):f2=n/np.sqrt(MX*MX-ATTR*f2)
    omega=np.sqrt(MX*MX-ATTR*f2);phi2=V*V-G*G/LAM*f2
    aa=2*LX*f2;bb=2*LAM*phi2;cc=4*G**4*phi2*f2
    scale=(ATTR*f2/(2*omega))**2
    def growth(y):
        k2=ATTR*f2*y
        def det(z):
            zz=z*scale
            return (((k2+aa+zz)*(k2+zz)+4*omega**2*zz)*(k2+bb+zz)-cc*(k2+zz))/(ATTR**2*f2*f2*bb)
        if det(0)>=0:return 0.
        return np.sqrt(scale*brentq(det,0,4,xtol=1e-14))
    opt=minimize_scalar(lambda y:-growth(y),bounds=(.1,1.9),method='bounded')
    gamma=growth(opt.x);k=np.sqrt(ATTR*f2*opt.x)
    h=np.sqrt(np.pi**2*10.75/90)*.02**2/MPL
    return dict(nQ_GeV3=n,f_squared_GeV2=f2,omega_GeV=omega,
                max_growth_GeV=gamma,Gamma_over_H=gamma/h,
                growth_time_seconds=HBAR/gamma,k_peak_GeV=k,
                charge_per_wavelength_cube_DIAGNOSTIC=n*(2*np.pi/k)**3,
                supplied_coherent_state=True,formation_charge_distribution_computed=False)


def main():
    fixed=json.loads((PARENT/'FIXED_YIELD_RESULTS.json').read_text())
    row=fixed['rows'][1];pop=row['joint_screen_population'];profile=fixed['refinement']
    fx=pop['late_X_particles_per_object']*MX/pop['late_mass_per_object_GeV']
    fd=pop['late_D_particles_per_object']*MD/pop['late_mass_per_object_GeV']
    theta=KAPPA*V*VH/(MH*MH-MS*MS)
    gx,gd=hidden(V,LAM);gv=theta**2*Widths().unit(MS)
    point=dict(version='C1-HP-03/v1',parent_commit='d71cc012f0e27ab947ca942675b441fe8979f9f6',
               g_X=G,lambda_phi=LAM,lambda_X=LX,h_D=HD,m_D0_GeV=G*V,lambda_D=.001,
               v_phi_GeV=V,kappa=KAPPA,YQ_supplied=YQ,Q_test=200000,
               scale_status='Calibrated to the previous conditional late total matter; not predicted.',
               m_X_GeV=MX,m_s_GeV=MS,m_D_GeV=MD,mass_condition_s_to_XXstar=bool(MS>2*MX),
               binding_necessary_condition=bool(LAM*LX<G**4),
               theta=theta,Gamma_s_GeV=gx+gd+gv,tau_s_seconds=HBAR/(gx+gd+gv),
               scalar_lifetime_in_gv_time_units=MX/(gx+gd+gv),
               conditional_object_fraction=pop['object_halo_fraction'],conditional_free_X_fraction=fx,
               conditional_free_D_fraction=fd)
    (ROOT/'PARAMETER_POINT.json').write_text(json.dumps(point,indent=2)+'\n')

    # Independently solve adjacent integer charges without changing the parent.
    seed=CSVSeed(profile);energies={200000:profile['energy_per_v_GeV']*V};neighbor_rows=[]
    for q in [199999,200001]:
        sol=solve(q,LAM,seed,length=120,tolerance=2e-9)
        # save() belongs to the parent module: reconstruct quantities only.
        from profiles import quantities
        r=quantities(sol,q,LAM);r['tolerance']=2e-9
        energies[q]=r['energy_per_v_GeV']*V;neighbor_rows.append(r)
    removal=energies[200000]-energies[199999]
    addition=energies[200001]-energies[200000]
    available=MX+removal
    data=np.loadtxt(PARENT/profile['profile_file'],delimiter=',',skiprows=1)
    rho,F,Fp,S,Sp=data.T
    c=cumulative_trapezoid(rho*rho*F*F,rho,initial=0)
    r90=np.interp(.9*c[-1],c,rho)
    capture=dict(scope='Cold one-antiparticle capture and emission of ONE FREE on-shell s only; not all capture or internal excitations.',
                 neighboring_profiles=neighbor_rows,energies_GeV={str(k):v for k,v in energies.items()},
                 chemical_potential_from_omega_GeV=profile['Omega']*MX,
                 removal_energy_GeV=removal,addition_energy_GeV=addition,
                 energy_available_for_cold_capture_emission_GeV=available,
                 free_s_mass_GeV=MS,deficit_GeV=MS-available,
                 required_incident_antiX_total_energy_GeV=MS-removal,
                 required_incident_antiX_kinetic_energy_GeV=MS-removal-MX,
                 energy_available_if_antiX_is_born_from_s_at_rest_GeV=MS/2+removal,
                 one_free_s_emission_deficit_for_that_newborn_daughter_GeV=MS/2-removal,
                 cold_single_free_s_channel_open=bool(available>=MS),
                 test_profile_charge_radius90_gv_units=r90,
                 scalar_decay_time_over_light_crossing_R90=MX/(gx+gd+gv)/r90,
                 daughter_at_rest_s_kinetic_GeV=MS/2-MX)
    # Tree annihilation is a separate free-particle channel, not object capture.
    annihilation=[]
    for velocity in [0.,1e-8,1e-4,1e-3]:
        W=2*MX*np.sqrt(1+velocity**2)
        for photon in [0.,10.,100.]:
            rate=ann(W,photon)*CM3S
            injected=fx*fx*rate/(2*MX)
            annihilation.append(dict(COM_particle_p_over_m=velocity,photon_coefficient=photon,
                                     W_GeV=W,sigma_v_cm3_per_s=rate,
                                     total_injected_p_ann_cm3_per_s_per_GeV=injected,
                                     ratio_to_named_Planck_32e28_screen=injected/3.2e-28))
    peak_rows=[]
    for photon in [0.,10.,100.]:
        # An explicit light-pole control also covers newly born warm pairs
        # whose invariant mass can approach m_s. No UV/global rate bound.
        gamma_s=gx+gd+gv
        grid=np.unique(np.r_[np.linspace(2*MX,2.,3001),
                             MS+gamma_s*np.linspace(-20,20,2001),MS])
        grid=grid[(grid>=2*MX)&(grid<=2.)]
        rates=np.array([ann(w,photon) for w in grid])
        i=int(np.argmax(rates));fit=minimize_scalar(lambda w:-ann(w,photon),
            bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),method='bounded',
            options={'xatol':1e-15})
        sigma=max(ann(MS,photon),-fit.fun)*CM3S
        peak_rows.append(dict(photon_coefficient=photon,range_W_GeV=[2*MX,2.],
                              maximum_sigma_v_cm3_per_s=sigma,
                              max_p_injected_if_ALL_DM_were_free_symmetric_X=sigma/(2*MX),
                              ratio_to_named_Planck_screen_if_ALL_DM_free=sigma/(2*MX)/3.2e-28))
    annihilation_record=dict(scope='Leading free XX* annihilation into SM at supplied cold abundance; not trapped-core emission, early resonant network or a new CMB likelihood.',
                             symmetric_X_fraction=fx,complex_particle_factor_in_energy_injection=.5,
                             named_p_ann_screen_cm3_per_s_per_GeV=3.2e-28,
                             planck_screen_status='All injected SM energy counted; efficiency/redshift weighting is not recomputed.',
                             rows=annihilation,light_resonance_envelope=peak_rows)
    # Collisionless free-streaming, with the conditional decay spectrum.
    formations=[('T20MeV',T0/.02*(GS0/10.75)**(1/3)),
                ('T1MeV',T0/.001*(GS0/10.75)**(1/3)),
                ('T1keV',T0/1e-6),('z100000',1/100001),('z10000',1/10001)]
    streaming=[]
    for label,af in formations:
        for gamma in [1.,10.,1000.]:
            energies,us,weights=decay_grid(gamma)
            distances=np.array([stream(af,u) for u in us])
            to_z5=np.array([stream(af,u,1/6) for u in us])
            p_rec=us*af*1101
            late_mass=pop['late_mass_per_object_GeV']
            new_mass=row['mass_GeV']+(late_mass-row['mass_GeV'])/gamma
            streaming.append(dict(formation=label,a_formation=af,parent_gamma=gamma,
                                  mean_daughter_gamma_initial=float(np.dot(weights,energies/MX)),
                                  mean_free_stream_Mpc_today=float(np.dot(weights,distances)),
                                  mean_free_stream_Mpc_z5=float(np.dot(weights,to_z5)),
                                  max_free_stream_Mpc_today=float(max(distances)),
                                  mean_p_over_m_at_z1100=float(np.dot(weights,p_rec)),
                                  X_fraction_of_calibrated_DM=fx/gamma,
                                  X_fraction_of_total_remaining_matter=fx*late_mass/gamma/new_mass,
                                  total_matter_ratio_at_fixed_original_YQ=new_mass/late_mass,
                                  fixed_energy_source_X_number_ratio_to_cold_s=1/gamma))
    u0=np.sqrt((MS/(2*MX))**2-1)
    af=formations[0][1];aeq=OR/OM
    radiation_exact=C_MPC_S*u0*af/(H0*np.sqrt(OR))*(np.arcsinh(aeq/(u0*af))-np.arcsinh(1/u0))
    radiation_numeric=stream(af,u0,aeq,True)
    a_ann=T0/.0005*(GS0/10.75)**(1/3)
    gR_late=2*(1+(7/8)*(4/11)**(4/3)*3.046)
    highT_H_factor=np.sqrt(10.75/gR_late*(GS0/10.75)**(4/3))
    entropy_step_stream=(stream(af,u0,a_ann,True)/highT_H_factor
                         +stream(a_ann,u0*af/a_ann))
    sigma=2*np.sqrt(2*LAM)*V**3/3
    highT=LAM/4+G*G/6+HD*HD/12
    mstar=(LAM*LX)**.25*V
    mass_bulk_with_daughters=mstar+2*(MX-mstar)/MS*(pop['hidden_X_branch']*MX+pop['hidden_D_branch']*MD)
    output=dict(point=point,point_specific_daughter_kinematics=dict(
                    p_star_GeV=np.sqrt(MS*MS/4-MX*MX),gamma_star=MS/(2*MX),
                    beta_star=np.sqrt(1-4*MX*MX/(MS*MS)),kinetic_per_daughter_GeV=MS/2-MX),
                capture=capture,annihilation=annihilation_record,
                free_streaming=dict(scope='Collisionless propagation in the stated adopted background; not a Boltzmann/structure/Lyman-alpha fit.',
                                    H0_per_s=H0,Omega_m=OM,Omega_r=OR,Omega_L=OL,
                                    radiation_analytic_control_relative=radiation_numeric/radiation_exact-1,
                                    reference_20MeV_instant_entropy_step_sensitivity_Mpc=entropy_step_stream,
                                    entropy_step_scope='Illustrative instantaneous e+- entropy change at .5 MeV, not a BBN background integration.',
                                    rows=streaming),
                new_point_bulk_fixed_YQ=dict(scope='Conditional ground-state bulk limit; not the radial simulation population.',
                    energy_floor_per_charge_GeV=mstar,
                    objects_only_DM_ratio=YQ*mstar/4.4e-10,
                    objects_plus_slow_s_decay_rest_matter_DM_ratio=YQ*mass_bulk_with_daughters/4.4e-10,
                    assumed_scalar_source='All free-carrier-to-bulk release as slow vacuum s, followed by nonreacting daughters.'),
                new_point_neutral_wall=dict(tension_GeV3=sigma,cuberoot_MeV=1000*sigma**(1/3),
                                           ratio_to_named_stable_network_tension_bound=sigma/.00081**3,
                                           scope='Unoccupied neutral tree kink; standard surviving network only.'),
                new_point_high_temperature_curvature=dict(coefficient=highT,
                    conditional_zero_curvature_T_hidden_GeV=V*np.sqrt(LAM/highT),
                    scope='Leading high-T common hidden equilibrium, zero chemical potential and relativistic species; not a full transition or actual history.'),
                new_point_linear_onset=linear_case(),
                detector_source='Parent FIXED_YIELD_RESULTS.json second row: already recomputed at THIS point, not old benchmark values.')
    (ROOT/'DAUGHTER_RESULTS.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['annihilation','free_streaming']},indent=2))
    print('cold annihilation',annihilation[1])
    for r in streaming:
        if r['parent_gamma']==1:print('stream',r)


if __name__=='__main__':main()
