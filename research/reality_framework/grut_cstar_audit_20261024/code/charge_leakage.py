#!/usr/bin/env python3
"""Leading forced radiation from explicitly specified U(1)-breaking probes.

This is not a quantum-gravity prediction for Wilson coefficients. The
unchanged Q=200000 stationary background is used only in weak response.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp,simpson

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_core_conversion_20261023'
sys.path.insert(0,str(PARENT/'code'))
from scattering import Background,G,V,MX

HBAR=6.582119569e-25
MRED=2.435e18
AGE=13.8e9*365.25*86400
QREF=200000.


def solve(bg,power,radius=40.,tol=2e-7,mesh=1001,scale=1.):
    nu=power*bg.omega
    frequencies=np.array([nu+bg.omega,nu-bg.omega,nu])
    k2=frequencies**2-np.array([1.,1.,6.])
    assert np.all(k2>0),'This implementation is scoped to all-open d>=5 probes.'
    momenta=np.sqrt(k2);r=np.linspace(0,radius,mesh)
    def rhs(x,y):
        z=y[:6]+1j*y[6:];potential=bg.potential(x,scale)
        potential[np.arange(3),np.arange(3)]-=k2[:,None]
        source=np.zeros_like(z[:3]);source[1]=x*bg.fields(x)[0]**(power-1)
        derivative=np.vstack([z[3:],np.einsum('ijn,jn->in',potential,z[:3])+source])
        return np.vstack([derivative.real,derivative.imag])
    def bc(left,right):
        l=left[:6]+1j*left[6:];z=right[:6]+1j*right[6:]
        residual=np.r_[l[:3],z[3:]-1j*momenta*z[:3]]
        return np.r_[residual.real,residual.imag]
    result=solve_bvp(rhs,bc,r,np.zeros((12,len(r))),tol=tol,bc_tol=tol*.02,max_nodes=40000)
    assert result.success,result.message
    outer=result.y[:3,-1]+1j*result.y[6:9,-1]
    amplitudes=outer*np.exp(-1j*momenta*radius)
    flux=momenta*abs(amplitudes)**2
    x=np.linspace(0,radius,16001);values=result.sol(x);u=values[:3]+1j*values[6:9]
    source=np.zeros_like(u);source[1]=x*bg.fields(x)[0]**(power-1)
    source_work=simpson(np.imag(np.sum(u.conj()*source,axis=0)),x=x)
    epsilon=power*2**(1-power/2)*V**(power-2)/(MX**2*MRED**(power-4))
    radial_fields=np.divide(u,x[None,:],out=np.zeros_like(u),where=x[None,:]>0)
    radial_fields[:,0]=values[3:6,0]+1j*values[9:12,0]
    number_rates=4*np.pi*MX/(G**2*HBAR)*epsilon**2*flux
    core_loss=np.dot(frequencies/bg.omega,number_rates)
    global_loss=power*sum(number_rates)
    torque_integral=4*np.pi/MX**3*simpson(x*x*(V*bg.fields(x)[0]/np.sqrt(2))**power,x=x)
    first_order_bound=4*torque_integral/(MRED**(power-4)*bg.omega*MX*QREF)
    max_coeff=np.sqrt(.01/(core_loss/QREF*AGE))
    out=dict(power=power,operator=f'c_{power} X^{power}/M_red^{power-4}+h.c.',radius=radius,tolerance=tol,initial_mesh=mesh,potential_scale=scale,
        nu=float(nu),frequencies_over_mX=frequencies.tolist(),momenta_over_mX=momenta.tolist(),
        amplitudes_real=amplitudes.real.tolist(),amplitudes_imag=amplitudes.imag.tolist(),
        unit_outgoing_wronskian_flux=flux.tolist(),unit_source_work=float(source_work),
        source_work_relative_residual=float(abs(sum(flux)-source_work)/max(sum(flux),1e-100)),
        maximum_BVP_residual=float(max(result.rms_residuals)),nodes=len(result.x),
        epsilon_for_c1=float(epsilon),number_rates_per_second_for_c1=number_rates.tolist(),
        maximum_unit_radial_field_amplitude=float(np.max(abs(radial_fields))),
        maximum_dimensionless_field_response_bound_for_c1=float(2*epsilon*np.max(abs(radial_fields))),
        core_charge_loss_per_second_for_c1=float(core_loss),global_charge_loss_per_second_for_c1=float(global_loss),
        initial_fractional_core_charge_loss_rate_per_second_for_c1=float(core_loss/QREF),
        frozen_profile_linear_charge_loss_fraction_over_reference_age_for_c1=float(core_loss/QREF*AGE),
        max_abs_coefficient_for_one_percent_initial_rate_age_target=float(max_coeff),
        bounded_first_order_charge_excursion_fraction_for_c1=float(first_order_bound),
        reference_age_seconds=AGE,reference_reduced_Planck_GeV=MRED,
        normalization='X=v psi/sqrt(2); psi_0=F exp(-i Omega tau); isotropic amplitudes, hence 4pi in flux. Source acts in anti-X channel; real-s channel uses sqrt(2) w.',
        scope='Ground-profile leading classical forced response. Rates scale as |c|^2. Unknown Wilson coefficients, hot cores, profile evolution and other quantum/multiparticle channels are not computed.')
    if scale==0:
        k=momenta[1];exact=-simpson(np.sin(k*x)*x*bg.fields(x)[0]**(power-1),x=x)/k
        out['free_amplitude_real_analytic']=float(exact)
        out['free_amplitude_relative_error']=float(abs(amplitudes[1]-exact)/max(abs(exact),1e-100))
        out['free_other_channel_amplitude_max']=float(max(abs(amplitudes[[0,2]])))
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    bg=Background(True);rows=[]
    settings=[(40.,2e-7,1001)] if args.pilot else [(30.,2e-7,801),(40.,2e-8,1401),(60.,2e-8,1801)]
    for power in [5,6]:
        for radius,tol,mesh in settings:
            row=solve(bg,power,radius,tol,mesh);rows.append(row)
            print(power,radius,'core rate',row['initial_fractional_core_charge_loss_rate_per_second_for_c1'],
                'age fraction',row['frozen_profile_linear_charge_loss_fraction_over_reference_age_for_c1'],
                'c target',row['max_abs_coefficient_for_one_percent_initial_rate_age_target'],flush=True)
    controls=[] if args.pilot else [solve(bg,power,40.,2e-8,1401,scale=0.) for power in [5,6]]
    background_controls=[] if args.pilot else [solve(Background(False),power,40.,2e-8,1401) for power in [5,6]]
    out=dict(scope='Explicit conditional symmetry-breaking probes, not quantum-gravity matching predictions or a population lifetime certificate.',
        parent_commit='ff49e5d0b4b40be1f7beb472c8bfa9d7070cdef2',background_profile_sha256=bg.hash,
        parent_scattering_sha256=hashlib.sha256((PARENT/'code/scattering.py').read_bytes()).hexdigest(),rows=rows,free_controls=controls,background_controls=background_controls)
    name='CHARGE_LEAKAGE_PILOT.json' if args.pilot else 'CHARGE_LEAKAGE_RESULTS.json'
    (ROOT/name).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
