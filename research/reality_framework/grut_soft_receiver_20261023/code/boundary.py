#!/usr/bin/env python3
"""Independent free radial Klein-Gordon test of the actual velocity damper.

All lengths/times are in m_X units. This is a continuum linear benchmark of
the split damper, not a new nonlinear formation run or a reflection estimate
for an unknown physical particle distribution.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

ROOT=Path(__file__).resolve().parents[1]


def solve(k,width,eta,power=3,rtol=2e-11,flux=False):
    omega=np.hypot(1.,k)
    def rhs(x,y):
        damping=eta*(x/width)**power
        return [y[1],-(k*k+1j*omega*damping)*y[0]]
    result=solve_ivp(rhs,[width,0.],[0j,1+0j],method='DOP853',
        rtol=rtol,atol=rtol*.01,max_step=min(.5,width/100),dense_output=flux)
    assert result.success,result.message
    u,p=result.y[:,-1]
    incident=(u+p/(1j*k))/2
    reflected=(u-p/(1j*k))/2
    amplitude=reflected/incident
    out=dict(k_over_m=float(k),omega_over_m=float(omega),
        reflected_energy=float(abs(amplitude)**2),
        reflection_real=float(amplitude.real),reflection_imag=float(amplitude.imag))
    if flux:
        x=np.linspace(0,width,8193);field=result.sol(x)[0]/incident
        loss=omega/k*simpson(eta*(x/width)**power*abs(field)**2,x=x)
        out['normalized_flux_loss_integral']=float(loss)
        out['flux_balance_absolute_error']=float(abs(1-abs(amplitude)**2-loss))
    return out


def discrete(k,width,eta,h,power=3):
    """Independent second-order BVP, including a second-order incident Robin BC."""
    n=round(width/h);h=width/n;x=np.linspace(0,width,n+1)
    omega=np.hypot(1,k);matrix=lil_matrix((n+1,n+1),dtype=complex)
    rhs=np.zeros(n+1,dtype=complex)
    matrix[0,0]=-3+2j*k*h;matrix[0,1]=4;matrix[0,2]=-1;rhs[0]=4j*k*h
    for j in range(1,n):
        matrix[j,j-1]=1;matrix[j,j]=-2+h*h*(k*k+1j*omega*eta*(x[j]/width)**power)
        matrix[j,j+1]=1
    matrix[n,n]=1
    field=spsolve(matrix.tocsr(),rhs)
    return field[0]-1


def main():
    cases=[('original',90.,.03),('half_strength',90.,.015),
        ('old_larger_box',120.,.03),('candidate_wide',360.,.1),
        ('candidate_very_wide',720.,.1)]
    momenta=np.unique(np.r_[np.geomspace(.01,1.,81),[.01,.02,.03,.05,.1,.2,.5,1.]])
    rows=[]
    for name,width,eta in cases:
        samples=[solve(k,width,eta) for k in momenta]
        rows.append(dict(name=name,width=width,eta=eta,power=3,samples=samples,
            maximum_reflection_on_sampled_band=max(x['reflected_energy'] for x in samples),
            passes_one_percent_gate_on_sampled_band=all(x['reflected_energy']<.01 for x in samples)))
    controls=[]
    for k in [.02,.1,.5]:
        reference=solve(k,90,.03,flux=True);r=complex(reference['reflection_real'],reference['reflection_imag'])
        grids=[.5,.25,.125,.0625,.03125,.015625]
        errors=[float(abs(discrete(k,90,.03,h)-r)) for h in grids]
        orders=[float(np.log2(a/b)) for a,b in zip(errors,errors[1:])]
        tighter=solve(k,90,.03,rtol=2e-13)
        controls.append(dict(kind='independent_grid_and_tolerance',k=k,
            h=grids,complex_amplitude_errors=errors,observed_orders=orders,
            reflected_energy_tolerance_difference=abs(reference['reflected_energy']-tighter['reflected_energy']),
            flux_balance_absolute_error=reference['flux_balance_absolute_error']))
    for k in [.01,.1,1.]:
        free=solve(k,90.,0.)
        expected=-np.exp(2j*k*90)
        measured=complex(free['reflection_real'],free['reflection_imag'])
        controls.append(dict(kind='zero_damping_wall',k=k,complex_amplitude_error=float(abs(measured-expected))))
        slab=solve(k,90.,.03,power=0,flux=True)
        q=np.sqrt(k*k+1j*np.hypot(1,k)*.03);impedance=-q/np.tan(q*90)
        exact=(1j*k-impedance)/(1j*k+impedance)
        measured=complex(slab['reflection_real'],slab['reflection_imag'])
        controls.append(dict(kind='constant_damping_slab',k=k,
            complex_amplitude_error=float(abs(measured-exact)),
            flux_balance_absolute_error=slab['flux_balance_absolute_error']))
    # Continuum complex-stretch reference, NOT a time-domain implementation.
    # s(x)=1+i sigma(x)/omega; cubic sigma has integral sigma_max*D/4.
    pml=[dict(k=float(k),width=90.,sigma_max=40.,
        ideal_wall_reflected_energy=float(np.exp(-k/np.hypot(1,k)*40*90))) for k in [.01,.02,.03,.05,.1,.2,.5,1.]]
    out=dict(scope='Free exterior linear continuum receiver benchmark. Actual nonlinear reflected-population fractions are not calculated.',
        equation='U_tt=U_xx-U-eta(x) U_t; eta(x)=eta_max (x/D)^3; U(D)=0',
        gate='Reflected energy < 0.01 at every sampled k/m in [0.01,1]. Passing samples would still not certify an unsampled continuum.',
        cases=rows,controls=controls,matched_layer_reference=dict(
            scope='Analytic frequency-domain design reference only. No new time-domain PML, discretization, stability or nonlinear source certificate.',samples=pml))
    (ROOT/'BOUNDARY_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for row in rows:
        print(row['name'],'maximum R_energy',row['maximum_reflection_on_sampled_band'],
            'gate',row['passes_one_percent_gate_on_sampled_band'],flush=True)
    print('Controls',json.dumps(controls),flush=True)


if __name__=='__main__':main()
