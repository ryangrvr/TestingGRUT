#!/usr/bin/env python3
"""Short full homogeneous X,Y,Sigma release, preserving the compact Gauss law.

The metric is prescribed radiation domination after a supplied Hubble-mass
release. Phi is held at zero in this high-field regime. No late yield is inferred.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]


def evolve(point,angle,end=100.,rtol=2e-10,atol=2e-12):
    p=point['parameters'];inf=point['inflation'];H=inf['assumed_H_inf_GeV']
    M=p['M_GeV']/H;hh=p['Higgs_amplitude_over_M']*M;mu=p['mu_over_M']*M
    k=p['kappa'];lx=p['lambda_X'];ly=p['lambda_Y'];ls=p['lambda_Sigma']
    amp=np.array(inf['full_radial_minimum']['amplitudes_over_H'])
    initial=np.array([amp[0]*np.exp(-1j*angle),-amp[1],-amp[2]],dtype=complex)
    def potential(z):
        x,y,s=z
        return lx*abs(x)**4+ly*abs(y)**4+ls*(abs(s)**2-hh*hh/2)**2+M*M*abs(y)**2+2*np.real(k*np.conj(y)*x**3+mu*np.conj(s)*y*y)
    def rhs(t,z):
        x,y,s=z[:3];vx,vy,vs=z[3:6];expansion=1/(1+2*t)
        force=np.array([2*lx*abs(x)**2*x+3*k*y*np.conj(x)**2,
            M*M*y+2*ly*abs(y)**2*y+k*x**3+2*mu*s*np.conj(y),
            2*ls*(abs(s)**2-hh*hh/2)*s+mu*y*y])
        vel=z[3:6];work=6*expansion*float(np.sum(abs(vel)**2))
        return np.r_[vel,-3*expansion*vel-force,complex(work)]
    z0=np.r_[initial,np.zeros(4,dtype=complex)]
    t=np.linspace(0,end,2001)
    sol=solve_ivp(rhs,[0,end],z0,t_eval=t,method='DOP853',rtol=rtol,atol=atol,max_step=.05)
    assert sol.success
    positions=sol.y[:3];velocity=sol.y[3:6];a=np.sqrt(1+2*t)
    charge=-2*np.imag(np.conj(positions)*velocity)*a**3
    gauge=np.array([1.,3.,6.])@charge
    kin=np.sum(abs(velocity)**2,axis=0)
    energy=kin+np.array([potential(sol.y[:3,i]) for i in range(len(t))])
    ledger=(energy+sol.y[6].real-energy[0])/energy[0]
    scale=max(float(np.max(np.sum(abs(charge)*np.array([1.,3.,6.])[:,None],axis=0))),1.)
    return dict(initial_X_phase_angle=angle,end_Hrel_times_t=end,scale_factor_end=float(a[-1]),rtol=rtol,atol=atol,
        initial_energy_over_H4=float(energy[0]),
        final_comoving_species_charge_over_H3=charge[:,-1].tolist(),
        max_abs_comoving_X_charge_over_H3=float(np.max(abs(charge[0]))),
        max_relative_Gauss_charge=float(np.max(abs(gauge))/scale),
        max_relative_scalar_stress_ledger_error=float(np.max(abs(ledger))),
        samples=dict(t=t[::20].tolist(),comoving_X_charge=charge[0,::20].tolist(),comoving_Y_charge=charge[1,::20].tolist(),comoving_Sigma_charge=charge[2,::20].tolist()))


def main():
    point=json.loads((ROOT/'POINT_RESULTS.json').read_text())
    rows=[evolve(point,a) for a in [0.,np.pi/24,-np.pi/24]]
    fine=evolve(point,np.pi/24,rtol=2e-12,atol=2e-14)
    mirror=max(abs(np.array(rows[1]['final_comoving_species_charge_over_H3'])+np.array(rows[2]['final_comoving_species_charge_over_H3'])))
    shift=max(abs(np.array(rows[1]['final_comoving_species_charge_over_H3'])-np.array(fine['final_comoving_species_charge_over_H3'])))
    out=dict(scope='Short full homogeneous three-field classical diagnostic, supplied radiation metric and Hubble-mass release, phi=0. Not inflationary fluctuations, a final asymmetry, or a cosmological history.',
        rows=rows,tight_control=fine,CP_conjugate_charge_residual=float(mirror),tolerance_final_charge_shift=float(shift),
        no_source_statement='Real constant couplings preserve the real-field subspace. Own-potential phase alignment and zero initial angular momenta yield exactly zero X,Y,Sigma charges for arbitrary real radial evolution.',
        missing_phase_mechanism='The inflaton-dependent phase mismatch assumed by the review is absent in 7737487. A selector/mismatch interaction or CP-odd preparation is an additional physical input.',
        a_late_charge_yield_has_not_been_computed=True)
    (ROOT/'RELEASE_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('FULL_RELEASE_ALIGNED_Q',rows[0]['final_comoving_species_charge_over_H3'])
    print('FULL_RELEASE_MISALIGNED_Q',rows[1]['final_comoving_species_charge_over_H3'])
    print('GAUSS',max(r['max_relative_Gauss_charge'] for r in rows),'STRESS',max(r['max_relative_scalar_stress_ledger_error'] for r in rows),'TOL_SHIFT',shift)


if __name__=='__main__':main()
