#!/usr/bin/env python3
"""Conditional phase-transfer diagnostic with prescribed radial/RD history."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad

ROOT=Path(__file__).resolve().parents[1]


def evolve(eta,ratio,delta0,Nend=8.,rtol=2e-10,atol=2e-13):
    # p=P/(H_i f_i^2), P=a^3 I delta_dot; eta=alpha6 f_i^4/(4 H_i^2).
    def rhs(N,z):
        delta,p=z
        return [p*(36*np.exp(N)+ratio**2*np.exp(-N)),eta*np.exp(-N)*np.sin(delta)]
    sol=solve_ivp(rhs,[0,Nend],[delta0,0.],method='DOP853',rtol=rtol,atol=atol,dense_output=True,max_step=.04)
    assert sol.success
    final=float(sol.y[1,-1])
    integral=quad(lambda N:np.exp(-N)*np.sin(sol.sol(N)[0]),0,Nend,epsabs=2e-11,epsrel=2e-10,limit=1000)[0]
    N=np.linspace(0,Nend,401);delta,p=sol.sol(N)
    # Reconstruct two independent phase rates in temporal gauge A0=0.
    theta_N=6*p*np.exp(N);sigma_N=-p*ratio**2*np.exp(-N)
    rel=6*theta_N-sigma_N
    residual=np.max(np.abs(rel-p*(36*np.exp(N)+ratio**2*np.exp(-N))))
    return dict(eta=eta,f_initial_over_Higgs=ratio,initial_relative_phase=delta0,N_end=Nend,
                p_final=final,comoving_X_charge_over_Hi_fi2=6*final,
                comoving_Sigma_charge_over_Hi_fi2=-final,
                remaining_absolute_p_change_bound=abs(eta)*np.exp(-Nend),
                source_integral=integral,source_integral_error=abs(final-eta*integral)/max(abs(eta),1e-100),
                reconstructed_phase_rate_error=float(residual),steps=len(sol.t),rtol=rtol,atol=atol)


def main():
    rows=[evolve(e,r,d) for e in [1e-6,1e-4,.01] for r in [.01,1.,6.,30.] for d in [0.,np.pi/2,-np.pi/2]]
    tight=evolve(.01,30.,np.pi/2,rtol=2e-12,atol=2e-15)
    coarse=next(r for r in rows if r['eta']==.01 and r['f_initial_over_Higgs']==30. and r['initial_relative_phase']>0)
    comparisons=[]
    for e in [1e-6,1e-4,.01]:
        base=next(r for r in rows if r['eta']==e and r['f_initial_over_Higgs']==.01 and r['initial_relative_phase']>0)
        for ratio in [1.,6.,30.]:
            full=next(r for r in rows if r['eta']==e and r['f_initial_over_Higgs']==ratio and r['initial_relative_phase']>0)
            comparisons.append(dict(eta=e,ratio=ratio,charge_ratio_to_large_reservoir=full['p_final']/base['p_final']))
    out=dict(scope='Conditional two-field phase EFT with prescribed f=f_i/a, h=constant, H=H_i/a^2. No radial, Friedmann, inflation, entropy or fragmentation solution.',
             definitions=dict(eta='alpha6 f_i^4/(4 H_i^2)',ratio='f_i/h',p='a^3 I delta_dot/(H_i f_i^2)',charge_yield='Y_X=6 H_i f_i^2 p_final/s_i; s_i must be supplied.'),
             rows=rows,tight_control=tight,tolerance_absolute_p_shift=abs(coarse['p_final']-tight['p_final']),
             finite_reservoir_comparisons=comparisons,net_gauge_charge_generated=0,
             initial_radial_amplitude_selected=False,initial_relative_phase_selected=False,entropy_history_selected=False,abundance_predicted=False)
    (ROOT/'PHASE_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('PHASE',len(rows),'conditional histories + tolerance control')


if __name__=='__main__':main()
