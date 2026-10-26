"""Load-bearing mathematical/numerical controls, separate from physical verdict."""
import json,math
import numpy as np
from background import ROOT,MH,HH,MU,KAP,LX,LY,LS,write


def force_canonical(q,a):
    z=(q[0::2]+1j*q[1::2])/math.sqrt(2);x,y,s=z
    f=np.array([2*LX*abs(x)**2*x+3*KAP*y*np.conj(x)**2,
                2*LY*abs(y)**2*y+MH**2*a*a*y+KAP*x**3+2*MU*a*s*np.conj(y),
                2*LS*(abs(s)**2-a*a*HH*HH/2)*s+MU*a*y*y])
    out=np.empty(6);out[0::2]=math.sqrt(2)*f.real;out[1::2]=math.sqrt(2)*f.imag
    return out


def analytic_hessian(q,a):
    x,y,s=(q[0::2]+1j*q[1::2])/math.sqrt(2)
    A=np.zeros((3,3),complex);B=np.zeros((3,3),complex)
    A[0,0]=4*LX*abs(x)**2;B[0,0]=2*LX*x*x+6*KAP*y*np.conj(x)
    A[0,1]=3*KAP*np.conj(x)**2;A[1,0]=np.conj(A[0,1])
    A[1,1]=4*LY*abs(y)**2+MH*MH*a*a;B[1,1]=2*LY*y*y+2*MU*a*s
    A[1,2]=2*MU*a*np.conj(y);A[2,1]=np.conj(A[1,2])
    A[2,2]=4*LS*abs(s)**2-LS*a*a*HH*HH;B[2,2]=2*LS*s*s
    matrix=np.empty((6,6))
    matrix[0::2,0::2]=(A+B).real;matrix[0::2,1::2]=-(A-B).imag
    matrix[1::2,0::2]=(A+B).imag;matrix[1::2,1::2]=(A-B).real
    return matrix


def hamiltonian_checks():
    rng=np.random.default_rng(2026102503)
    T=np.zeros((6,6))
    for i,charge in enumerate([1.,3.,6.]):T[2*i,2*i+1]=-charge;T[2*i+1,2*i]=charge
    fdmax=wardmax=symmax=legendremax=constraintmax=0.
    for _ in range(18):
        a=float(rng.uniform(1,30));q=rng.normal(size=6)*[60,60,2,2,28*a,28*a]
        matrix=analytic_hessian(q,a);force=force_canonical(q,a)
        finite=np.empty((6,6))
        for j in range(6):
            step=np.zeros(6);h=1e-4*max(1.,abs(q[j]));step[j]=h
            finite[:,j]=(-force_canonical(q+2*step,a)+8*force_canonical(q+step,a)
                          -8*force_canonical(q-step,a)+force_canonical(q-2*step,a))/(12*h)
        fdmax=max(fdmax,float(np.max(abs(finite-matrix))/max(np.max(abs(matrix)),1)))
        g=T@q;gp=T@rng.normal(size=6);scale=max(float(np.linalg.norm(T@force)),1.)
        wardmax=max(wardmax,float(np.linalg.norm(matrix@g-T@force)/scale))
        symmax=max(symmax,float(np.max(abs(matrix-matrix.T))/max(np.max(abs(matrix)),1)))
        fluct=rng.normal(size=6);momentum=rng.normal(size=6);wave=.005;e=1e-7
        D=wave**2+e*e*(g@g);c=e*e/D;bar=g@momentum-gp@fluct
        velocity=momentum+e*e/wave**2*g*bar
        current=g@velocity-gp@fluct;A0=-e*current/D
        recovered=velocity+e*A0*g
        constraintmax=max(constraintmax,float(np.max(abs(recovered-momentum))))
        Om=matrix+wave**2*np.eye(6)
        lag=.5*(velocity@velocity-fluct@Om@fluct)-.5*c*current*current
        hamilton=.5*(momentum@momentum+fluct@Om@fluct)+.5*e*e/wave**2*bar*bar
        legendremax=max(legendremax,abs(float(momentum@velocity-lag-hamilton))/max(abs(float(hamilton)),1.))
    return dict(passed=fdmax<1e-7 and wardmax<1e-10 and symmax<1e-12 and legendremax<1e-12 and constraintmax<1e-10,
                canonical_hessian_fd_relative_error=fdmax,gauge_Ward_relative_error=wardmax,
                hessian_symmetry_relative_error=symmax,Legendre_relative_error=legendremax,
                Schur_momentum_recovery_absolute_error=constraintmax,samples=18)


def main():
    b=json.loads((ROOT/'BACKGROUND_RESULTS.json').read_text())
    v=json.loads((ROOT/'VECTOR_RESULTS.json').read_text())
    cc=json.loads((ROOT/'COUPLED_CONTROL.json').read_text())
    cb=json.loads((ROOT/'COUPLED_BAND.json').read_text())
    rel=json.loads((ROOT/'RELIC_REPLAY.json').read_text())
    budget=json.loads((ROOT/'BUDGET_RESULTS.json').read_text())
    bgmax=max(r['finish']['max_work'] for r in b['runs'])
    gaussmax=max(r['finish']['max_gauss'] for r in b['runs'])
    vals=[]
    for run in b['runs']:
        vals.append(next(w for w in run['windows'] if w['start']==30 and w['end']==40)['observables'][0])
    change=max(vals)/min(vals)-1
    groups={}
    groups['prepared_scalar_background']=dict(passed=bgmax<1e-8 and gaussmax<1e-9 and abs(change)<1e-8,
        work_max=bgmax,gauss_max=gaussmax,charge_tolerance_relative_spread=abs(change),
        conditional_gradient=b['preparation']['conditional_gradient_norm'])
    ab=max(r['absolute_Wronskian_error'] for r in v['modes']);scaled=max(r['scaled_Wronskian_error'] for r in v['modes'])
    changes=[]
    for start in range(0,len(v['controls']),4):
        ref,end,tight,long=v['controls'][start:start+4]
        changes.append(dict(k_over_kstar=ref['k_over_kstar'],
                            endpoint=abs(end['occupation']/ref['occupation']-1),
                            tolerance_and_start=abs(tight['occupation']/ref['occupation']-1),
                            longer_background=abs(long['occupation']/ref['occupation']-1)))
    groups['Proca_mode_controls']=dict(passed=ab<1e-5 and scaled<1e-8 and max(x['endpoint'] for x in changes)<.003
        and max(x['tolerance_and_start'] for x in changes)<1e-5 and max(x['longer_background'] for x in changes)<1e-4,
        absolute_Wronskian_max=ab,scaled_Wronskian_max=scaled,changes=changes)
    sp=v['spectra'];quad=abs(sp['actual_fine']['N_A_band']/sp['actual_coarse']['N_A_band']-1)
    profile=abs(sp['smooth']['N_A_band']/sp['actual_coarse']['N_A_band']-1)
    groups['band_integration_and_profile']=dict(passed=quad<.001 and profile<.001
        and all(s['negative_integral']==0 for s in sp.values()),quadrature_change=quad,smooth_profile_change=profile)
    groups['Coulomb_Hamiltonian_algebra']=hamiltonian_checks()
    comparison=abs(cb['N_A_common_band']/sp['actual_coarse']['N_A_band']-1)
    late_changes=[];tight_changes=[]
    loose,tight=cc['rows'][:2]
    for x in [.3,1.,3.]:
        aa=[p['late_common_mode_occupation'] for p in loose['projections'] if p['k_over_kstar']==x]
        bb=[p['late_common_mode_occupation'] for p in tight['projections'] if p['k_over_kstar']==x]
        late_changes.append(max(aa)/min(aa)-1)
        tight_changes.append(max(abs(c/d-1) for c,d in zip(aa,bb)))
    coarse_dict={round(p['k_over_kstar'],12):p for p in v['modes'][:97:2]}
    pointwise=max(abs(p['late_common_mode_occupation']/coarse_dict[round(p['k_over_kstar'],12)]['occupation']-1)
                  for p in cb['points'])
    groups['coupled_common_mode_transfer']=dict(passed=comparison<.005 and max(late_changes)<.003
        and max(tight_changes)<1e-5 and pointwise<.01,
        band_relative_to_Proca=comparison,max_matching_time_change=max(late_changes),
        tighter_common_mode_change=max(tight_changes),maximum_pointwise_change=pointwise,
        scope='Initial common mode only; late projection controls do not constitute all initial-mode production.')
    r0,r1=budget['rows'];dep_ratio=r1['deposit_R_per_vector']/r0['deposit_R_per_vector']
    nr=r1['N_A_band']/sp['actual_coarse']['N_A_band']
    groups['energy_deposit_and_gauge_scaling']=dict(passed=abs(dep_ratio*math.sqrt(10)-1)<1e-12
        and abs(nr*math.sqrt(10)-1)<.01,deposit_ratio=dep_ratio,number_ratio=nr,
        energy_ratio=r1['R_vector_band_only']/r0['R_vector_band_only'])
    er={r['kind']:r for r in rel['rows']};a=er['ke_plus_s'];z=er['ke_plus_s_refined']
    evolution_change=abs(a['f']/z['f']-1)
    groups['review_code_ke_replay']=dict(passed=evolution_change<.03 and abs(z['energy_work_cumulative_residual'])<.005,
        f_coarse=a['f'],f_refined=z['f'],fractional_grid_change=evolution_change,
        refined_cumulative_work_residual=z['energy_work_cumulative_residual'],
        source='Actual replay of reviewer f7d641a code, not an independent relic implementation.')
    result=dict(groups=groups,groups_passed=int(sum(bool(g['passed']) for g in groups.values())),groups_total=len(groups),
                physical_results=dict(specified_point_fails_corrected_proxy=bool(z['f']>.01),
                    finite_band_vector_R=r0['R_vector_band_only'],common_mode_coupled_R=cb['N_A_common_band']*r0['deposit_R_per_vector'],
                    original_isentropic_proxy_f=er['original_proxy']['f'],
                    full_kinetic_action_exclusion=False,new_viable_e_point_demonstrated=False,
                    complete_initial_physical_mode_covariance_evolved=False,
                    foundational_closure_prediction_derived=False))
    write(ROOT/'VERIFICATION_RESULTS.json',result)
    print(json.dumps(result,indent=2))
    if result['groups_passed']!=result['groups_total']:raise SystemExit(1)


if __name__=='__main__':main()
