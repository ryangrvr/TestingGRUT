"""New late-response tolerance and radiation-clock checks."""
from timing_common import *
from source_clocks import energy_tail,Clock

def main():
    groups=read('RESPONSE_CHECKS.json')['groups'];exe=build();rows=[]
    for beta,decay in [(3.,.5),(5.,1.)]:
        times=np.arange(0,240.0001,.5)
        a,v,_=integrate(exe,f'tight_b{beta:g}',beta,decay,[2.13],times,rtol=2e-12)
        b,w,_=integrate(exe,f'normal_b{beta:g}',beta,decay,[2.13],times,rtol=2e-10)
        av,_,_=responses(beta,a,v);bw,_,_=responses(beta,b,w)
        err=float(np.max(abs(av-bw)/(1+abs(av))))
        aa=np.exp(a[:,15]);ca=crossings(times,aa,.03*av,[2.13],1.)
        cb=crossings(times,np.exp(b[:,15]),.03*bw,[2.13],1.)
        same=(ca is None and cb is None) or (ca is not None and cb is not None and ca['tau_bracket']==cb['tau_bracket'])
        rows.append(dict(beta=beta,response_scaled_error=err,crossing_bracket_same=bool(same)))
    groups['late_tolerance_and_crossing']=dict(passed=max(r['response_scaled_error'] for r in rows)<1e-5 and all(r['crossing_bracket_same'] for r in rows),rows=rows)
    gamma_test=.1;C=.8;offset=-.7
    direct=quad(lambda u:np.sqrt(2*C/gamma_test*(u-gamma_test*offset))*np.exp(-u),0,np.inf,epsabs=1e-12,epsrel=2e-12)[0]
    exact=energy_tail(0,gamma_test,C,offset);error=abs(exact/direct-1)
    legacy=[]
    for key in ['e1e-7_band','e1e-6_band']:
        x=K.SRC['vector_sources'][key];g=x['parent']['gamma']
        for q in x['quantiles']:
            fraction=q['number_fraction'];want=q['analytic_a']
            got=np.sqrt(1-2*np.log1p(-fraction)/g)
            legacy.append(float(abs(got/want-1)))
    assert len(legacy)==10
    refinements=[]
    for row in read('RESPONSE_RESULTS.json')['rows']:
        if not row['tag'].startswith('quiet'):continue
        z=np.load(ROOT/(row['tag']+'.npz'));a=np.exp(z['bg'][:,15]);gam=K.SRC['mediator_debris']['YL_three_body']['parents']['gamma_L']
        clock4=Clock(z['times'],a,row['late_C'],row['matched_radiation_tau_offset'],gam,order=4)
        clock8=Clock(z['times'],a,row['late_C'],row['matched_radiation_tau_offset'],gam,order=8)
        refinements.append(float(abs(clock4.energy_early(134.)/clock8.energy_early(134.)-1)))
    groups['decay_energy_clock_integral']=dict(passed=error<1e-10 and max(legacy)<1e-10 and max(refinements)<1e-10,
        analytic_vs_independent_quadrature=error,legacy_quantile_max_relative=max(legacy),
        piecewise_Gauss4_vs8_energy_relative=max(refinements))
    result=read('RESPONSE_RESULTS.json');checks=[]
    for r in result['rows']:
        parts=r['controls']['check'].split();f=float(parts[3]);g=float(parts[4])
        checks.append(dict(tag=r['tag'],Friedmann=f,Gauss=g,C_late_range=r['late_C_relative_range']))
    groups['extended_background_constraints']=dict(passed=all(r['Friedmann']<1e-6 and r['Gauss']<1e-5 and r['C_late_range']<1e-6 for r in checks),rows=checks)
    quartic=read('QUARTIC_CHECKS.json');groups['auxiliary_quartic_extended_band']=quartic
    save('VERIFICATION_RESULTS.json',dict(groups=groups,independent_reviewer_ran=False,physical_f_sym=None))
    print((ROOT/'VERIFICATION_RESULTS.json').read_text(),flush=True)
    assert all(g['passed'] for g in groups.values())

if __name__=='__main__':main()
