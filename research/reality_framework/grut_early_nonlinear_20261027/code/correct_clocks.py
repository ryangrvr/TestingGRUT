"""Repair source-clock diagnostics at the newly resolved warning samples."""
from en_common import *
from source_clocks import Clock

def main():
    old=T.read('RESPONSE_RESULTS.json')['rows'];rows=[]
    decay=T.K.SRC['mediator_debris']['cascade_vectors']['e1e-7']['Gamma_YH_to_YL_A_GeV_estimate']
    for r in read('DENSE_CORRECTIONS.json')['rows']:
        b=r['beta'];tag='quiet3' if b==3 else 'quiet5';prior=next(p for p in old if p['tag']==tag)
        data=np.load(ROOT/f'dense_b{b:g}.npz');times=data['times'];a=data['a']
        warning=r['preparation_phase_003_amplitude_warning'];tau=warning['tau_bracket'][1]
        gam=T.K.SRC['mediator_debris']['YL_three_body']['parents']['Gamma_L_GeV']/T.K.H_I
        clock=Clock(times,a,prior['late_C'],prior['matched_radiation_tau_offset'],gam)
        rows.append(dict(beta=b,preparation_phase_003_warning=warning,
            YL_rest_clock_number_fraction_at_warning=clock.number_fraction(tau),
            YL_prescribed_conformal_injection_fraction_at_warning=clock.energy_fraction(tau),
            YH_cascade_fraction_at_warning=float(-np.expm1(-decay/T.K.H_I*tau)),
            YH_cascade_direct_X_heating=False))
    save('CORRECTED_CLOCKS.json',dict(rows=rows,early_YH_cascade_Gamma_GeV=decay,
        scope='Fixed benchmark clocks with corrected dense warning times. Actual populations, in-medium widths and X-sector heat deposition are not recomputed.'))

if __name__=='__main__':main()
