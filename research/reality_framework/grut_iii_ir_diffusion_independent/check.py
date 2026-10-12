#!/usr/bin/env python3
"""GRUT III v13-v14 independently derived late-time IR coefficient.

No imports from GRUT_III; derives a second representation from the large-|w|
hypergeometric asymptotic and checks it against a finite-time derivative of
the full unexpanded two-point function. Uses the REAL, symmetrized part of
the analytically continued correlator (branch choice affects its imaginary
component). THIS IS NOT a detector-rate or generic-state KMS test.

Run: python check.py   (requires mpmath)
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 55
PI = mp.pi

def delta(eps):
    return eps/(1-eps)

def W_closed_real(eps, t, earlier):
    """Real part, analytic-continuation fixed, H(0)=a(0)=1."""
    d = delta(eps)
    H = lambda s: 1/(1+eps*s)
    a = lambda s: mp.power(1+eps*s,1/eps)
    tau = lambda s: -1/((1-eps)*a(s)*H(s))
    w = 1 + (tau(t)-tau(earlier))**2/(4*tau(t)*tau(earlier))
    expr = ((1-eps)**2*H(t)*H(earlier)*mp.gamma(3+d)*
            mp.gamma(-d)*mp.hyp2f1(3+d,-d,2,w)/(16*PI**2))
    return mp.re(expr)

def D_from_hypergeom_asymptotic(eps):
    """Independent derivation, using the leading (-w)^delta hypergeom term.

    2F1(3+d,-d;2;w) ~ Gamma(3+2d)/
    [Gamma(3+d)Gamma(2+d)]*(-w)^d.

    After multiplying H(t)H(t'), taking the real part across w>1, the
    late-t correlator is H(t')^2 * C(eps), with
    C=(1-eps)^2/(16pi^2)*Gamma(-d)*Gamma(3+2d)/
        Gamma(2+d)*4^(-d)*cos(pi*d).
    As dH(t')^2/dt'=-2eps H(t')^3, D=-2eps*C*H^3.
    """
    d = delta(eps)
    return (-eps*(1-eps)**2/(8*PI**2) * mp.gamma(-d) *
            mp.gamma(3+2*d)/mp.gamma(2+d) *
            mp.power(4,-d)*mp.cos(PI*d))

def D_gamma_reduced(eps):
    """Target expression tested, NOT used to derive the asymptotic answer."""
    d = delta(eps)
    return (mp.sqrt(PI)/(2*PI**3)*(1-eps)**3 *
            mp.gamma(mp.mpf(3)/2+d)*mp.gamma(1-d)*mp.cos(PI*d))

def run():
    rows=[]
    for label in ("0.05","0.1","0.2","0.3"):
        eps=mp.mpf(label)
        derived=D_from_hypergeom_asymptotic(eps)
        target=D_gamma_reduced(eps)
        independent=mp.diff(lambda early:W_closed_real(eps,mp.mpf(10000),early),
                            mp.mpf(0))
        algebra_rel=abs(derived/target-1)
        finite_rel=abs(independent/target-1)
        assert algebra_rel<mp.mpf("1e-48"), (eps,algebra_rel)
        assert finite_rel<mp.mpf("1e-10"), (eps,finite_rel)
        rows.append(dict(epsilon=label,
            derived=mp.nstr(derived,28),
            closed_gamma=mp.nstr(target,28),
            late_time_direct_derivative=mp.nstr(independent,28),
            relative_error_direct=mp.nstr(finite_rel,8)))
    # Derivative at epsilon=0 computed from the gamma form only AFTER the
    # independent late-time asymptotic establishes the formula.
    coeff=-3+mp.digamma(mp.mpf(3)/2)-mp.digamma(1)
    expected=-(1+2*mp.log(2))
    assert abs(coeff-expected)<mp.mpf("1e-50")
    assert abs(D_gamma_reduced(mp.mpf("1e-9")) /
               (1/(4*PI**2)) - (1+expected*mp.mpf("1e-9"))) < mp.mpf("1e-16")
    # Scope warning: this object cannot be a universally positive diffusion
    # variance for every constant epsilon. The formula has a first zero at 1/3.
    assert abs(D_gamma_reduced(mp.mpf(1)/3))<mp.mpf("1e-50")
    assert D_gamma_reduced(mp.mpf("0.3"))>0
    assert D_gamma_reduced(mp.mpf("0.4"))<0
    out=dict(
        assumptions=("power-law FRW; massless minimally coupled free scalar; "
                     "given analytic-continuation/IR prescription; real "
                     "late-time two-point function, H(0)=1; fixed epsilon"),
        late_time_controls=rows,
        first_order_coefficient=mp.nstr(coeff,35),
        negative_at_epsilon_0_4=mp.nstr(D_gamma_reduced(mp.mpf(".4")),20),
        valid_as_positive_noise_for_all_epsilon=False,
        exact_match_found_in_literature=False,
        warning=("Fixed mean-time spectral detailed balance is not global KMS "
                 "for a nonstationary state; phi-phi correlator and UDW "
                 "switching must be specified separately.")
    )
    print(json.dumps(out,indent=2))
    Path(__file__).with_name("results.json").write_text(
        json.dumps(out,indent=2)+"\n")
    print("GRUT_III_INDEPENDENT_IR_DIFFUSION_ASYMPTOTIC_PASS")

if __name__=="__main__":
    run()
