#!/usr/bin/env python3
"""Independent O(epsilon^2) expansion of GRUT III v11 conformal power-law correlator.

The test is for the exact local worldline two-point function at fixed finite
complex separation, not the globally/windowed Wigner spectrum. No new GRUT law.
Run: python check_v11_powerlaw_expansion.py   (requires mpmath)
"""
import mpmath as mp

mp.mp.dps = 65


def exact_kernel(x, eps, H=1):
    """Conformal Wightman kernel; midpoint H, x=H*(t-t'), epsilon=1/p."""
    e = mp.mpf(eps)
    up, um = 1 + e*x/2, 1 - e*x/2
    return (-H**2 * (1-e)**2 /
            (4*mp.pi**2 * (up*um)**(1/e)
             * (up**(1-1/e) - um**(1-1/e))**2))


def de_sitter_kernel(x, H=1):
    """Here x=H_mid*Delta; argument H changes while H_mid is held fixed."""
    return -H**2 / (16*mp.pi**2 * mp.sinh(H*x/2)**2)


def first_order_coefficient(x):
    return x*mp.coth(x/2) - 2


def nonthermal_second_order_coefficient(x):
    return x**2/4 - x**3*mp.coth(x/2)/12


def test():
    xs = [mp.mpc(1, -1), mp.mpc(mp.mpf('0.75'), -mp.pi),
          mp.mpc(2, -mp.pi), mp.mpc(3, -2)]
    es = [mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]
    for x in xs:
        c1 = first_order_coefficient(x)
        c2 = nonthermal_second_order_coefficient(x)
        errors = []
        for eps in es:
            w = exact_kernel(x, eps)
            w0 = de_sitter_kernel(x)
            w_eff = de_sitter_kernel(x, 1-eps)
            first_est = (w/w0 - 1) / eps
            second_est = (w/w_eff - 1) / eps**2
            errors.append(abs(second_est-c2))
            # First-order and second-order coefficient errors both vanish
            # linearly in epsilon; tolerances accommodate all sampled x.
            assert abs(first_est-c1) < 50*eps, (x, eps, 'first')
            assert abs(second_est-c2) < 10*eps, (x, eps, 'second')
        assert errors[2] < errors[1] < errors[0], (x, errors)
        print('x=', mp.nstr(x,6),
              'analytic_second_order=', mp.nstr(c2,12),
              'error_at_eps=.005=', mp.nstr(errors[-1],6))
    print('PASS: first-order correction is an effective H shift;')
    print('      first non-temperature-shift shape appears at O(epsilon^2)')


if __name__ == '__main__':
    test()
