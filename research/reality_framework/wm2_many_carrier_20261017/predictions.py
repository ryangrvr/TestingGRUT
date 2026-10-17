"""Frozen tau=32 free-channel holdouts with exact rational series enclosures.

The enclosure certifies these mathematical probability evaluations, not D1's
all-size proof or a device. Independent quadrature and a 64-site propagation
are ordinary-double calculation controls. No experimental data are present.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from math import factorial, log, sqrt
import json
import resource
import time
import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm


def bessel_interval(order, z, last=24):
    assert Fraction(0) < z <= 1
    term = (z / 2) ** order / factorial(order)
    total = term
    for k in range(1, last + 1):
        term *= (z * z / 4) / (k * (order + k))
        total += term if k % 2 == 0 else -term
    omitted = term * (z * z / 4) / ((last + 1) * (order + last + 1))
    if (last + 1) % 2:
        return total - omitted, total
    return total, total + omitted


def outward(value, lower):
    mode = ROUND_FLOOR if lower else ROUND_CEILING
    with localcontext() as context:
        context.prec = 100
        context.rounding = mode
        result = Decimal(value.numerator) / Decimal(value.denominator)
        return str(result.quantize(Decimal('1e-50'), rounding=mode))


def run():
    started = time.perf_counter()
    tau = Fraction(32)
    out = {'status': 'CALCULATED_NOT_MEASURED',
           'target_time_tau': 32, 'prediction_charter_commit': '907272a7720a7fd8e99eb5b1fc5f67e77cf58e17',
           'certificate': 'exact rational alternating-series bounds, 25 terms, outward decimal endpoints',
           'independent_controls': 'ordinary double spectral quadrature and 64-site matrix propagation',
           'experimental_observations': 0, 'species': {}}
    for name, hopping in [('A', Fraction(1, 64)), ('B', Fraction(1, 96))]:
        z = 2 * hopping * tau
        intervals = []
        for r in range(3):
            low1, high1 = bessel_interval(r, z)
            low2, high2 = bessel_interval(r + 2, z)
            lo, hi = low1 + low2, high1 + high2
            assert lo > 0
            intervals.append((lo * lo, hi * hi))
        remainder = (1 - sum(hi for lo, hi in intervals),
                     1 - sum(lo for lo, hi in intervals))
        assert remainder[0] > 0
        intervals.append(remainder)
        H = np.diag(np.full(63, -float(hopping)), 1) + np.diag(np.full(63, -float(hopping)), -1)
        U = expm(-1j * 32 * H)
        physical = abs(U[:, 0]) ** 2
        finite = [*map(float, physical[:3]), float(sum(physical[3:]))]
        spectral = []
        quadrature_errors = []
        for r in range(3):
            integrand = lambda k: 2 / np.pi * np.sin((r + 1) * k) * np.sin(k)
            re, re_err = quad(lambda k: integrand(k) * np.cos(float(z) * np.cos(k)),
                              0, np.pi, epsabs=1e-13, epsrel=1e-13)
            im, im_err = quad(lambda k: integrand(k) * np.sin(float(z) * np.cos(k)),
                              0, np.pi, epsabs=1e-13, epsrel=1e-13)
            spectral.append(re * re + im * im)
            quadrature_errors.append(re_err + im_err)
        spectral.append(1 - sum(spectral))
        rows = []
        for event, (lo, hi), fp, integral in zip(['0', '1', '2', '3_or_more'], intervals, finite, spectral):
            midpoint = float((lo + hi) / 2)
            assert hi - lo < Fraction(1, 10 ** 60)
            assert abs(fp - midpoint) < 2e-14
            assert abs(integral - midpoint) < 2e-13
            rows.append({'recorded_position': event,
                         'probability_lower': outward(lo, True),
                         'probability_upper': outward(hi, False),
                         'rounded_probability': round(midpoint, 9),
                         'finite_chain_probability': fp,
                         'spectral_integral_probability': integral})
        # Ratio of successive positive Dyson terms is <=z/65 after order64.
        boundary = (z ** 64 / factorial(64)) / (1 - z / 65)
        claimed_bound = Fraction(801, 100) / 10 ** 90 if name == 'A' else Fraction(428, 100) / 10 ** 101
        assert boundary < claimed_bound
        out['species'][name] = {'hopping': str(hopping), 'z': str(z), 'probabilities': rows,
            'finite_chain_length': 64,
            'probability_boundary_error_upper': outward(boundary, False),
            'boundary_error_upper_scientific': float(boundary),
            'quadrature_estimated_amplitude_errors': quadrature_errors}
    M, alpha, channels = 10000, .01, 8
    radius = sqrt(log(2 * channels / alpha) / (2 * M))
    out['future_test_rule'] = {'independent_trials_per_species': M,
        'family_error_alpha': alpha, 'channels': channels,
        'Hoeffding_radius': radius,
        'required_pretest_relative_time_scale_error': .001,
        'required_pretest_dimensionless_Hamiltonian_norm_error': 1e-5,
        'required_pretest_preparation_trace_distance': .001,
        'required_pretest_readout_operator_error': .001,
        'calibration_values_are_measurements': False,
        'Hamiltonian_error_coefficient_with_time_calibration': 66,
        'A_maximum_systematic_budget_without_geometry': .00466,
        'B_maximum_systematic_budget_without_geometry': .003993333333333333,
        'decision': 'Reject if any heldout probability exceeds its fixed systematic plus family-wise statistical allowance; do not refit.'}
    elapsed = time.perf_counter() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    assert elapsed < 60 and peak < 1024
    out['elapsed_seconds'] = elapsed
    out['peak_resident_MiB'] = peak
    return out


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, allow_nan=False))
