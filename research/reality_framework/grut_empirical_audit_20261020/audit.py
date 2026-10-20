"""Same-author rederivation, archive integrity and retrospective distance test.

No old prediction routine is imported. No scalar evolution is executed.
Commands: audit | predict | compare | verify. See PROTOCOL.json for decisions.
"""
import csv
from decimal import Decimal, localcontext
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import random
import subprocess
import sys
import tarfile

P = Path(__file__).resolve().parent
ROOT = P.parents[2]
BASE = P.parent / 'grut_3d_empirical_20261020'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, data):
    p.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inputs():
    protocol = read(P / 'PROTOCOL.json')
    return protocol, read(P / protocol['calibration_file'])


def guard():
    frozen = read(P / 'FREEZE.json')
    public = read(P / 'PUBLIC_FREEZE.json')
    require(sha(P / 'FREEZE.json') == public['freeze_sha256'], 'Freeze changed')
    for name, expected in frozen['source_sha256'].items():
        require(sha(P / name) == expected, 'Frozen source changed: ' + name)
        relative = (P / name).relative_to(ROOT).as_posix()
        raw = subprocess.check_output(['git', 'show', public['commit'] + ':' + relative], cwd=ROOT)
        require(hashlib.sha256(raw).hexdigest() == expected, 'Public source differs: ' + name)
    for name, expected in frozen['baseline_input_sha256'].items():
        require(sha(ROOT / name) == expected, 'Baseline input changed: ' + name)
        raw = subprocess.check_output(['git', 'show', frozen['base_commit'] + ':' + name], cwd=ROOT)
        require(hashlib.sha256(raw).hexdigest() == expected, 'Baseline not pinned: ' + name)
    return frozen, public


def parameters(protocol):
    d = protocol['parameters']
    return d, min(1.0, (d['a'] * d['b']) ** 0.25)


def floor_values(protocol, calibration, conservative=False):
    d, mu = parameters(protocol)
    n = protocol['sensitivity']['width_multiplier'] if conservative else 0.0
    wb = calibration['omega_b'] - n * calibration['omega_b_sigma']
    mb = calibration['mean_mass_per_baryon_GeV_convention']
    if conservative:
        mb *= 1 + calibration['mean_mass_systematic_fraction']
    wh = wb * mu * d['mX_GeV'] * abs(d['r_YQ0_over_YB0']) / mb
    return {'omega_b': wb, 'baryon_mass_GeV': mb, 'omega_hidden_floor': wh,
            'H_floor_km_s_Mpc': 100 * math.sqrt(wb + wh)}


def envelope(z, H, c):
    # expm1 avoids loss of significance at small z.
    return 2 * c / H / (1 + z) * (-math.expm1(-0.5 * math.log1p(z)))


def decimal_envelope(z, H, c):
    with localcontext() as context:
        context.prec = 60
        z, H, c = (Decimal(str(x)) for x in (z, H, c))
        return float(2 * c / H / (1 + z) * (1 - 1 / (1 + z).sqrt()))


def simpson_integral(z, steps):
    h = z / steps
    vals = [(1 + i * h) ** -1.5 for i in range(steps + 1)]
    return h / 3 * math.fsum([vals[0], vals[-1],
                             4 * math.fsum(vals[1:-1:2]),
                             2 * math.fsum(vals[2:-1:2])])


def baseline_integrity():
    sm = read(BASE / 'SOURCE_MANIFEST.json')
    source_count = 0
    for name, item in sm['files'].items():
        p = BASE / name
        require(p.stat().st_size == item['bytes'] and sha(p) == item['sha256'],
                'Old packet file changed: ' + name)
        source_count += 1
    for name, expected in sm['preserved_parents'].items():
        require(sha(ROOT / name) == expected, 'Old parent changed: ' + name)
    rm = read(BASE / 'RAW_MANIFEST.json')
    found = {}
    for archive, expected in rm['archives'].items():
        p = BASE / archive
        require(p.stat().st_size == expected['bytes'] and sha(p) == expected['sha256'],
                'Old raw archive changed: ' + archive)
        with tarfile.open(p, 'r:gz') as tar:
            for item in tar:
                require(item.isfile() and item.name.startswith('raw/'), 'Unexpected archive entry')
                name = item.name[4:]
                require(name in rm['raw_files'] and name not in found, 'Missing/duplicate raw key')
                data = tar.extractfile(item).read()
                actual = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
                require(actual == rm['raw_files'][name], 'Old archived raw changed: ' + name)
                found[name] = actual
    require(set(found) == set(rm['raw_files']), 'Raw manifest membership mismatch')
    return {'packet_files_verified': source_count, 'parent_files_verified': len(sm['preserved_parents']),
            'archives_verified': len(rm['archives']), 'raw_files_verified': len(found),
            'raw_bytes_verified': sum(x['bytes'] for x in found.values()),
            'old_integrity_check_count': read(BASE / 'VERIFICATION.json')['count'],
            'independent_scientific_replication': False}


def numerical_audit():
    protocol, cal = inputs()
    d, mu = parameters(protocol)
    m = d['mX_GeV']
    g = d['gX']
    v = d['v_GeV']
    checks = {}
    old_point = read(P / protocol['parameter_file'])['empirical_point']
    checks['parameters_preserved'] = all(math.isclose(d[k], old_point[k], rel_tol=1e-14)
                                         for k in ['a', 'b', 'gX', 'v_GeV', 'mX_GeV'])
    checks['physical_dimensionless_couplings'] = (
        math.isclose(d['lambda_phi'] / g**2, d['a'], rel_tol=1e-14) and
        math.isclose(d['lambda_X'] / g**2, d['b'], rel_tol=1e-14))
    checks['mX_equals_gv'] = math.isclose(m, g*v, rel_tol=1e-14)
    physical_mstar = v * min(g, (d['lambda_phi'] * d['lambda_X'])**0.25)
    checks['dimensionful_floor_route'] = math.isclose(physical_mstar, mu*m, rel_tol=1e-14)
    # A rotating canonical state, with independently chosen phase, checks sign and factor two.
    angle, f, omega = .381, 2.3, .77
    q1, q2 = f*math.cos(angle), -f*math.sin(angle)
    dq1, dq2 = omega*q2, -omega*q1
    charge = q2*dq1 - q1*dq2
    checks['canonical_charge_normalization'] = math.isclose(charge, omega*f*f, rel_tol=1e-14)
    checks['physical_E_Q_conversion'] = math.isclose((v/g)/(1/g**2), m, rel_tol=1e-14)
    # Sample arbitrary velocities, gradients and both charges; analytic proof is in DERIVATION.
    rng = random.Random(20261020)
    minimum_gap = math.inf
    square_error = 0.0
    for _ in range(12000):
        q1, q2, p1, p2, S = [rng.uniform(-3, 3) for _ in range(5)]
        F2 = q1*q1 + q2*q2
        potential = d['a']/4*(S*S-1)**2 + .5*S*S*F2 + d['b']/4*F2*F2
        q = q2*p1-q1*p2
        kinetic = .5*(p1*p1+p2*p2)
        positive_gradients = rng.random()*4
        minimum_gap = min(minimum_gap, kinetic+potential+positive_gradients-mu*abs(q))
        for sign in [-1, 1]:
            left = kinetic + .5*mu*mu*F2 - sign*mu*q
            right = .5*(p1-sign*mu*q2)**2 + .5*(p2+sign*mu*q1)**2
            square_error = max(square_error, abs(left-right)/(1+abs(left)+abs(right)))
    checks['sampled_local_energy_bound'] = minimum_gap >= -1e-12
    checks['square_completion_both_signs'] = square_error < 1e-14
    # Saturating homogeneous state is a negative control for an incorrectly doubled floor.
    F2 = math.sqrt(d['a']/d['b'])
    n = mu*F2
    E = .5*mu*mu*F2 + d['a']/4 + d['b']/4*F2*F2
    checks['sharp_coexistence_equality'] = abs(E-mu*n) < 1e-14
    checks['wrong_factor_two_detected'] = E-2*mu*n < -0.1
    regime_values = []
    for a, b in [(1, .1), (.1, .1), (1, 1), (2, 2), (.1, 40)]:
        svals = [0, .25, .9, 1, 1.2, 3]
        candidates = []
        for s in svals:
            y = max(1e-14, math.sqrt(a/b)*abs(s-1))
            candidates.append(s+a*(s-1)**2/(2*y)+b*y/2)
        expected = min(1, math.sqrt(a*b))
        err = abs(min(candidates)-expected)
        regime_values.append({'a': a, 'b': b, 'ratio_minimum_error': err})
    checks['potential_ratio_both_binding_regimes'] = max(x['ratio_minimum_error'] for x in regime_values) < 1e-12
    nominal = floor_values(protocol, cal)
    conservative = floor_values(protocol, cal, True)
    with localcontext() as context:
        context.prec = 60
        D = lambda x: Decimal(str(x))
        muD = min(Decimal(1), (D(d['a'])*D(d['b'])).sqrt().sqrt())
        floorD = D(cal['omega_b'])*muD*D(m)*D(abs(d['r_YQ0_over_YB0']))/D(cal['mean_mass_per_baryon_GeV_convention'])
        HD = Decimal(100)*(D(cal['omega_b'])+floorD).sqrt()
    decimal_error = abs(nominal['H_floor_km_s_Mpc']/float(HD)-1)
    checks['decimal60_floor_conversion'] = decimal_error < protocol['numerical_requirements']['float_decimal_relative_tolerance']
    old = read(BASE / 'GEOMETRIC_PREDICTION.json')
    old_errors = {
        'omega_floor': abs(nominal['omega_hidden_floor']/old['omega_hidden_floor']-1),
        'H0_floor': abs(nominal['H_floor_km_s_Mpc']/old['H0_lower_bound_km_s_Mpc']-1),
        'conservative_H0': abs(conservative['H_floor_km_s_Mpc']/old['H0_lower_bound_BBN_3sigma_mass_allowance']-1),
    }
    checks['old_exclusion_arithmetic_reproduced'] = max(old_errors.values()) < protocol['numerical_requirements']['old_floor_relative_tolerance']
    # SI and entropy route, independent of old prediction code.
    c, hbar, G = cal['c_SI'], cal['h_SI']/(2*math.pi), cal['G_SI']
    GeV_J = cal['eV_J']*1e9
    GeV_kg = GeV_J/c**2
    H100_s = 100e3/(1e6*cal['parsec_m'])
    rho_unit = 3*H100_s**2/(8*math.pi*G)/GeV_kg/1e6
    thermal_wavenumber = cal['kB_SI']*cal['T_gamma_K']/(hbar*c)
    entropy = 2*math.pi**2/45*(43/11)*thermal_wavenumber**3/1e6
    nB = cal['omega_b']*rho_unit/cal['mean_mass_per_baryon_GeV_convention']
    Y = nB/entropy
    via_entropy = mu*m*abs(d['r_YQ0_over_YB0'])*Y*entropy/rho_unit
    checks['SI_entropy_ratio_cancellation'] = math.isclose(via_entropy, nominal['omega_hidden_floor'], rel_tol=1e-14)
    checks['yield_unit_route_matches_baseline'] = math.isclose(Y, read(BASE / 'EMPIRICAL_PREDICTION.json')['charge_yield_YQ_equals_late_YB'], rel_tol=1e-13)
    quadrature_errors = []
    zmax = protocol['sensitivity']['z_grid_max']
    exact = 2*(-math.expm1(-.5*math.log1p(zmax)))
    for steps in [16, 32, 64, 128]:
        quadrature_errors.append({'steps': steps, 'relative_error': abs(simpson_integral(zmax, steps)/exact-1)})
    checks['distance_analytic_vs_quadrature'] = quadrature_errors[-1]['relative_error'] < protocol['numerical_requirements']['analytic_quadrature_relative_tolerance']
    decimal_D_errors = [abs(envelope(z, nominal['H_floor_km_s_Mpc'], c/1000)/decimal_envelope(z, nominal['H_floor_km_s_Mpc'], c/1000)-1)
                        for z in [1e-8, .001, .01, zmax]]
    checks['distance_float_vs_decimal60'] = max(decimal_D_errors) < protocol['numerical_requirements']['float_decimal_relative_tolerance']
    errors = []
    for z in [.002, .01, .04]:
        for vp in [-1500, -750, 0, 750, 1500]:
            doppler = math.sqrt((1+vp/(c/1000))/(1-vp/(c/1000)))
            observed = (1+z)*doppler-1
            errors.append(abs((1+observed)/doppler-1-z))
    checks['exact_Doppler_inverse'] = max(errors) < protocol['numerical_requirements']['Doppler_inverse_absolute_tolerance']
    checks['gradient_FRW_factors_positive'] = all(A**-2 > 0 for A in [1e-10, 1e-3, 1, 100])
    wb, wh = nominal['omega_b'], nominal['omega_hidden_floor']
    # Demonstrates why the nonnegative-Lambda hypothesis is load-bearing.
    compensation = -.8*wh
    checks['negative_Lambda_counterexample_to_unqualified_H0_bound'] = wb+wh+compensation < wb+wh
    obs = read(P.parent / 'grut_candidate_1_20261018/code/results.json')['observables']
    M_GeV = v/g*obs['Ebar']
    M_kg = M_GeV*GeV_kg
    R_m = obs['R90_bar']/m*(hbar/GeV_J)*c
    compactness = 2*G*M_kg/(R_m*c*c)
    optical_w = 4*G*M_kg*(cal['eV_J']/hbar)/c**3
    checks['reference_gravity_is_weak_only_for_reference_object'] = compactness < 1e-20
    checks['reference_optical_geometric_lensing_not_accessible'] = optical_w < 1e-20
    # No calibrated scale/yield are introduced by these numerical checks.
    for name, value in checks.items():
        require(value, 'Audit check failed: ' + name)
    return {
        'scope': 'Same-author analytic rederivation and implementation checks, not scientific replication or all-orders quantum verification.',
        'checks': checks, 'check_count': len(checks), 'all_checks_pass': all(checks.values()),
        'mu': mu, 'mstar_GeV': mu*m, 'nominal': nominal, 'conservative': conservative,
        'decimal60_H_floor': str(HD), 'decimal_floor_relative_error': decimal_error,
        'old_arithmetic_relative_errors': old_errors,
        'square_completion_max_normalized_residual': square_error,
        'sampled_energy_minimum_gap': minimum_gap, 'sample_count': 12000,
        'potential_ratio_regimes': regime_values,
        'distance_quadrature_convergence': quadrature_errors,
        'distance_max_decimal_relative_error': max(decimal_D_errors),
        'Doppler_max_inverse_absolute_error': max(errors),
        'units': {'rho_critical_over_h2_GeV_cm3': rho_unit, 'visible_entropy_cm3_convention': entropy,
                  'nB_cm3': nB, 'YB_entropy_convention': Y,
                  'entropy_route_omega_floor': via_entropy,
                  'entropy_note': 'Ideal standard photon+neutrino entropy convention for unit checking; cancels from primary bound.'},
        'reference_object_accessibility_only': {'M_GeV': M_GeV, 'M_kg': M_kg, 'R90_m': R_m,
                                               'Schwarzschild_radius_over_R90': compactness,
                                               'optical_1eV_wave_parameter': optical_w,
                                               'interpretation': 'One old prepared object, not a cosmological mass function. No detector event rate is inferred.'},
        'uncertainty_separation': {
            'BBN_fractional_sigma': cal['omega_b_sigma']/cal['omega_b'],
            'H_floor_BBN_linearized_sigma_km_s_Mpc': .5*nominal['H_floor_km_s_Mpc']*cal['omega_b_sigma']/cal['omega_b'],
            'baryon_mass_allowance_fraction': cal['mean_mass_systematic_fraction'],
            'hidden_fraction_of_minimum_total': wh/(wb+wh),
            'unquantified': ['EFT loop/running/higher-operator changes', 'full cosmological realization and BBN consistency', 'strong-binding/backreaction and collapse', 'raw maser posterior tails/covariances and alternative disk/systematic models'],
        },
        'decision': 'Old arithmetic and leading conditional bound survive. Correct the characterization of the published H0 inference. No independent scalar scale or hidden yield calibration exists in this packet.'}


def prediction():
    # Deliberately no test file, old comparison file, maser number or fitted H0 here.
    protocol, cal = inputs()
    d, mu = parameters(protocol)
    return {'protocol_id': protocol['id'], 'protocol_sha256': sha(P / 'PROTOCOL.json'),
            'calibration_sha256': sha(P / protocol['calibration_file']),
            'source_sha256': sha(P / 'audit.py'),
            'parameters': d, 'mu': mu, 'mstar_GeV': mu*d['mX_GeV'],
            'nominal': floor_values(protocol, cal),
            'conservative': floor_values(protocol, cal, True),
            'observable': protocol['observable'],
            'claim': 'Conditional upper distance envelope and necessary mass-yield constraints. No independently calibrated positive prediction or new data exposure.'}


def curve_text(pred):
    protocol, cal = inputs()
    stream = io.StringIO()
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(['z', 'D_A_max_nominal_Mpc', 'D_A_max_conservative_Mpc'])
    count = protocol['sensitivity']['z_grid_points']
    for i in range(count):
        z = protocol['sensitivity']['z_grid_max']*i/(count-1)
        writer.writerow([format(z, '.17g'), format(envelope(z, pred['nominal']['H_floor_km_s_Mpc'], cal['c_SI']/1000), '.17g'),
                         format(envelope(z, pred['conservative']['H_floor_km_s_Mpc'], cal['c_SI']/1000), '.17g')])
    return stream.getvalue()


def compare(pred):
    protocol, cal = inputs()
    data = read(P / protocol['test_file'])
    c = cal['c_SI']/1000
    nwidth = protocol['sensitivity']['width_multiplier']
    conservative = pred['conservative']
    rows = []
    for row in data['rows']:
        z = row['v_km_s']/c
        Dnom = envelope(z, pred['nominal']['H_floor_km_s_Mpc'], c)
        Dlow = row['D_Mpc']-nwidth*row['D_minus_Mpc']
        require(Dlow > 0, 'Nonpositive lower-distance corner; need a different declared estimator')
        corners = {}
        for label, sigma in [('primary', protocol['sensitivity']['peculiar_sigma_km_s']),
                             ('stress', protocol['sensitivity']['stress_peculiar_sigma_km_s'])]:
            beta = -nwidth*sigma/c
            observed_high_z = (row['v_km_s']+nwidth*row['v_sigma_km_s'])/c
            zhigh = (1+observed_high_z)/math.sqrt((1+beta)/(1-beta))-1
            require(0 < zhigh < .1, 'Envelope monotonic-corner domain exceeded')
            Dmax = envelope(zhigh, conservative['H_floor_km_s_Mpc'], c)
            Hmax = envelope(zhigh, 1.0, c)/Dlow
            product_limit = conservative['baryon_mass_GeV']/pred['mu']*((Hmax/100)**2/conservative['omega_b']-1)
            corners[label] = {'peculiar_sigma_km_s': sigma, 'cosmological_z_max': zhigh,
                              'D_A_max_Mpc': Dmax, 'observed_lower_over_upper_envelope': Dlow/Dmax,
                              'incompatible_with_declared_corners': Dlow>Dmax,
                              'necessary_mX_abs_YQ_over_YB_upper_GeV': product_limit}
        rows.append({'name': row['name'], 'z_observed_optical': z, 'D_median_Mpc': row['D_Mpc'],
                     'D_lower_sensitivity_Mpc': Dlow, 'D_A_max_nominal_Mpc': Dnom,
                     'nominal_D_ratio': row['D_Mpc']/Dnom, 'corners': corners})
    primary = [r['name'] for r in rows if r['corners']['primary']['incompatible_with_declared_corners']]
    stress = [r['name'] for r in rows if r['corners']['stress']['incompatible_with_declared_corners']]
    return {'dataset': data['dataset_id'], 'test_sha256': sha(P / protocol['test_file']),
            'prediction_sha256': sha(P / 'PREDICTION.json'), 'rows': rows,
            'primary_incompatible_rows': primary, 'stress_incompatible_rows': stress,
            'minimum_primary_necessary_product_limit_GeV': min(r['corners']['primary']['necessary_mX_abs_YQ_over_YB_upper_GeV'] for r in rows),
            'minimum_stress_necessary_product_limit_GeV': min(r['corners']['stress']['necessary_mX_abs_YQ_over_YB_upper_GeV'] for r in rows),
            'outcome': 'A specified preparation is excluded within the stated physical hypotheses and published-summary sensitivity treatment.' if primary else 'Necessary bound remains compatible within declared corners; no positive identification.',
            'global_model_outcome': 'The model cannot yet make a discriminating independently fixed empirical prediction because physical scale and primordial yield inputs are missing.',
            'statistical_limit': 'Retrospective summary-based interval sensitivity, not a full likelihood, calibrated p-value, blind test, or independent replication. The old H0 result uses the same masers and is not multiplied into this comparison.',
            'no_refit': 'mX=100 GeV and r=1 retained; product limits are not used to select a successful point.'}


def main():
    guard()
    command = sys.argv[1] if len(sys.argv)>1 else 'verify'
    if command == 'audit':
        baseline = baseline_integrity()
        results = numerical_audit()
        write(P / 'BASELINE_AUDIT.json', baseline)
        write(P / 'EXCLUSION_AUDIT.json', results)
        write(P / 'ENVIRONMENT.json', {'python': sys.version, 'platform': platform.platform(),
                                      'dependencies': 'Python standard library only',
                                      'decimal_precision': 60, 'random_seed': 20261020,
                                      'client_timezone': 'America/Chicago',
                                      'commands': ['python audit.py audit', 'python audit.py predict', 'python audit.py compare', 'python audit.py verify']})
        print(json.dumps({'baseline': baseline, 'checks': results['check_count'],
                          'nominal': results['nominal'], 'conservative': results['conservative']}))
    elif command == 'predict':
        pred = prediction()
        write(P / 'PREDICTION.json', pred)
        (P / 'DISTANCE_ENVELOPE.csv').write_text(curve_text(pred))
        print(json.dumps(pred))
    elif command == 'compare':
        pred = read(P / 'PREDICTION.json')
        require(pred == prediction(), 'Prediction differs from frozen inputs/source')
        result = compare(pred)
        write(P / 'COMPARISON.json', result)
        print(json.dumps(result, indent=2))
    elif command == 'verify':
        checks = {'baseline_packet_and_all_raw_hashes': baseline_integrity() == read(P / 'BASELINE_AUDIT.json'),
                  'audit_recomputed': numerical_audit() == read(P / 'EXCLUSION_AUDIT.json'),
                  'prediction_recomputed': prediction() == read(P / 'PREDICTION.json'),
                  'comparison_recomputed': compare(read(P / 'PREDICTION.json')) == read(P / 'COMPARISON.json'),
                  'distance_curve_recomputed': curve_text(read(P / 'PREDICTION.json')) == (P / 'DISTANCE_ENVELOPE.csv').read_text()}
        for name, value in checks.items():
            require(value, 'Verification failed: '+name)
        write(P / 'VERIFICATION.json', {'checks': checks, 'all_checks_pass': all(checks.values()),
                                        'scope': 'Same-author source/integrity and alternate numerical checks; no independent scientific replication or evidence of fields existing.'})
        print(json.dumps(checks))
    else:
        raise SystemExit('audit | predict | compare | verify')


if __name__ == '__main__':
    main()
