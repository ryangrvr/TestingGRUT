"""Exact controls for published consistency mechanisms, not a GRUT law search."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json
import math
from pathlib import Path


def encoded(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {str(k): encoded(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encoded(v) for v in x]
    return x


def determinant(matrix):
    n = len(matrix)
    if n == 0:
        return F(1)
    if n == 1:
        return matrix[0][0]
    return sum((-1)**j * matrix[0][j] * determinant([
        row[:j] + row[j+1:] for row in matrix[1:]
    ]) for j in range(n))


def principal_minors(matrix):
    return [determinant([[matrix[i][j] for j in ids] for i in ids])
            for size in range(1, len(matrix)+1)
            for ids in combinations(range(len(matrix)), size)]


def rank(matrix):
    a = [[F(x) for x in row] for row in matrix]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][c]
        a[r] = [x/scale for x in a[r]]
        for i in range(len(a)):
            if i != r:
                scale = a[i][c]
                a[i] = [x-scale*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def alternating_tensor(n, seeds, scale=F(1)):
    values = {}
    for seed, value in seeds.items():
        for p in permutations(range(3)):
            inversions = sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
            values[tuple(seed[i] for i in p)] = value*scale*(-1)**inversions
    return lambda a, b, c: values.get((a, b, c), F(0))


def jacobi_cells(n, f):
    return [sum(f(b, c, e)*f(a, e, d) + f(c, a, e)*f(b, e, d)
                + f(a, b, e)*f(c, e, d) for e in range(n))
            for a, b, c, d in product(range(n), repeat=4)]


def run():
    # A+B -> A+B, on-shell massive matter, signature (+---).
    momenta = [[5, 0, 0, 3], [5, 0, 0, -3], [5, 3, 0, 0], [5, -3, 0, 0]]
    signs = [-1, -1, 1, 1]
    species = [0, 1, 0, 1]
    masses2 = [p[0]**2-sum(x*x for x in p[1:]) for p in momenta]
    assert masses2 == [16]*4
    assert all(sum(s*p[k] for s, p in zip(signs, momenta)) == 0 for k in range(4))
    gravity = [[sum(s*p[k] for s, p, label in zip(signs, momenta, species)
                    if label == j) for j in range(2)] for k in range(4)]
    assert rank(gravity) == 1
    assert all(sum(row) == 0 for row in gravity)
    assert any(row[0] != 0 for row in gravity)
    photon = []
    for charges in [(1, 2), (-1, 3), (0, 5), (7, -2)]:
        residual = sum(s*charges[label] for s, label in zip(signs, species))
        assert residual == 0
        photon.append({'charges': charges, 'ward_residual': residual})
    universal = []
    for coupling in [F(0), F(1, 3), F(2)]:
        residual = [sum(coupling*x for x in row) for row in gravity]
        assert all(x == 0 for x in residual)
        universal.append({'coupling': coupling, 'ward_residual': residual})

    lie = []
    for scale in [F(0), F(1, 3), F(2)]:
        cells = jacobi_cells(3, alternating_tensor(3, {(0, 1, 2): F(1)}, scale))
        assert all(x == 0 for x in cells)
        lie.append({'scale': scale, 'cells': len(cells), 'nonzero_residuals': 0})
    bad = jacobi_cells(5, alternating_tensor(5, {(0, 1, 2): F(1), (0, 3, 4): F(1)}))
    assert any(x != 0 for x in bad)

    # Explicit nonnegative two-point variables. This is not a wormhole model.
    distributions = []
    for mu, variance in product([F(1, 2), F(1), F(3)], [F(0), F(1, 4), F(1), F(4), F(100)]):
        height = mu+variance/mu
        probability = mu/height
        assert 0 < probability <= 1
        moments = [F(1)] + [probability*height**n for n in range(1, 6)]
        assert moments[1] == mu and moments[2]-mu**2 == variance
        k3 = moments[3]-3*mu*moments[2]+2*mu**3
        assert k3 == variance**2/mu-mu*variance
        shifted_determinant = mu*moments[3]-moments[2]**2
        assert shifted_determinant == mu*k3+mu**2*variance-variance**2 == 0
        h0 = [[moments[i+j] for j in range(3)] for i in range(3)]
        h1 = [[moments[i+j+1] for j in range(3)] for i in range(3)]
        minors = principal_minors(h0)+principal_minors(h1)
        assert all(x >= 0 for x in minors)
        distributions.append({'mean': mu, 'variance': variance,
                              'positive_height': height, 'height_probability': probability,
                              'third_cumulant': k3, 'moments_0_through_5': moments,
                              'principal_minors': minors})

    # Gaussian moment proposals are not presumed to have nonnegative support.
    gaussian = []
    for variance in [F(1, 4), F(1), F(4), F(100)]:
        mu, k3 = F(1), F(0)
        m2, m3 = mu**2+variance, mu**3+3*mu*variance+k3
        d = mu*m3-m2**2
        assert d == variance*(1-variance)
        gaussian.append({'variance': variance, 'shifted_determinant': d,
                         'nonnegative_support_moment_test_passes': d >= 0})
    assert gaussian[2]['shifted_determinant'] < 0

    strictly_positive = []
    for variance in [F(1, 4), F(1), F(4), F(100)]:
        mu, low = F(1), F(1, 2)
        high = mu+variance/(mu-low)
        p = (mu-low)/(high-low)
        moments = [(1-p)*low**n+p*high**n for n in range(6)]
        assert moments[0] == 1 and moments[1] == mu
        assert moments[2]-mu**2 == variance
        k3 = moments[3]-3*mu*moments[2]+2*mu**3
        assert k3 == variance*(variance/(mu-low)-(mu-low))
        h0 = [[moments[i+j] for j in range(3)] for i in range(3)]
        h1 = [[moments[i+j+1] for j in range(3)] for i in range(3)]
        minors = principal_minors(h0)+principal_minors(h1)
        assert all(x >= 0 for x in minors)
        strictly_positive.append({'mean': mu, 'variance': variance,
                                  'low': low, 'high': high, 'high_probability': p,
                                  'third_cumulant': k3, 'principal_minors': minors})

    # 1-loop determinant ignored; exact rational proxy for the flux kernel.
    flux = []
    for attenuation, weight in [(F(1, 4), F(1)), (F(1, 4), F(2)), (F(1, 4), F(3)), (F(1, 4), F(4)), (F(1, 4), F(5))]:
        outward, inward = attenuation*weight, attenuation/weight
        converges = outward < 1 and inward < 1
        closed = (1-attenuation**2)/((1-outward)*(1-inward)) if converges else None
        assert converges == (F(1, 4) < weight < F(4))
        if converges:
            assert closed == 1+outward/(1-outward)+inward/(1-inward)
        flux.append({'exp_minus_C': attenuation, 'exp_x': weight,
                     'both_geometric_tails_converge': converges, 'closed_sum': closed})

    return encoded({
        'status': 'EXACT_INTERACTION_CONSISTENCY_CONTROLS_PASS',
        'scope': 'Known Ward/Jacobi/moment/geometric-series controls; no GPI evaluation or GRUT candidate',
        'soft_control': {'momenta': momenta, 'signs': signs, 'species': species,
                         'mass_squared': masses2, 'gravity_ward_matrix': gravity,
                         'gravity_rank': 1, 'gravity_null_direction': [1, 1],
                         'photon_fixtures': photon, 'universal_coupling_fixtures': universal},
        'jacobi_control': {'su2_fixtures': lie, 'antisymmetric_non_lie_cells': len(bad),
                           'antisymmetric_non_lie_nonzero_residuals': sum(x != 0 for x in bad)},
        'positive_moment_control': {'two_point_fixtures': distributions,
                                    'strictly_positive_fixtures': strictly_positive,
                                    'principal_minors_checked': (len(distributions)+len(strictly_positive))*14,
                                    'gaussian_moment_proposals': gaussian},
        'geometric_flux_proxy': flux,
        'four_dimensional_unit_charge_bound_coefficients': {
            'ads_fS_over_reduced_Mpl': math.pi/2*math.sqrt(2/3),
            'flat_fS_over_reduced_Mpl': math.pi/2*math.sqrt(3/2),
            'kind': 'Floating evaluation of published coefficients, not measured values or interval certificates'},
        'new_grut_law': False, 'card_2_opened': False, 'cards_remaining': 2,
        'independent_review': 'PENDING', 'experiment_performed': False
    })


def render_tables(result):
    moment = result['positive_moment_control']
    lines = ['# Generated exact controls', '', 'Generated from `results.json` by `controls.py`.', '',
             '| Control | Disclosed coverage |', '|---|---|',
             f"| Photon charge Ward identity | {len(result['soft_control']['photon_fixtures'])} charge assignments |",
             f"| Universal gravity Ward identity | rank {result['soft_control']['gravity_rank']}; {len(result['soft_control']['universal_coupling_fixtures'])} common-coupling assignments |",
             f"| SU(2) Jacobi identity | {sum(x['cells'] for x in result['jacobi_control']['su2_fixtures'])} exact cells |",
             f"| Antisymmetry without Jacobi | {result['jacobi_control']['antisymmetric_non_lie_nonzero_residuals']} nonzero cells of {result['jacobi_control']['antisymmetric_non_lie_cells']} (expected rejection) |",
             f"| Nonnegative and strictly positive distributions | {len(moment['two_point_fixtures'])+len(moment['strictly_positive_fixtures'])} fixtures; {moment['principal_minors_checked']} principal minors |",
             f"| Gaussian moment proposals | {len(moment['gaussian_moment_proposals'])} variances, including expected rejections |",
             f"| Geometric flux proxy | {len(result['geometric_flux_proxy'])} interior/boundary/exterior fixtures |", '',
             '| Fixed-mean hostile control | Mean | Variance | Support | High-value probability | Third cumulant |',
             '|---|---|---|---|---|---|']
    for x in moment['two_point_fixtures']:
        if x['mean'] == '1' and x['variance'] == '100':
            lines.append(f"| Nonnegative | {x['mean']} | {x['variance']} | 0, {x['positive_height']} | {x['height_probability']} | {x['third_cumulant']} |")
    for x in moment['strictly_positive_fixtures']:
        if x['variance'] == '100':
            lines.append(f"| Strictly positive | {x['mean']} | {x['variance']} | {x['low']}, {x['high']} | {x['high_probability']} | {x['third_cumulant']} |")
    lines += ['', 'These are exact mathematical controls, not quantum-gravity realizations or experimental data.', '']
    return '\n'.join(lines)


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('results.json').write_text(json.dumps(result, indent=2)+'\n')
    Path(__file__).with_name('CONTROL_TABLES.md').write_text(render_tables(result))
    print(result['status'])
    moment = result['positive_moment_control']
    print(f"{len(moment['two_point_fixtures'])+len(moment['strictly_positive_fixtures'])} positive-distribution fixtures; {moment['principal_minors_checked']} principal minors; coverage in CONTROL_TABLES.md")
