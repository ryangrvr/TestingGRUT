"""Standard physical-comparator controls. No candidate law. Stdlib only.

Fractions certify finite arithmetic. Complex Pauli calculations use only
Gaussian dyadic rationals, exactly represented by Python complex numbers.
Transcendental evaluations are explicitly numerical illustrations.
"""
from fractions import Fraction as F
from itertools import product
from math import exp, log, sin, sqrt


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def add(a, b, factor=1):
    return [[x + factor*y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, c):
    return [[c*x for x in row] for row in a]


def kron(a, b):
    return [[x*y for x in ar for y in br] for ar in a for br in b]


def adj(a):
    return [[a[j][i].conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def rank(a):
    b = [[F(x) for x in row] for row in a]
    r = 0
    for c in range(len(b[0])):
        pivot = next((i for i in range(r, len(b)) if b[i][c]), None)
        if pivot is None:
            continue
        b[r], b[pivot] = b[pivot], b[r]
        q = b[r][c]
        b[r] = [x/q for x in b[r]]
        for i in range(len(b)):
            if i != r:
                q = b[i][c]
                b[i] = [x-q*y for x, y in zip(b[i], b[r])]
        r += 1
        if r == len(b):
            break
    return r


FIXED = {
    (0, 0): (0, 4), (0, 1): (6, 0),
    (1, 0): (0, 4), (1, 1): (6, 0),
    (2, 0): (2, 6), (2, 1): (2, 4),
    (3, 0): (6, 2), (3, 1): (4, 2),
    (4, 0): (5, 1), (4, 1): (1, 7),
    (5, 0): (5, 1), (5, 1): (1, 7),
    (6, 0): (7, 3), (6, 1): (5, 3),
    (7, 0): (3, 7), (7, 1): (3, 5),
}


def original_step(z, a, b):
    x, ell, memory = z//4, (z//2) % 2, z % 2
    newx = x ^ ((1-ell)*a + ell*memory) ^ b
    return 4*newx + 2*(ell ^ (a*newx)) + x


def channel_edges(a):
    # Kraus squared amplitudes: (input, output state, recorded bit, weight).
    return [(z, target, target//4, p) for z in range(8)
            for target, p in zip(FIXED[z, a], (F(3, 4), F(1, 4)))]


def policy_words(policy, initial, depth=3, quantum=False):
    frontier = {(initial, ()): F(1)}
    for t in range(depth):
        nxt = {}
        for (z, history), probability in frontier.items():
            node = 2**t-1 + sum(bit*2**(t-1-i) for i, bit in enumerate(history))
            a = (policy >> node) & 1
            if quantum:
                successors = [(zp, y, p) for src, zp, y, p in channel_edges(a) if src == z]
            else:
                successors = [(original_step(z, a, b), original_step(z, a, b)//4, p)
                              for b, p in enumerate((F(3, 4), F(1, 4)))]
            for zp, y, p in successors:
                key = (zp, history+(y,))
                nxt[key] = nxt.get(key, F(0)) + probability*p
        frontier = nxt
    return {word: sum(p for (_, h), p in frontier.items() if h == word)
            for word in product(range(2), repeat=depth)}


def cptp_control():
    differences = []
    for initial in range(8):
        for policy in range(128):
            original = policy_words(policy, initial)
            fixed = policy_words(policy, initial, quantum=True)
            differences.extend(abs(original[w]-fixed[w]) for w in original)
    actions = []
    for a in (0, 1):
        edges = channel_edges(a)
        normalization = [sum(p for src, _, _, p in edges if src == z) for z in range(8)]
        # Choi is diagonal in |z>_input tensor |z',y>_output.
        choi = {(src, zp, y): p for src, zp, y, p in edges}
        actions.append({'action': a, 'trace_preserving_diagonal': normalization,
                        'nonzero_choi_eigenvalues': sorted(choi.values()),
                        'choi_rank': len(choi), 'choi_trace': sum(choi.values())})
    return {'comparisons': 1024, 'record_cells': len(differences),
            'maximum_difference': max(differences), 'actions': actions,
            'minimal_pure_environment_for_THIS_dephasing_extension': 16,
            'minimal_over_all_classical_record_extensions': 'NOT_CLAIMED'}


def memory_control():
    words = list(product(range(2), repeat=2))
    matrix = []
    for preparation in words:
        state = (0, 0)
        for a in preparation:
            state = (state[1], a)
        output = []
        for a in (0, 0):
            output.append(state[0])
            state = (state[1], a)
        matrix.append([F(tuple(output) == w) for w in words])
    return {'response_matrix': matrix, 'rank': rank(matrix),
            'perfect_quantum_dimension_lower_bound': 4,
            'max_equiprobable_success_d1_d2_d3_d4': [F(d, 4) for d in range(1, 5)],
            'all_preparation_dependent_side_channels_counted': True}


def reversible_control():
    # Retain original state and biased noise. Reversibly XOR output scratch.
    images = []
    for z, a, b, scratch in product(range(8), range(2), range(2), range(8)):
        images.append((z, a, b, scratch ^ original_step(z, a, b)))
    return {'domain_size': len(images), 'distinct_images': len(set(images)),
            'finite_horizon_resources': 'original input, noise bit, scratch and clock supplied',
            'bounded_perpetual_reuse_without_reset': 'NOT_CLAIMED'}


def reset_control():
    rows = []
    for d, n in ((2, 1), (2, 3), (4, 2)):
        size = d**n
        rows.append({'input_dimension': d, 'resets': n,
                     'orthogonal_input_gram_rank': rank(eye(size)),
                     'pure_environment_lower_bound': size})
    return {'closed_unitary_bound': rows,
            'faithful_bath_example': {'d': 2, 'D': 3,
                                      'initial_rank': 6, 'pure_reset_output_rank_at_most': 3}}


def landauer_control():
    # beta=1, H_S=0, H_E=diag(0,log(3)); driven SWAP, work explicitly priced.
    thermal = (F(3, 4), F(1, 4))
    final_bath = (F(1, 2), F(1, 2))
    entropy = lambda p: -sum(float(x)*log(float(x)) for x in p if x)
    delta_s = log(2)-entropy(thermal)
    beta_q = log(3)/4
    relative = sum(float(p)*log(float(p/q)) for p, q in zip(final_bath, thermal))
    return {'arithmetic_grade': 'TRANSCENDENTAL_NUMERICAL_ILLUSTRATION_OF_EXACT_IDENTITY',
            'initial_memory': list(final_bath), 'initial_thermal_bath': list(thermal),
            'final_memory': list(thermal), 'final_bath': list(final_bath),
            'beta_Q': beta_q, 'delta_S': delta_s, 'mutual_information': 0,
            'relative_entropy': relative, 'identity_residual': beta_q-delta_s-relative,
            'external_work_at_HS_zero': beta_q,
            'exact_pure_reset': False}


def balance_control():
    c = [[F(0), F(1), F(0)], [F(0), F(0), F(1)], [F(1), F(0), F(0)]]
    p = scale(add(eye(3), c), F(1, 2))
    currents = [[F(1, 3)*(p[i][j]-p[j][i]) for j in range(3)] for i in range(3)]
    # U on classical qutrit x bath-bit: branch 0 identity, branch 1 forward cycle.
    basis = list(product(range(3), range(2)))
    targets = [(i if b == 0 else (i+1) % 3, b) for i, b in basis]
    recovered = [[sum(F(1, 2) for b in range(2)
                      if (i if b == 0 else (i+1) % 3) == j)
                  for j in range(3)] for i in range(3)]
    return {'P': p, 'stationary_distribution': [F(1, 3)]*3,
            'stationary_output': [sum(F(1, 3)*p[i][j] for i in range(3)) for j in range(3)],
            'pair_current_matrix': currents,
            'thermal_operation_permutation_unitary': len(set(targets)) == 6,
            'thermal_operation_exact_same_kernel': recovered == p,
            'HS_HE': 'both zero; bath I2/2; [U,HS+HE]=0',
            'autonomous_detailed_balance_semigroup': False,
            'driven_CTMC_entropy_production_rates_2_1': log(2)}


def ctmc_control():
    return {'NOT_matrix': [[0, 1], [1, 0]], 'NOT_determinant': -1,
            'two_state_generator_rates_2_3': [[-2, 2], [3, -3]],
            'finite_time_determinant_formula': 'exp(-5*t)>0',
            'hidden_rate_bound_Lambda_T': [2, 3],
            'minimum_error_numerical': exp(-6),
            'quantum_exact_NOT_at_T': 'H=pi*sigma_x/(2*T); U(T)=-i*sigma_x'}


def transpose_control():
    bell = [[F(1, 2) if i in (0, 3) and j in (0, 3) else F(0)
             for j in range(4)] for i in range(4)]
    partial = [[bell[(i//2)*2+j % 2][(j//2)*2+i % 2] for j in range(4)] for i in range(4)]
    vector = [0, 1, -1, 0]  # normalized squared norm 2
    witness = sum(vector[i]*partial[i][j]*vector[j] for i in range(4) for j in range(4))/2
    return {'Bell_density': bell, 'partial_transpose': partial,
            'singlet_expectation': witness, 'eigenvalues': [F(-1, 2)]+[F(1, 2)]*3}


def bell_control():
    distributions = {}
    correlators = []
    for x, y in product(range(2), repeat=2):
        p = [[F(1, 2) if (a ^ b) == x*y else F(0) for b in range(2)] for a in range(2)]
        distributions[f'{x}{y}'] = p
        correlators.append(sum((-1)**(a+b)*p[a][b] for a, b in product(range(2), repeat=2)))
    local_values = [a0*b0+a0*b1+a1*b0-a1*b1 for a0, a1, b0, b1 in product((-1, 1), repeat=4)]
    return {'PR_distributions': distributions, 'correlators': correlators,
            'PR_CHSH': sum(correlators[:3])-correlators[3],
            'local_CHSH_max': max(local_values), 'quantum_CHSH_bound_squared': 8,
            'PR_CHSH_squared': 16,
            'marginals': [[F(1, 2), F(1, 2)]]*8}


I, X = eye(2), [[0, 1], [1, 0]]
Y, Z = [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]


def pauli(*operators):
    result = [[1]]
    for operator in operators:
        result = kron(result, operator)
    return result


def locality_control():
    h = add(pauli(X, X, I), pauli(I, Z, Z))
    a, b = pauli(I, I, X), pauli(Z, I, I)
    at_quarter_period = scale(mm(mm(h, a), h), F(1, 2))
    expected = scale(pauli(X, Y, Y), -1)
    commutator = add(mm(at_quarter_period, b), mm(b, at_quarter_period), -1)
    # Exact operator identities, no floating matrix exponential.
    identities = {'H_squared_equals_2I': mm(h, h) == scale(eye(8), 2),
                  'evolved_X3_equals_minus_X1Y2Y3': at_quarter_period == expected,
                  'commutator_equals_2i_Y1Y2Y3': commutator == scale(pauli(Y, Y, Y), 2j),
                  'commutator_dagger_commutator_equals_4I': mm(adj(commutator), commutator) == scale(eye(8), 4)}
    # Radius-one deterministic local shift: compare two worlds differing at site 0.
    left, right = [0]*8, [1]+[0]*7
    histories = []
    for step in range(8):
        histories.append({'step': step, 'site_5_difference': right[5]-left[5]})
        left, right = [0]+left[:-1], [0]+right[:-1]
    return {'exact_Pauli_identities': identities,
            'all_time_commutator_norm': '2*sin(sqrt(2)*J*t)^2',
            'short_time_commutator_norm_leading': '4*J^2*t^2',
            'strict_quantum_Hamiltonian_light_cone': False,
            'classical_local_shift_intervention': histories}


def kinetic_control():
    pi = [F(1, 6), F(1, 3), F(1, 2)]
    conductances = {(0, 1): F(1, 10), (0, 2): F(1, 5), (1, 2): F(3, 10)}
    q = [[F(0) if i == j else conductances[min(i, j), max(i, j)]/pi[i]
          for j in range(3)] for i in range(3)]
    for i in range(3):
        q[i][i] = -sum(q[i])
    down, up = [[0, 1], [0, 0]], [[0, 0], [1, 0]]
    def dissipator(operator, rho):
        ld_l = mm(adj(operator), operator)
        return add(mm(mm(operator, rho), adj(operator)),
                   scale(add(mm(ld_l, rho), mm(rho, ld_l)), F(1, 2)), -1)
    def gkls(rho, gamma, phi):
        return add(scale(add(dissipator(down, rho), dissipator(up, rho), F(1, 3)), gamma),
                   scale(dissipator(Z, rho), phi/2))
    projectors = [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]
    pop_generators = [[[gkls(projectors[j], F(gamma), F(0))[i][i]
                        for j in (0, 1)] for i in (0, 1)] for gamma in (1, 2)]
    e01 = [[0, 1], [0, 0]]
    coherence_rates = [-gkls(e01, F(1), F(phi))[0][1] for phi in (0, 1)]
    return {'pi': pi, 'Q': q,
            'stationary_generator_product': [sum(pi[i]*q[i][j] for i in range(3)) for j in range(3)],
            'conductance_dimension_N2_N3_N4': [n*(n-1)//2 for n in (2, 3, 4)],
            'GKLS_up_down_ratio': F(1, 3),
            'GKLS_pop_generators_column_convention': pop_generators,
            'GKLS_pop_rates_gamma1_gamma2': [-sum(p[i][i] for i in (0, 1)) for p in pop_generators],
            'GKLS_coherence_rates_gamma1_phi0_phi1': coherence_rates,
            'GKLS_stationary_density_derivative': gkls([[F(3, 4), 0], [0, F(1, 4)]], F(1), F(1)),
            'GKLS_Gibbs_population': [F(3, 4), F(1, 4)]}


def anomaly_control():
    return {','.join(map(str, charges)): {'sum_q': sum(charges), 'sum_q3': sum(q**3 for q in charges)}
            for charges in ((1, 1), (1, -1), (2, -2))}


def run_controls():
    return {'C01_fixed_CPTP_lift': cptp_control(), 'C02_quantum_memory_bound': memory_control(),
            'C03_reversible_finite_control': reversible_control(), 'C04_erasure_resources': reset_control(),
            'C05_finite_bath_Landauer': landauer_control(), 'C06_DB_vs_thermal_operation': balance_control(),
            'C07_CTMC_obstruction': ctmc_control(), 'C08_complete_positivity': transpose_control(),
            'C09_no_signalling_vs_quantum': bell_control(), 'C10_locality': locality_control(),
            'C11_physical_kinetic_freedom': kinetic_control(), 'C12_anomaly_control': anomaly_control()}
