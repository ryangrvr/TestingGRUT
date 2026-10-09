"""Issue #5: exact conventional controls, not a GRUT law or law search."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial


# Independent fixed transition table: (state id, intervention) -> successors
# under noise b=0 and b=1, with probabilities 3/4 and 1/4 respectively.
# id = 4*x + 2*lambda + memory. No equation evaluator is used by this table.
TABLE = {
    (0, 0): (0, 4), (0, 1): (6, 0),
    (1, 0): (0, 4), (1, 1): (6, 0),
    (2, 0): (2, 6), (2, 1): (2, 4),
    (3, 0): (6, 2), (3, 1): (4, 2),
    (4, 0): (5, 1), (4, 1): (1, 7),
    (5, 0): (5, 1), (5, 1): (1, 7),
    (6, 0): (7, 3), (6, 1): (5, 3),
    (7, 0): (3, 7), (7, 1): (3, 5),
}


def adaptive_step(state, action):
    x, law, memory = state
    for noise, weight in ((0, F(3, 4)), (1, F(1, 4))):
        xp = x ^ ((1-law)*action + law*memory) ^ noise
        lp = law ^ (action*xp)
        yield (xp, lp, x), xp, weight


def fixed_step(state, action, table=TABLE):
    for successor, weight in zip(table[state, action], (F(3, 4), F(1, 4))):
        yield successor, successor//4, weight


def records(initial, horizon, policy, step):
    """Exact marginal law of the (action, observation) record, no simulation."""
    paths = {(initial, ()): F(1)}
    for _ in range(horizon):
        following = defaultdict(F)
        for (state, history), mass in paths.items():
            pa = policy(history)
            for action, paction in ((0, 1-pa), (1, pa)):
                if not paction:
                    continue
                for successor, observation, p in step(state, action):
                    following[successor, history+((action, observation),)] += mass*paction*p
        paths = following
    marginal = defaultdict(F)
    for (_, history), mass in paths.items():
        marginal[history] += mass
    return dict(marginal)


def deterministic_policy(bits):
    def policy(history):
        observations = [y for _, y in history]
        node = (1 << len(observations))-1
        index = sum(y << (len(observations)-i-1) for i, y in enumerate(observations))
        return F(bits[node+index])
    return policy


def equivalence_certificate():
    cells = comparisons = 0
    maximum = F(0)
    for bits in product((0, 1), repeat=7):
        policy = deterministic_policy(bits)
        for initial in product((0, 1), repeat=3):
            encoded = 4*initial[0]+2*initial[1]+initial[2]
            left = records(initial, 3, policy, adaptive_step)
            right = records(encoded, 3, policy, fixed_step)
            assert sum(left.values()) == sum(right.values()) == 1
            comparisons += 1
            for h in left.keys() | right.keys():
                cells += 1
                maximum = max(maximum, abs(left.get(h, 0)-right.get(h, 0)))
    return {"deterministic_policies": 128, "initial_states": 8,
            "horizon": 3, "comparisons": comparisons,
            "record_cells": cells, "maximum_exact_difference": maximum}


def rank(matrix):
    a = [[F(v) for v in row] for row in matrix]
    lead = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(lead, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[lead], a[pivot] = a[pivot], a[lead]
        divisor = a[lead][column]
        a[lead] = [v/divisor for v in a[lead]]
        for i in range(len(a)):
            if i != lead:
                scale = a[i][column]
                a[i] = [x-scale*y for x, y in zip(a[i], a[lead])]
        lead += 1
        if lead == len(a):
            break
    return lead


def shift_register(action_word, initial=(0, 0)):
    state = initial
    observed = []
    for action in action_word:
        observed.append(state[0])
        state = (state[1], action)
    return state, tuple(observed)


def response_matrix():
    words = list(product((0, 1), repeat=2))
    matrix = []
    for preparation in words:
        state, _ = shift_register(preparation)
        _, response = shift_register((0, 0), state)
        matrix.append([int(response == word) for word in words])
    determinant = 0
    for permutation in permutations(range(4)):
        inversions = sum(permutation[i] > permutation[j] for i in range(4) for j in range(i+1, 4))
        term = (-1)**inversions
        for i, j in enumerate(permutation):
            term *= matrix[i][j]
        determinant += term
    return {"preparations": words, "probe_actions": (0, 0), "responses": words,
            "matrix": matrix, "rank": rank(matrix), "determinant": determinant,
            "excluded_classical_hidden_state_bounds": [1, 2, 3],
            "fixed_comparator_states": 4}


def matadd(a, b, scale=1):
    return tuple(tuple(x+scale*y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def matmul(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(a)))
                       for j in range(len(a))) for i in range(len(a)))


def smolin_step(previous, current):
    commutator = matadd(matmul(current, previous), matmul(previous, current), -1)
    return matadd(matadd(current, current), matadd(previous, commutator, -1), -1)


def matrix_fixture():
    x0 = ((0, 1, 0), (-1, 0, 0), (0, 0, 0))
    x1 = ((0, 1, 1), (-1, 0, 0), (-1, 0, 0))
    matrices = [x0, x1]
    pair = (x0, x1)
    for _ in range(4):
        new = smolin_step(matrices[-2], matrices[-1])
        # The autonomous state is a pair of matrices; this is the full law.
        pair = (pair[1], smolin_step(*pair))
        matrices.append(new)
        assert pair == tuple(matrices[-2:])
    return {"N": 3, "matrices": matrices,
            "autonomous_pair_exact_match": True,
            "antisymmetric": all(all(a[i][j] == -a[j][i] for i in range(3)
                                      for j in range(3)) for a in matrices),
            "traces": [sum(a[i][i] for i in range(3)) for a in matrices]}


def contraction(target, initial, steps):
    x = F(initial)
    for _ in range(steps):
        x = (x+target)/2
    return x


def bistable(initial, steps):
    x = F(initial)
    for _ in range(steps):
        x += x*(1-x*x)/4
    return x


def selection(weight, p, steps):
    p = F(p)
    for _ in range(steps):
        p = weight*p/(weight*p+1-p)
    return p


def urn_law(word, alpha=1, beta=1):
    success, failure = alpha, beta
    probability = F(1)
    for y in word:
        probability *= F(success if y else failure, success+failure)
        success += y
        failure += 1-y
    return probability


def beta_mixture_uniform(word):
    """Integral_0^1 p^k(1-p)^(n-k) dp, by the exact beta identity."""
    n, k = len(word), sum(word)
    return F(factorial(k)*factorial(n-k), factorial(n+1))


def urn_fixture():
    cells = 0
    maximum = F(0)
    joint = []
    for n in range(1, 7):
        total = F(0)
        for word in product((0, 1), repeat=n):
            left, right = urn_law(word), beta_mixture_uniform(word)
            total += left
            maximum = max(maximum, abs(left-right))
            cells += 1
            if n == 3:
                joint.append({"word": word, "urn_probability": left,
                              "static_beta_mixture_probability": right})
        assert total == 1
    return {"enumerated_word_cells": cells, "maximum_exact_difference": maximum,
            "three_trial_joint": joint, "pair_covariance": F(1, 12),
            "after_do_first_success": {"urn_next_success": F(2, 3),
                                      "static_mixture_next_success": F(1, 2)},
            "after_observe_first_success": {"both_next_success": F(2, 3)}}


def run_controls():
    targets = (F(1, 137), F(2, 137))
    initials = (F(-1), F(0), F(1), F(10))
    contractions = [{"target": a, "initial": x, "steps": 6,
                     "value": contraction(a, x, 6),
                     "error": abs(contraction(a, x, 6)-a),
                     "closed_form_error": abs(x-a)/64}
                    for a in targets for x in initials]
    p0 = F(1, 3)
    replication = [{"relative_fitness": w, "initial": p, "steps": 5,
                    "value": selection(F(w), p, 5),
                    "closed_form": F(w**5)*p/(1-p+F(w**5)*p)}
                   for w in (1, 2) for p in (F(0), p0, F(1))]
    m = 0
    for action in [1, 0, 1, 1]:
        m = 2*m+action
    q, r = F(1, 4), F(1, 2)
    stationary = q/(q+r)
    p_one = lambda memory: sum(w for _, y, w in adaptive_step((0, 1, memory), 0) if y)
    return {
        "C01_kernel": {"states": len({s for s, _ in TABLE}), "interventions": 2, "rows": len(TABLE),
                         "row_sums": [sum(p for _, _, p in fixed_step(s, a)) for s, a in TABLE],
                         "noise_weights": [F(3, 4), F(1, 4)]},
        "C02_all_feedback_policies": equivalence_certificate(),
        "C03_hidden_memory_witness": {"visible_state": 0, "intervention": 0,
            "law_label": 1, "P_next_one_memory_zero": p_one(0),
            "P_next_one_memory_one": p_one(1)},
        "C04_rank_bound": response_matrix(),
        "C05_state_law_unification": matrix_fixture(),
        "C06_unbounded_clock": {"definition": "y_n=1 iff n is a power of two, n>=1",
            "one_times_through_64": [1 << k for k in range(7)],
            "successive_gaps": [1 << k for k in range(6)],
            "claim_grade": "finite prefix only; nonperiodicity is proved in PROOFS.md"},
        "C07_growing_memory": {"actions": [1, 0, 1, 1],
            "binary_prefix_integer": m, "update": "m'=2*m+a",
            "state_domain": "nonnegative integers; finite coordinate, unbounded capacity"},
        "C08_contraction_selection": contractions,
        "C09_basin_dependence": {"fixed_points": [-1, 0, 1],
            "negative_initial_value_after_3": bistable(F(-1, 2), 3),
            "positive_initial_value_after_3": bistable(F(1, 2), 3),
            "zero_initial_value_after_3": bistable(0, 3)},
        "C10_replication": replication,
        "C11_stochastic_constant": {"mutation_prob_0_to_1": q,
            "mutation_prob_1_to_0": r, "stationary_P_1": stationary,
            "stationary_flux_each_direction": (1-stationary)*q},
        "C12_precedence_vs_static_latent": urn_fixture(),
        "C13_zero_to_one": {"birth_rate_per_step": F(1, 4), "steps": 3,
            "P_born": 1-F(3, 4)**3, "zero_rate_P_born": 0,
            "literal_empty_domain_supports_normalized_probability": False},
        "C14_objective_price": {"branches": [0, 1, 2],
            "cost_A": [0, 1, 4], "cost_B": [4, 1, 0],
            "minimizer_A": 0, "minimizer_B": 2},
        "C15_drift": {"initial_dimensionless_parameter": F(1, 137),
            "increment": F(1, 1000000), "steps": 4,
            "parameter_after_4": F(1, 137)+F(4, 1000000),
            "fixed_augmented_update": "(x,c)->(c,c+delta); record y=x"},
    }
