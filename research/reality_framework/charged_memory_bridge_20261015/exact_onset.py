"""Unregulated integer check of the first charged amplitude and its return.

No rotor cutoff, floating arithmetic, permutation quotient or cached basis is
used. One-way evolution removes P-to-complement hops only; it is a diagnostic
triangular generator, not an alternative unitary physical model.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial
from pathlib import Path
import importlib.util
import json
import sys

P = Path(__file__).resolve().parent
RF = P.parent
sys.path.insert(0, str(RF/'cp12_ground_connectivity'))
from exact_destination import sea, transitions

spec = importlib.util.spec_from_file_location('charged_record_onset',
    RF/'constructive_sprint_20261015/protected_geometry/neutral_record.py')
record = importlib.util.module_from_spec(spec)
spec.loader.exec_module(record)


@lru_cache(None)
def target(state):
    # Independent directed-degree definition, checked against the inherited
    # ring-walk implementation on every state visited by this recurrence.
    g, bits, E = state
    outgoing, incoming = [0]*4, [0]*4
    valid = bool(bits&15)
    adjacency = [set() for _ in range(4)]
    for e, (a, b) in enumerate(combinations(range(4), 2)):
        if not (g>>e&1):
            continue
        adjacency[a].add(b); adjacency[b].add(a)
        ea, eb = E[2*e:2*e+2]
        if ea*eb >= 0:
            valid = False
        for v, field in [(a, ea), (b, eb)]:
            outgoing[v] += int(field > 0)
            incoming[v] += int(field < 0)
    seen, todo = {0}, [0]
    for a in todo:
        for b in adjacency[a]-seen:
            seen.add(b); todo.append(b)
    answer = valid and outgoing == incoming == [1]*4 and len(seen) == 4
    assert answer == (bool(bits&15) and record.circulation(state))
    return answer


def advance(vector, one_way=False):
    out = {}
    for state, amplitude in vector.items():
        inside = target(state)
        for destination, coefficient in transitions(4, state, 120):
            if one_way and inside and not target(destination):
                continue
            out[destination] = out.get(destination, 0) + amplitude*coefficient
    return {s: a for s, a in out.items() if a}


def main():
    parents = [sum(1<<e for e in es) for es in combinations(range(6), 4)]
    full = {sea(4, g): 1 for g in parents}
    forward = full.copy()
    rows, selected_full, selected_forward = [], {}, {}
    for k in range(12):
        fp = {s: a for s, a in full.items() if target(s)}
        gp = {s: a for s, a in forward.items() if target(s)}
        rows.append(dict(order=k, reachable_states=len(full),
            target_states=len(fp), full_and_first_transfer_equal=fp==gp))
        if k in [9, 10, 11]:
            selected_full[k], selected_forward[k] = fp, gp
        if k < 11:
            full, forward = advance(full), advance(forward, True)
    assert all(r['target_states']==0 for r in rows[:9])
    assert rows[9]['target_states'] == 24
    assert all(r['full_and_first_transfer_equal'] for r in rows[:11])
    assert not rows[11]['full_and_first_transfer_equal']
    def dot(a, b):
        return sum(x*b.get(s, 0) for s, x in a.items())
    f9, f10, f11 = [selected_full[k] for k in [9, 10, 11]]
    g11 = selected_forward[11]
    p18 = Fraction(dot(f9, f9), 15*120**18*factorial(9)**2)
    common = Fraction(dot(f10, f10), 15*120**20*factorial(10)**2)
    p20 = common - Fraction(2*dot(f9, f11),
        15*120**20*factorial(9)*factorial(11))
    g20 = common - Fraction(2*dot(f9, g11),
        15*120**20*factorial(9)*factorial(11))
    assert p18 > 0 and p20-g20 < 0
    exited = {s: a for s, a in advance(f9).items() if not target(s)}
    return_norm_ratio = Fraction(dot(exited, exited), 120**2*dot(f9, f9))
    assert return_norm_ratio == 6
    assert p20-g20 == -Fraction(2, 10*11)*return_norm_ratio*p18
    result = dict(status='PASS', N=4, Q=4, normalized_parent_count=15,
        event='All four active edges in one directed spanning circulation; nv>=1.',
        exact_integer_full_unregulated=True, rows=rows,
        independent_directed_degree_event_check=True,
        first_amplitude_order=9, first_return_amplitude_difference_order=11,
        leading_probability_order=18, leading_probability_coefficient=str(p18),
        probability_order20_coefficient=str(p20),
        first_transfer_probability_order20_coefficient=str(g20),
        full_minus_first_transfer_order20_coefficient=str(p20-g20),
        leading_vector_Cdagger_norm_square_over_norm_square=str(return_norm_ratio),
        general_first_return_identity_verified=True,
        probability_remainder='O(T^22): real H and real preparation imply even probability.',
        scope='A new finite-control coefficient; does not establish a large-N charged phase.')
    (P/'EXACT_ONSET.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
