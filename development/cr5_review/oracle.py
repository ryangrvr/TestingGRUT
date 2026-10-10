"""Separate exact arithmetic implementation of issued G2-14R items 4–10.

Inputs represent ALREADY PROVED all-branch enclosures. This module neither
constructs nor approves a physical limit theorem or an experimental lock.
"""
from fractions import Fraction as F


def exact(x):
    if type(x) not in (int, F, str):
        raise ValueError('Exact rational input required')
    return F(x)


def box(values):
    if not values:
        raise ValueError('Nonempty finite frozen chart required')
    result = {}
    for key, pair in values.items():
        center, error = map(exact, pair)
        if error < 0:
            raise ValueError('Negative enclosure error')
        result[key] = (center, error)
    return result


def reference(probe, bounds):
    probe = box(probe)
    if set(probe) != set(bounds):
        raise ValueError('Incomplete whole-chart error bound')
    result = {}
    for k, (center, error) in probe.items():
        extra = exact(bounds[k])
        if extra < 0:
            raise ValueError('Nonnegative proven bound required')
        result[k] = (center, error + extra)
    return result


def decide(current, limit, tolerance):
    current, limit = box(current), box(limit)
    if set(current) != set(limit) or set(current) != set(tolerance):
        raise ValueError('Frozen chart mismatch')
    upper, lower = [], []
    for key in current:
        c, e = current[key]
        r, f = limit[key]
        t = exact(tolerance[key])
        if t < 0:
            raise ValueError('Negative physical tolerance')
        upper.append(abs(c-r) + e + f <= t)
        lower.append(max(F(0), abs(c-r)-e-f) > t)
    return 'EXCLUDE' if all(upper) else 'INCLUDE' if any(lower) else 'UNRESOLVED'


def nr17(decisions, holds):
    if len(decisions) != len(holds) or any(type(h) is not bool for h in holds):
        raise ValueError('Incomplete independently certified relation decisions')
    if any(d not in ('EXCLUDE', 'INCLUDE', 'UNRESOLVED') for d in decisions):
        raise ValueError('Unknown decision')
    if 'UNRESOLVED' in decisions:
        return None
    removed = [d == 'EXCLUDE' and h for d, h in zip(decisions, holds)]
    return [False] * len(removed) if removed and all(removed) else removed


def q3(first, second):
    if any(d not in ('EXCLUDE', 'INCLUDE', 'UNRESOLVED') for d in (first, second)):
        raise ValueError('Unknown chain decision')
    if first == second == 'EXCLUDE':
        return 'EXCLUDE'
    if 'UNRESOLVED' in (first, second):
        return 'UNRESOLVED'
    return 'INCLUDE'


def measure(weights, excluded, hits):
    if len(weights) != len(excluded) or len(weights) != len(hits):
        raise ValueError('Shared domain mismatch')
    weights = list(map(exact, weights))
    if any(w < 0 for w in weights) or sum(weights) <= 0:
        raise ValueError('Reference measure requires nonnegative weights and positive total')
    if any(type(d) is not bool for d in excluded) or any(type(h) is not bool for h in hits):
        raise ValueError('Complete shared mask and appearance decisions required')
    retained = sum(w for w, d in zip(weights, excluded) if not d)
    numerator = sum(w for w, d, h in zip(weights, excluded, hits) if not d and h)
    return None if retained == 0 else numerator / retained


def transport(matrix, values):
    vals = list(box(values).values())
    rows = [list(map(exact, row)) for row in matrix]
    if not rows or any(len(row) != len(vals) for row in rows):
        raise ValueError('Transport chart mismatch')
    return [(sum(x*c for x, (c, _) in zip(row, vals)),
             sum(abs(x)*e for x, (_, e) in zip(row, vals))) for row in rows]


def controls():
    # Boundary inequalities, ambiguity, independent chains, hostile VOID semantics.
    assert decide({'x': (1, 0)}, {'x': (0, 0)}, {'x': 1}) == 'EXCLUDE'
    assert decide({'x': (1, '1/4')}, {'x': (0, 0)}, {'x': 1}) == 'UNRESOLVED'
    assert decide({'x': ('3/2', '1/4')}, {'x': (0, 0)}, {'x': 1}) == 'INCLUDE'
    assert nr17(['EXCLUDE', 'EXCLUDE'], [True, True]) == [False, False]
    assert q3('EXCLUDE', 'UNRESOLVED') == 'UNRESOLVED'
    assert measure([1, 1], [True, True], [True, True]) is None
    assert measure([1, 1], [False, False], [False, True]) == F(1, 2)
    return 7
