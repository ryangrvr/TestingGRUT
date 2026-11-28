"""Bound the continuum negative index by a verified nodal comparison form."""
import numpy as np
from mpmath import iv, mp

from intervals import I, ZERO


def matrix_bounds(potentials):
    """Singleton centers and an upper bound on the maximum matrix row error."""
    centers = [(p.lo / 2 + p.hi / 2).astype(float) for p in potentials]
    errors = [(p - I(c)).magnitude() for p, c in zip(potentials, centers)]
    if len(errors) == 1:
        return centers, errors[0]
    a, b, c = errors
    return centers, I(np.maximum((a+b).hi, (c+b).hi))


def down_double(x):
    return np.nextafter(np.asarray(x, dtype=float), -np.inf)


def lower_cells(potentials, nodes, h, ell, level):
    centers, error = matrix_bounds(potentials)
    barrier = I(ell*(ell+1)) / I(nodes[1:]).square() if ell else ZERO
    diag = [0, 2] if len(centers) == 3 else [0]
    comparison = list(centers)
    for index in diag:
        comparison[index] = down_double((I(centers[index]) - error + barrier - level).lo)
    m = I(abs(comparison[0]))
    if len(comparison) == 3:
        m = I(np.maximum((m + abs(comparison[1])).hi,
                         (I(abs(comparison[2])) + abs(comparison[1])).hi))
    floor = I(9) / h.square() - m
    if np.any(floor.lo <= 0):
        raise ArithmeticError("Nodal remainder has no positive Poincare floor")
    correction = m.square() / floor
    for index in diag:
        comparison[index] = down_double((I(comparison[index]) - correction).lo)
    return comparison


def entries(a, convert):
    return [[convert(x) for x in row] for row in a]


def transpose(a):
    return [list(column) for column in zip(*a)]


def product(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def combine(a, b, sign=1):
    return [[x+sign*y for x, y in zip(left, right)] for left, right in zip(a, b)]


def inverse_block(d):
    if len(d) == 1:
        return [[1/d[0][0]]]
    a, b, c = d[0][0], d[0][1], d[1][1]
    determinant = a*c-b*b
    return [[c/determinant, -b/determinant], [-b/determinant, a/determinant]]


def pivot_index(d):
    a = d[0][0]
    if len(d) == 1:
        if a < 0:
            return 1
        if a > 0:
            return 0
    else:
        b, c = d[0][1], d[1][1]
        determinant = a*c-b*b
        if determinant < 0:
            return 1
        if determinant > 0:
            if a+c < 0:
                return 2
            if a+c > 0:
                return 0
    raise ArithmeticError("A point factor's pivot sign is unresolved")


def norm_upper(a):
    return max(sum(abs(x) for x in row).b for row in a)


def mesh_form(nodes, comparison):
    mesh = [iv.mpf(float(x)) for x in nodes]
    elements = []
    for j, (left, right) in enumerate(zip(mesh[:-1], mesh[1:])):
        h = right-left
        kinetic = 1/h
        mass = [h*iv.mpf(float(v[j]))/6 for v in comparison]
        if len(mass) == 1:
            elements.append(([[kinetic+2*mass[0]]], [[-kinetic+mass[0]]]))
        else:
            a, b, c = mass
            elements.append(([[kinetic+2*a, 2*b], [2*b, kinetic+2*c]],
                             [[-kinetic+a, b], [b, -kinetic+c]]))
    diagonal = [combine(d, elements[i+1][0]) if i+1 < len(elements) else d
                for i, (d, _) in enumerate(elements)]
    return diagonal, [off for _, off in elements[1:]]


def factorize(diagonal, off):
    pivots, coupling = [], []
    for i, d in enumerate(diagonal):
        point = entries(d, lambda x: mp.mpf(x.mid))
        if i:
            t = coupling[-1]
            point = combine(point, product(product(t, pivots[-1]), transpose(t)), -1)
        if len(point) == 2:
            point[1][0] = point[0][1]
        pivots.append(point)
        if i < len(off):
            try:
                coupling.append(product(entries(off[i], lambda x: mp.mpf(x.mid)), inverse_block(point)))
            except ZeroDivisionError as err:
                raise ArithmeticError("Point factor has a zero pivot") from err
    return pivots, coupling


def triangular_bound(couplings, dimension):
    """Bound row sums of the inverse of a unit block bidiagonal matrix."""
    row_sums = [iv.mpf(1)] * dimension
    maximum = iv.mpf(1)
    for t in couplings:
        row_sums = [(1+sum(abs(x)*y for x, y in zip(row, row_sums))).b for row in t]
        maximum = max(maximum, *row_sums)
    return maximum


def verified_factor_count(diagonal, off):
    """Certify the exact form's inertia by an outward residual and Neumann test."""
    point_d, point_t = factorize(diagonal, off)
    d = [entries(block, iv.mpf) for block in point_d]
    t = [entries(block, iv.mpf) for block in point_t]
    negative = sum(pivot_index(pivot) for pivot in d)
    inverse_d = max(norm_upper(inverse_block(pivot)) for pivot in d)
    residual_diag, residual_lower, residual_upper = [iv.mpf(0) for _ in range(3)]
    for i, pivot in enumerate(d):
        computed = pivot
        if i:
            td = product(t[i-1], d[i-1])
            computed = combine(pivot, product(td, transpose(t[i-1])))
            residual = combine(off[i-1], td, -1)
            residual_lower = max(residual_lower, norm_upper(residual))
            residual_upper = max(residual_upper, norm_upper(transpose(residual)))
        residual_diag = max(residual_diag, norm_upper(combine(diagonal[i], computed, -1)))
    residual = (residual_diag+residual_lower+residual_upper).b
    dimension = len(d[0])
    row_norm = triangular_bound(t, dimension)
    col_norm = triangular_bound((transpose(block) for block in reversed(t)), dimension)
    product_bound = (residual*row_norm*inverse_d*col_norm).b
    if not product_bound < 1:
        raise ArithmeticError("Factorization residual times verified inverse norm reaches one")
    return negative, float(np.nextafter(float(product_bound), np.inf))


def count(nodes, h, potentials, ell, level, exterior_floor):
    if exterior_floor.lo <= level:
        raise ArithmeticError("Exterior vacuum reaches the spectral test level")
    comparison = lower_cells(potentials, nodes, h, ell, level)
    mp.dps = iv.dps = 55
    return verified_factor_count(*mesh_form(nodes, comparison))


def trial_form(trial, h, potentials, level):
    """Enclose an upper form for a compact affine Dirichlet trial."""
    centers, error = matrix_bounds(potentials)
    aa, ab, bb = [I(c) for c in centers]
    aa, bb = aa+error-level, bb+error-level
    left = [I(column[:-1]) for column in trial]
    right = [I(column[1:]) for column in trial]

    def bilinear(x, y):
        return aa*x[0]*y[0] + ab*(x[0]*y[1]+x[1]*y[0]) + bb*x[1]*y[1]

    kinetic = sum(((a-b).square()/h for a, b in zip(left, right)), ZERO)
    mass = h*(bilinear(left,left)+bilinear(left,right)+bilinear(right,right))/3
    return (kinetic+mass).sum_cells()
