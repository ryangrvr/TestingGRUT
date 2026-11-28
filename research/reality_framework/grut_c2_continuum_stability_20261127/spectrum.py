"""Bound the continuum negative index by a verified nodal comparison form."""
import numpy as np
from mpmath import iv

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


def interval_pivot_count(nodes, comparison):
    """Enclose every scalar or two-by-two LDL pivot of the lower form."""
    iv.dps = 55
    mesh = [iv.mpf(float(x)) for x in nodes]
    elements = []
    for j, (left, right) in enumerate(zip(mesh[:-1], mesh[1:])):
        h = right-left
        kinetic = 1/h
        mass = [h*iv.mpf(float(v[j]))/6 for v in comparison]
        if len(mass) == 1:
            elements.append((kinetic+2*mass[0], -kinetic+mass[0]))
        else:
            a, b, c = mass
            elements.append(((kinetic+2*a, 2*b, kinetic+2*c),
                             (-kinetic+a, b, -kinetic+c)))
    negative = 0
    if len(comparison) == 1:
        schur = iv.mpf(0)
        for i, (diagonal, off) in enumerate(elements):
            pivot = diagonal - schur
            if i+1 < len(elements):
                pivot += elements[i+1][0]
            if pivot < 0:
                negative += 1
            elif not pivot > 0:
                raise ArithmeticError(f"Scalar pivot {i} includes zero")
            if i+1 < len(elements):
                schur = elements[i+1][1]**2/pivot
        return negative

    schur = tuple(iv.mpf(0) for _ in range(3))
    for i, (diagonal, off) in enumerate(elements):
        adjacent = elements[i+1][0] if i+1 < len(elements) else (0, 0, 0)
        a, b, c = [x+y-z for x, y, z in zip(diagonal, adjacent, schur)]
        determinant = a*c-b*b
        if determinant < 0:
            negative += 1
        elif determinant > 0:
            if a+c < 0:
                negative += 2
            elif not a+c > 0:
                raise ArithmeticError(f"Block pivot {i} trace includes zero")
        else:
            raise ArithmeticError(f"Block pivot {i} determinant includes zero")
        if i+1 < len(elements):
            u, v, w = c/determinant, -b/determinant, a/determinant
            x, y, z = elements[i+1][1]
            schur = (u*x*x+2*v*x*y+w*y*y,
                     u*x*y+v*(x*z+y*y)+w*y*z,
                     u*y*y+2*v*y*z+w*z*z)
    return negative


def count(nodes, h, potentials, ell, level, exterior_floor):
    if exterior_floor.lo <= level:
        raise ArithmeticError("Exterior vacuum reaches the spectral test level")
    comparison = lower_cells(potentials, nodes, h, ell, level)
    return interval_pivot_count(nodes, comparison)


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
