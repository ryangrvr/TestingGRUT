"""Outward interval operations and cell polynomials for the C2 verifier."""
from itertools import zip_longest
from math import comb

import numpy as np


LD = np.longdouble


def down(x):
    return np.nextafter(x, LD(-np.inf))


def up(x):
    return np.nextafter(x, LD(np.inf))


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=LD)
        self.hi = self.lo if hi is None else np.asarray(hi, dtype=LD)

    @staticmethod
    def of(value):
        return value if isinstance(value, Interval) else Interval(value)

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __add__(self, value):
        b = self.of(value)
        return Interval(down(self.lo + b.lo), up(self.hi + b.hi))

    __radd__ = __add__

    def __sub__(self, value):
        return self + (-self.of(value))

    def __rsub__(self, value):
        return self.of(value) + (-self)

    def __mul__(self, value):
        b = self.of(value)
        products = [a * c for a in (self.lo, self.hi) for c in (b.lo, b.hi)]
        return Interval(down(np.minimum.reduce(products)), up(np.maximum.reduce(products)))

    __rmul__ = __mul__

    def __truediv__(self, value):
        b = self.of(value)
        if np.any((b.lo <= 0) & (b.hi >= 0)):
            raise ArithmeticError("An interval divisor contains zero")
        return self * Interval(down(LD(1) / b.hi), up(LD(1) / b.lo))

    def square(self):
        at_zero = (self.lo <= 0) & (self.hi >= 0)
        return Interval(np.where(at_zero, 0, down(np.minimum(self.lo * self.lo, self.hi * self.hi))),
                        up(np.maximum(self.lo * self.lo, self.hi * self.hi)))

    def sqrt(self):
        if np.any(self.lo < 0):
            raise ArithmeticError("Negative lower endpoint for square root")
        return Interval(down(np.sqrt(self.lo)), up(np.sqrt(self.hi)))

    def magnitude(self):
        return Interval(np.maximum(abs(self.lo), abs(self.hi)))

    def __getitem__(self, index):
        return Interval(self.lo[index], self.hi[index])

    def sum_cells(self):
        lo, hi = self.lo, self.hi
        while lo.size > 1:
            n = lo.size // 2
            a = down(lo[:2*n:2] + lo[1:2*n:2])
            b = up(hi[:2*n:2] + hi[1:2*n:2])
            lo = np.concatenate((a, lo[2*n:]))
            hi = np.concatenate((b, hi[2*n:]))
        return Interval(lo[0], hi[0])

    def display(self):
        return [float(np.nextafter(float(self.lo), -np.inf)),
                float(np.nextafter(float(self.hi), np.inf))]


I = Interval
ZERO = I(0)


def add(a, b):
    return [x + y for x, y in zip_longest(a, b, fillvalue=ZERO)]


def scale(a, factor):
    return [x * factor for x in a]


def multiply(a, b):
    result = [ZERO] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] = result[i+j] + x * y
    return result


def hessian(a, h):
    return [(i * (i-1) * a[i]) / h.square() for i in range(2, len(a))]


def hermite(values, first, second, h):
    """Cells use the same nodal values and first two radial derivatives."""
    v, d, dd = I(values), I(first), I(second)
    a, b, c = v[:-1], h * d[:-1], (h.square() * dd[:-1]) / 2
    ev = v[1:] - (a + b + c)
    ed = h * d[1:] - (b + 2 * c)
    edd = h.square() * dd[1:] - 2 * c
    return [a, b, c, 10 * ev - 4 * ed + edd / 2,
            -15 * ev + 7 * ed - edd, 6 * ev - 3 * ed + edd / 2]


def hull(coeffs):
    degree = len(coeffs) - 1
    bernstein = [
        sum((coeffs[i] * (I(comb(j, i)) / comb(degree, i)) for i in range(j+1)), ZERO)
        for j in range(degree + 1)
    ]
    return I(np.minimum.reduce([b.lo for b in bernstein]), np.maximum.reduce([b.hi for b in bernstein]))


def integral(coeffs, h):
    return (h * sum((coefficient / (i+1) for i, coefficient in enumerate(coeffs)), ZERO)).sum_cells()


def residual_norm(defects, h):
    cell_bound = h * sum((hull(poly).magnitude().square() for poly in defects), ZERO)
    return I(0, cell_bound.sum_cells().hi).sqrt()
