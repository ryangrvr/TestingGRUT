#!/usr/bin/env python3
"""V0-1: symbolic checks of the bath elimination, kernel, free force, slip and free-force laws."""
from sympy import (symbols, Function, cos, sin, integrate, diff, simplify, sqrt, pi, Rational,
                   exp, expand, trigsimp, Integral, Symbol)

t, s, u, w, c, q0, xi, eta, T, phi = symbols("t s u omega c q0 xi eta T phi", real=True)
w = Symbol("omega", positive=True)
q = Function("q")

# (a) bath solution given system path: x(t) = x0 cos wt + p0/w sin wt + (c/w) int_0^t sin w(t-s) q(s) ds
x0, p0 = symbols("x0 p0", real=True)
X = x0 * cos(w * t) + p0 / w * sin(w * t) + c / w * Integral(sin(w * (t - s)) * q(s), (s, 0, t))
ode = diff(X, t, 2) + w**2 * X - c * q(t)
print("(a) x'' + w^2 x - c q(t) =", simplify(ode.doit()))
print("    x(0) - x0 =", simplify(X.subs(t, 0).doit() - x0),
      "; x'(0) - p0 =", simplify(diff(X, t).doit().subs(t, 0) - p0))

# (b) integration by parts identity, checked on concrete paths (exact)
for qq in [s**3 - 2 * s + 1, exp(s) * s, exp(-2 * s) + s**4]:
    lhs = c / w * integrate(sin(w * (t - s)) * qq, (s, 0, t))
    rhs = (c / w**2) * (qq.subs(s, t) - qq.subs(s, 0) * cos(w * t)) \
        - (c / w**2) * integrate(cos(w * (t - s)) * diff(qq, s), (s, 0, t))
    print("(b) IBP identity residual for q =", qq, ":", simplify(lhs - rhs))

# (c) force on system: c (x - c q/w^2) = F(t) - int gamma(t-s) q'(s) ds with xi = x0 - c q0/w^2
Fj = c * ((x0 - c * q0 / w**2) * cos(w * t) + p0 / w * sin(w * t))
print("(c) free force per mode F_j(t) =", Fj)

# (g) thermal shifted Gibbs: xi ~ N(0,T/w^2), eta ~ N(0,T) indep -> Cov = c^2 [T/w^2 cos cos + T/w^2 sin sin]
cov = c**2 * (T / w**2 * cos(w * t) * cos(w * s) + T / w**2 * sin(w * t) * sin(w * s))
print("(g) Cov F_j(t)F_j(s) - T c^2/w^2 cos w(t-s) =", simplify(trigsimp(cov - T * c**2 / w**2 * cos(w * (t - s)))))

# random-phase: xi = sqrt(2T)/w cos phi, eta = -sqrt(2T) sin phi
FR = c * (sqrt(2 * T) / w * cos(phi) * cos(w * t) - sqrt(2 * T) * sin(phi) / w * sin(w * t))
print("    random-phase F_j(t) - c sqrt(2T)/w cos(wt+phi) =", simplify(trigsimp(FR - c * sqrt(2 * T) / w * cos(w * t + phi))))
E = lambda f: simplify(integrate(expand(f), (phi, 0, 2 * pi)) / (2 * pi))
print("    random-phase mean:", E(FR))
covR = E(FR * FR.subs(t, s))
print("    random-phase Cov - T c^2/w^2 cos w(t-s):", simplify(trigsimp(covR - T * c**2 / w**2 * cos(w * (t - s)))))
F0 = FR.subs(t, 0)
k4 = simplify(E(F0**4) - 3 * E(F0**2) ** 2)
print("    random-phase 4th cumulant at equal times (t=0):", k4, "  (Gaussian: 0)")

# (d) product prep: x0 ~ N(0,T/w^2) => xi = x0 - c q0/w^2, F = F_thermal - (c^2/w^2) cos(wt) q0
print("(d) product-prep mean of free force per mode:", simplify(c * (-c * q0 / w**2) * cos(w * t)),
      " = -gamma_j(t) q0")
