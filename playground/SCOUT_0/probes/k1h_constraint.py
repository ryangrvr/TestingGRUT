# K1-H step 1: independent symbolic re-derivation of the exact K1 condition from Linder (2020) Eqs. for luminal
# Horndeski (PRIMARY-TEXT formulas as supplied by the external audit; ' = d/d ln a; R = m_p^2/M_*^2, R' = -alpha_M R).
import sympy as sp
aB, aM, aBp, R = sp.symbols('alpha_B alpha_M alpha_Bp R')
A = aB + 2*aM; D = (2 - aB)*A + 2*aBp
mu = R*((2 + 2*aM)*A + 2*aBp)/D
Sg = R*((2 + aM)*A + 2*aBp)/D
lock = sp.factor(sp.together(2*Sg - mu - 1))
num = sp.expand(sp.numer(lock))
target = sp.expand(2*(R - 1)*aBp + (2*(R - 1) + aB)*(aB + 2*aM))
print("2*Sigma - mu =", sp.simplify(2*Sg - mu), "   (= R(2A + 2aB')/D)")
print("numerator of 2Sigma - mu - 1 :", sp.factor(num))
print("numerator == -(audit form 2(R-1)aB' + [2(R-1)+aB](aB+2aM)) :", sp.simplify(num + target) == 0, " (same zero set)")
print("at R = 1:", sp.factor(target.subs(R, 1)))
print("Only Run (aB = aB' = 0): 2Sigma - mu =", sp.simplify((2*Sg - mu).subs({aB: 0, aBp: 0})))
print("No Run (aM = 0): Sigma - mu =", sp.simplify((Sg - mu).subs(aM, 0)), " (no slip)")
eps = sp.symbols('eps')
r, m = sp.symbols('r m', positive=True)
expr = ((Sg - 1)/(mu - 1)).subs({R: 1, aBp: 0}).subs({aB: eps*r*m, aM: eps*m}, simultaneous=True)
small = sp.limit(sp.simplify(expr), eps, 0)
print("small-alpha, R=1, aB'=0: (Sigma-1)/(mu-1) ->", sp.simplify(small), " (the (1+r)/(2+r) form used in K1-OCC is exact only at R=1, aB'=0, leading order)")
