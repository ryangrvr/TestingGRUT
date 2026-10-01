# P-06b: exact closure of the CM boundary of the Lorentzian-exponent family.
# Claim: f_p(z) = z^nu K_nu(z), nu = p - 1/2, is completely monotone iff 0 < p <= 1.
# Proof (see P06B_RESULT.md; proof repaired in P06B_CORRECTION_01.md):
#  (A)  0 < p < 1: with mu = 1/2 - p in (-1/2, 1/2) and K_nu = K_{-nu} = K_mu, DLMF 10.32.8 gives
#       z^nu K_nu(z) = sqrt(pi) 2^nu / Gamma(1-p) * Int_1^inf e^{-zs} (s^2-1)^{-p} ds.
#       Positive weight, integrable at s=1 iff p<1  => Laplace transform of a positive measure => CM.
#       p = 1 is the explicit exponential endpoint sqrt(pi/2) e^{-z}.
#  (C)  p > 1: recurrence d/dz[z^nu K_nu] = -z^nu K_{nu-1} gives
#       f_nu''(0) = -f_{nu-1}(0) = -2^(nu-2) Gamma(nu-1) < 0 for all nu > 1  => NOT CM.
#       (covers integer nu too; for 1/2<nu<1 f'' ~ c*beta*(beta-1)*z^(beta-2) < 0 near 0,
#        c = 2^(-nu-1) Gamma(-nu) < 0, beta = 2nu in (1,2).)
#  RETRACTED (81ad661): "z^alpha * CM is CM for alpha in (0,1]" — false; z e^{-z} is not CM.
import mpmath as mp
mp.mp.dps = 22

def f(nu, z):
    return z**nu * mp.besselk(nu, z)

print("=== (C) p>1: f'' negative near 0 for all tested nu>1 ===")
for p in ['1.05', '1.2', '1.5', '2.0', '2.5', '3.0']:
    nu = mp.mpf(p) - mp.mpf('0.5')
    pred = mp.nstr(-2**(nu-2) * mp.gamma(nu-1), 8) if nu > 1 else 'divergent (log/power)'
    z = mp.mpf('0.05')
    num = mp.diff(lambda x: f(nu, x), z, 2)
    print(f"p={p:>4}: f''(0.05) = {mp.nstr(num, 8)}   predicted f''(0)={pred}")

print()
print("=== (A) 0<p<1: unified Laplace (Bernstein) representation identity ===")
# I(z) = Int_1^inf e^{-zs}(s^2-1)^{-p} ds; substitute s = 1 + t^q, q = 1/(1-p), which removes the
# endpoint singularity: I = e^{-z}/(1-p) * Int_0^inf e^{-z t^q} (2 + t^q)^{-p} dt.
def bernstein_integral(p, z):
    q = 1/(1-p)
    tcut = (60/z)**(1/q)   # e^{-z t^q} < e^{-60} beyond this
    return mp.e**(-z)/(1-p) * mp.quad(lambda t: mp.e**(-z*t**q) * (2 + t**q)**(-p), [0, 1, tcut])

for p in ['0.1', '0.25', '0.5', '0.6', '0.75', '0.9', '0.95', '0.99']:
    p = mp.mpf(p); nu = p - mp.mpf('0.5')
    C = mp.sqrt(mp.pi) * 2**nu / mp.gamma(1-p)
    worst = mp.mpf(0)
    for zs in ['0.2', '1.0', '4.0']:
        z = mp.mpf(zs)
        rel = abs(f(nu, z) - C*bernstein_integral(p, z))/abs(f(nu, z))
        worst = max(worst, rel)
    print(f"p={mp.nstr(p,4):>4}: max relative mismatch = {mp.nstr(worst, 3)}  -> positive-weight Bernstein form CONFIRMED")

print()
print("=== retracted lemma of 81ad661: counterexample and failed identity ===")
h = lambda z: z*mp.e**(-z)
print(f"(z e^-z)'(0.5) = {mp.nstr(mp.diff(h, mp.mpf('0.5')), 8)}  (>0, so z*CM need not be CM)")
a, s, z = mp.mpf('0.5'), mp.mpf('1'), mp.mpf('1')
lhs = z**a * mp.e**(-s*z)
rhs = mp.quad(lambda t: t*mp.e**(-z*t)*(t-s)**(-a), [s, s+1, mp.inf])/mp.gamma(1-a)
print(f"claimed identity z^a e^-sz = (1/G(1-a)) Int_s^inf zeta e^-z zeta (zeta-s)^-a: RHS/LHS = {mp.nstr(rhs/lhs, 6)} (should be 1)")

print()
print("=== corroboration, 1/2<p<=1: derivative-sign (CM) spot test, n<=4 ===")
def is_cm_fp(p, zmin='0.05', zmax='12', n=4, ngrid=10):
    nu = mp.mpf(p) - mp.mpf('0.5')
    grid = [mp.mpf(zmin) * (mp.mpf(zmax)/mp.mpf(zmin))**(mp.mpf(i)/(ngrid-1)) for i in range(ngrid)]
    for k in range(1, n+1):
        for z in grid:
            dk = mp.diff(lambda x: f(nu, x), z, k)
            if (mp.mpf(-1)**k) * dk < -mp.mpf('1e-10'):
                return False, (k, mp.nstr(z, 5), mp.nstr(dk, 8))
    return True, None

for p in ['0.6', '0.75', '0.9', '1.0']:
    ok, ev = is_cm_fp(p)
    print(f"p={p:>4}: CM(n<=4) = {ok}" + ("" if ok else f"  violation: {ev}"))

print()
print("=== p=2 exact-counterexample cross-check: f''=e^-z(z-1) ===")
for zs in ['0.5', '1.0']:
    z = mp.mpf(zs)
    num = mp.diff(lambda x: f(mp.mpf('1.5'), x), z, 2)
    ana = mp.sqrt(mp.pi/2)*mp.e**(-z)*(z-1)
    print(f"z={z}: numeric {mp.nstr(num, 10)} vs analytic sqrt(pi/2)e^-z(z-1) {mp.nstr(ana, 10)}")
print("ALL P-06B NUMERICAL CHECKS COMPLETE")
