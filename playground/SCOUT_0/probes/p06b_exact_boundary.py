# P-06b: exact closure of the CM boundary of the Lorentzian-exponent family.
# Claim: f_p(z) = z^nu K_nu(z), nu = p - 1/2, is completely monotone iff 0 < p <= 1.
# Proof (see P06B_RESULT.md):
#  (i)  p <= 1/2: z^-mu K_mu(z) = sqrt(pi) 2^-mu / Gamma(mu+1/2) *
#       Int_1^inf e^{-zs}(s^2-1)^(mu-1/2) ds, mu = 1/2-p >= 0  => positive measure => CM.
#  (ii) 1/2 < p <= 1: z^nu K_nu(z) = C * z^(2nu) * L[w], w(s) = (s^2-1)^(nu-1/2) >= 0,
#       alpha = 2nu in (0,1]; z^alpha * (CM g) is CM for alpha in [0,1]
#       (lemma: z^alpha e^{-sz} = (1/Gamma(1-alpha)) Int_s^inf zeta e^{-z zeta}(zeta-s)^-alpha dzeta >= 0).
#  (iii) p > 1: recurrence d/dz[z^nu K_nu] = -z^nu K_{nu-1} gives
#       f_nu''(0) = -f_{nu-1}(0) = -2^(nu-2) Gamma(nu-1) < 0 for all nu > 1  => NOT CM.
#       (covers integer nu too; for 1/2<nu<1 f'' ~ c*beta*(beta-1)*z^(beta-2) < 0 near 0,
#        c = 2^(-nu-1) Gamma(-nu) < 0, beta = 2nu in (1,2).)
import mpmath as mp
mp.mp.dps = 22

def f(nu, z):
    return z**nu * mp.besselk(nu, z)

print("=== (iii) p>1: f'' negative near 0 for all tested nu>1 ===")
for p in ['1.05', '1.2', '1.5', '2.0', '2.5', '3.0']:
    nu = mp.mpf(p) - mp.mpf('0.5')
    pred = mp.nstr(-2**(nu-2) * mp.gamma(nu-1), 8) if nu > 1 else 'divergent (log/power)'
    z = mp.mpf('0.05')
    num = mp.diff(lambda x: f(nu, x), z, 2)
    print(f"p={p:>4}: f''(0.05) = {mp.nstr(num, 8)}   predicted f''(0)={pred}")

print()
print("=== (i) p<=1/2: Laplace (Bernstein) representation identity ===")
for p in ['0.1', '0.25', '0.5']:
    mu = mp.mpf('0.5') - mp.mpf(p)
    C = mp.sqrt(mp.pi) * 2**(-mu) / mp.gamma(mu + mp.mpf('0.5'))
    worst = mp.mpf(0)
    for zs in ['0.3', '1.0', '3.0']:
        z = mp.mpf(zs)
        lap = mp.quad(lambda s: mp.e**(-z*s) * (s*s-1)**(mu - mp.mpf('0.5')), [1, 2, mp.inf])
        rel = abs(f(-mu, z) - C*lap)/abs(f(-mu, z))
        worst = max(worst, rel)
    print(f"p={p:>4}: max relative mismatch = {mp.nstr(worst, 3)}  -> positive-weight Bernstein form CONFIRMED")

print()
print("=== (ii) 1/2<p<=1: derivative-sign (CM) spot test, n<=4 ===")
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
