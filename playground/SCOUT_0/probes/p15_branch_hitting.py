# P-15: two-branch absorbing diffusion dp = b(p) dt + sigma(p) dW on [0,1]; selection = absorption, branch-1 weight
# h(p) = P_p[p_t -> 1]. Exact via scale function; Monte Carlo ONLY to verify exact hitting probabilities.
# FROZEN BEFORE EVALUATION: D = 1; control B drift b = lam p(1-p)(2p-1) with lam = 2 (c = lam/D = 2);
# control D: gamma = 4, p0 = 1/2; control E: kappa = 1/4. Requires numpy, sympy.
import numpy as np, sympy as sp
rng = np.random.default_rng(15)
D, LAM, GAM, P0, KAP = 1.0, 2.0, 4.0, 0.5, 0.25

print("=== PRE-FLIGHT: does the S2 C-B variable define a branch population? ===")
Kb = np.diag([2.3]*22 + [1.3]) - np.eye(23, k=1) - np.eye(23, k=-1)
print(f"  C-B potential V = x.K_b.x/2 + beta sum x^4: min eig K_b = {np.linalg.eigvalsh(Kb).min():.4f} > 0, beta >= 0 => V strictly convex")
print("  => unique equilibrium x = 0, no two attractors, no bounded [0,1] coordinate, no canonical map to |alpha|^2: UNFORMULABLE-FROM-S2")

print("\n=== CORE: backward equation b h' + (1/2) sigma^2 h'' = 0, h(0)=0, h(1)=1 ===")
p, y, c = sp.symbols('p y c', real=True)
b, s2 = sp.Function('b'), sp.Function('s2')
h = sp.Function('h')
print("  Born h(p) = p  =>  L p = b(p) * 1 + (1/2) sigma^2 * 0 = b(p)  =>  Born for all p in (0,1)  <=>  b == 0 on (0,1)")
print("  scale function s'(p) = exp(-int 2b/sigma^2);  h(p) = s(p)/s(1) when s(0+), s(1-) finite and exit is a.s.")

print("\n=== A: driftless Wright-Fisher-type, sigma^2 = 2 D p(1-p) ===")
T = sp.Function('T')
pp = sp.symbols('pp', positive=True)
Tc = -(pp*sp.log(pp) + (1-pp)*sp.log(1-pp))/D
ok = sp.simplify(D*pp*(1-pp)*sp.diff(Tc, pp, 2) + 1) == 0
lim = (sp.limit(Tc, pp, 0, '+'), sp.limit(Tc, pp, 1, '-'))
print("  h(p) = p exactly (s' = 1).  Mean absorption time E[tau] = -(p ln p + (1-p) ln(1-p))/D: solves D p(1-p) T'' = -1:",
      ok, ", T(0+), T(1-) =", lim, " (finite: exit boundaries, absorption in finite time)")
print("  A': sigma = 4 sqrt(kappa) p(1-p) (QSD diffusive unraveling): speed density ~ 1/p^2 near 0 is non-integrable -> natural")
print("      boundaries: absorption only ASYMPTOTIC (p_t -> {0,1} a.s. by bounded-martingale convergence), still h(p) = p")

print("\n=== B: frozen nonlinear drift b = lam p(1-p)(2p-1), same sigma^2 = 2 D p(1-p); c = lam/D = 2 ===")
cval = LAM/D
sprime = sp.exp(cval*y*(1-y))          # exp(-int 2b/sigma^2) = exp(-int c(2y-1)) = exp(c y(1-y))
import mpmath as mp
mp.mp.dps = 30
sfun_ = lambda t: mp.e**(cval*t*(1-t))
S1 = mp.quad(sfun_, [0, 1])
hB = lambda pv: float(mp.quad(sfun_, [0, pv])/S1)   # exact scale-function integral (erfi form), real quadrature
pts = [0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9]
print("  exact h_B(p):", ", ".join(f"h({q})={hB(q):.5f}" for q in pts))
print("  h_B(p) - p   :", ", ".join(f"{hB(q)-q:+.5f}" for q in pts), " -> SELECTION-NON-BORN (majority amplified)")

print("\n=== Monte Carlo verification of exact hitting probabilities (Euler-Maruyama, absorbing clip) ===")
def hit(p0, bfun, sfun, n=20000, dt=2e-4, tmax=60.0):
    x = np.full(n, p0); alive = np.ones(n, bool); up = np.zeros(n, bool); t = 0.0
    while alive.any() and t < tmax:
        xa = x[alive]
        xa = xa + bfun(xa)*dt + sfun(xa)*np.sqrt(dt)*rng.standard_normal(xa.size)
        x[alive] = xa
        done = (x <= 0) | (x >= 1)
        up |= alive & (x >= 1); alive &= ~done; t += dt
    return up.mean(), alive.mean()
sigA = lambda x: np.sqrt(2*D*np.clip(x*(1-x), 0, None))
for q in (0.25, 0.6):
    fa, la = hit(q, lambda x: 0*x, sigA)
    fb, lb = hit(q, lambda x: LAM*x*(1-x)*(2*x-1), sigA)
    se = np.sqrt(q*(1-q)/20000)
    print(f"  p0={q}: A MC={fa:.4f} exact={q:.4f} (|d|/se={abs(fa-q)/se:.1f}) | B MC={fb:.4f} exact={hB(q):.4f} (|d|/se={abs(fb-hB(q))/se:.1f})  unabsorbed: {la:.1e}, {lb:.1e}")

print("\n=== C: noise off (sigma = 0): p' = lam p(1-p)(2p-1) ===")
print("  fixed points 0, 1/2 (unstable for lam>0), 1; p0>1/2 -> 1, p0<1/2 -> 0, approach exponential (asymptotic, never in finite time)")
print("  'probability of branch 1' = mu({p0 > 1/2}) requires a SUPPLIED measure mu on initial conditions: DETERMINISTIC-BASINS / MEASURE-REQUIRED")

print("\n=== D: noise without absorption: b = gamma(p0 - p), sigma^2 = 2 D p(1-p), gamma = 4, p0 = 1/2 ===")
a_ = GAM*P0/D; b_ = GAM*(1-P0)/D
print(f"  stationary density ~ p^({a_}-1)(1-p)^({b_}-1) = Beta({a_},{b_}); exponents >= 1 => both boundaries unattainable")
fd, ld = hit(0.25, lambda x: GAM*(P0-x), sigA, n=4000, tmax=20.0)
print(f"  MC from p0=0.25 to t=20: absorbed fraction {1-ld:.4f} (-> 0): NO-SELECTION (stochastic, fluctuating, no definite outcome)")

print("\n=== E: same dephasing Lindblad (L = sqrt(kappa) sigma_z), two unravelings ===")
n, dt, tt = 4000, 1e-3, 4.0
steps = int(tt/dt)
ph = np.zeros(n); pq = np.full(n, 0.3)
for _ in range(steps):
    ph += 2*np.sqrt(KAP)*np.sqrt(dt)*rng.standard_normal(n)                       # random-unitary (phase kicks)
    pq += 4*np.sqrt(KAP)*pq*(1-pq)*np.sqrt(dt)*rng.standard_normal(n); pq = np.clip(pq, 0, 1)   # QSD
coh_ru = abs(np.mean(np.exp(1j*ph)))*np.sqrt(0.3*0.7); coh_qsd = np.mean(np.sqrt(pq*(1-pq)))
print(f"  ensemble |rho_01|(t=4): random-unitary {coh_ru:.4f}, QSD {coh_qsd:.4f}, Lindblad exact {np.sqrt(0.21)*np.exp(-2*KAP*tt):.4f}")
print(f"  populations: random-unitary p_t = 0.3 on every trajectory; QSD: mean {pq.mean():.4f}, fraction within 0.05 of 0 or 1: {np.mean((pq<0.05)|(pq>0.95)):.3f}")
print("  -> identical decoherence (same ensemble master equation); selection present only in the QSD unraveling")

print("\n=== F: decomposition independence <=> Born (mixtures with the same rho_11) ===")
mix = 0.5*hB(0.2) + 0.5*hB(0.6)
print(f"  B: pure p=0.4 -> h={hB(0.4):.5f};  50/50 mixture of p=0.2, 0.6 (same rho_11 = 0.4) -> {mix:.5f}: outcome frequencies depend on the decomposition")
print("  A: h(p)=p -> both give 0.4. h decomposition-independent for all mixtures <=> h affine <=> (with h(0)=0, h(1)=1) h(p)=p")
