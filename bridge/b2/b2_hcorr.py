"""B2 -- H_corr versus GRUT preparation / environment.  Independent finite-N re-implementation of the canonical S6 parent
(read-only evidence: S6_1_THEOREM_01.md, S6_1_COARSEGRAINED_ARROW_CHARTER_01.md, S6_OWNER_RULING_02.md).
Model: q'' = -K q, K = 2.3 I - T (half-line adjacency), system site 1, bath {2,...,N}; classical Gaussian ensembles.
N = 600, all times << recurrence time; finite-N numbers are cross-checks of the N = infinity closed forms, never adjudicating.
Admissibility of a classical Gaussian state: covariance positive definite (certified by min eigenvalue). Fixed throughout:
D (= K), Sigma (site 1 | bath), access (retained 2x2 block of site 1 and the site-2 coupling observable J), times, resolution."""
import numpy as np
np.set_printoptions(precision=5, suppress=True)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

N = 600
K = 2.3 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
lam, V = np.linalg.eigh(K); om = np.sqrt(lam)
Kinv = V @ np.diag(1 / lam) @ V.T
KBB = K[1:, 1:]; KBBinv = np.linalg.inv(KBB)
r_inf = (2.3 - np.sqrt(1.29)) / 2; rN = Kinv[0, 0]
print(f"N = {N}; spectrum [{lam[0]:.4f}, {lam[-1]:.4f}]; r_N = {rN:.10f} vs r = {r_inf:.10f}")
Z = np.zeros((N, N)); I = np.eye(N)

def M_apply(t, vq, vp, transpose=False):
    """apply M(t) = [[C, S], [-K S, C]] (or its transpose [[C, -K S], [S, C]]) to (vq, vp) columns"""
    c = np.cos(om * t)[:, None]; s = (np.sin(om * t) / om)[:, None]; ks = (om * np.sin(om * t))[:, None]
    aq, ap = V.T @ vq, V.T @ vp
    if not transpose: return V @ (c * aq + s * ap), V @ (-ks * aq + c * ap)
    return V @ (c * aq - ks * ap), V @ (s * aq + c * ap)
ROWS = {"q1": 0, "q2": 1, "p1": N, "p2": N + 1}
def rows_at(Sig, t, names):
    """Sigma(t)[names, :] = R Sigma0 M^T with R = rows of M -> computed as (M (Sigma0 R^T))^T"""
    E = np.zeros((2 * N, len(names)))
    for k, nme in enumerate(names): E[ROWS[nme], k] = 1
    Rq, Rp = M_apply(t, E[:N], E[N:], transpose=True)               # columns = M^T e_r  (rows of M)
    W = Sig @ np.vstack([Rq, Rp])                                     # Sigma0 R^T
    Wq, Wp = M_apply(t, W[:N], W[N:])
    return np.vstack([Wq, Wp]).T                                      # len(names) x 2N
def local(Sig, t):
    R3 = rows_at(Sig, t, ["q1", "p1", "p2"])
    Q11, C11, P11 = R3[0, 0], R3[0, N], R3[1, N]
    Q12, J = R3[0, 1], R3[2, 0]                                       # J = <p2 q1>
    return Q11, P11, C11, Q12, J, R3[:2]
def E1(Q11, P11): return 1.15 * Q11 + 0.5 * P11
def XJ(Sig, T):
    a, b = local(Sig, 0.0), local(Sig, T)
    return E1(a[0], a[1]) - E1(b[0], b[1]) + b[3] - a[3]              # LS-2.1, general (Q12(0) need not vanish)
def KL(S, Rf):
    Ri = np.linalg.inv(Rf); return 0.5 * (np.trace(Ri @ S) - 2 - np.log(np.linalg.det(S) / np.linalg.det(Rf)))
def Dt(Sig, t, Tb):
    Q11, P11, C11 = local(Sig, t)[:3]; return KL(np.array([[Q11, C11], [C11, P11]]), Tb * np.diag([rN, 1.0]))
def CSB(Sig, t):
    R = local(Sig, t)[5]; idx = [i for i in range(2 * N) if i not in (0, N)]
    return np.linalg.norm(R[:, idx])
def gibbs(T): return np.block([[T * Kinv, Z], [Z, T * I]])
def s6_product(Ts, Tb):
    Q = np.zeros((N, N)); Q[0, 0] = Ts / 2.3; Q[1:, 1:] = Tb * KBBinv
    P = Tb * I.copy(); P[0, 0] = Ts
    return np.block([[Q, Z], [Z, P]])
def mineig(S): return np.linalg.eigvalsh(S).min()
def entropy(S): return 0.5 * np.linalg.slogdet(S)[1]                  # Gaussian entropy up to a constant
TS = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0]
TLONG = 200.0

# =============================================================================================
hdr("B2-0 the S6 preparation as Sigma0 = Sigma_G + Delta0 (consistency of the re-implementation)")
for (Ts, Tb) in [(2, 1), (10, 1), (0.5, 1), (0.1, 1)]:
    S0 = s6_product(Ts, Tb); D0 = S0 - gibbs(Tb)
    rk = np.linalg.matrix_rank(D0, tol=1e-9)
    xj = XJ(S0, TLONG); xj_th = (Ts - Tb) + 0.5 * Tb * r_inf ** 2
    print(f"  (Ts,Tb)=({Ts},{Tb}): rank(Delta0) = {rk} (theorem: <= 3); X_J(200) = {xj:.5f} vs closed form {xj_th:.5f};"
          f" D(0) = {Dt(S0, 0, Tb):.5f}")

# =============================================================================================
hdr("B2-1 EQUAL-TEMPERATURE CONTROL: factorized product at T vs global Gibbs at T (same D, Sigma, A, T)")
T = 1.0; SA = s6_product(T, T); SB = gibbs(T)
print(f"  admissible: min eig State A {mineig(SA):.3e}, State B {mineig(SB):.3e}")
for nm, S in [("A factorized", SA), ("B Gibbs", SB)]:
    devs = [np.abs(np.array(local(S, t)[:5]) - np.array(local(S, 0)[:5])).max() for t in TS]
    print(f"  {nm}: max |local(t) - local(0)| over t in {TS}: {max(devs):.3e}; X_J(200) = {XJ(S, TLONG):.6f}; D(0) = {Dt(S, 0, T):.5f}")
print(f"  closed form for A: X_J(inf) = (Ts - Tb) + T r^2/2 = {0.5 * T * r_inf ** 2:.6f}")
DA = SA - SB
print("  what distinguishes A from B (Delta = A - B, all temperatures equal):")
print(f"    site-1 q-variance: A {SA[0,0]:.5f} (bare T/2.3) vs B {SB[0,0]:.5f} (dressed T r)  -> marginal is interaction-dressed in B")
print(f"    S-B cross q-covariance <q1 q2>: A {SA[0,1]:.5f} vs B {SB[0,1]:.5f} (= T r^2 = {T*r_inf**2:.5f}); |cross block| A {np.linalg.norm(SA[0,1:N]):.3f} vs B {np.linalg.norm(SB[0,1:N]):.5f}")
print(f"    bath q-block: |A_BB - B_BB| = {np.linalg.norm(DA[1:N,1:N]):.5f} (rank {np.linalg.matrix_rank(DA[1:N,1:N], tol=1e-9)}); p-sector identical: {np.abs(DA[N:,N:]).max():.1e}")
print(f"    interaction energy <E_int> = -<q1 q2>: A {-SA[0,1]:.5f} vs B {-SB[0,1]:.5f}; rank(A - B) = {np.linalg.matrix_rank(DA, tol=1e-9)}")
print("  -> TEMPERATURE DIFFERENCE NOT NECESSARY for the S6 transient. The mismatch is entirely interaction-induced")
print("     (missing S-B correlation + undressed marginals), all in the q-sector.")

# marginal-matched product: Gibbs marginals, zero cross -> isolates H_cross
SM = SB.copy(); SM[0, 1:N] = 0; SM[1:N, 0] = 0
print(f"  ISOLATING H_cross: product of the Gibbs MARGINALS (S1(0) = S_ref exactly, bath marginal = Gibbs), cross = 0:")
print(f"    min eig {mineig(SM):.3e}; D(0) = {Dt(SM, 0, T):.2e}; D(t) at t = 1, 2, 5, 10, 20: "
      + ", ".join(f"{Dt(SM, t, T):.2e}" for t in (1, 2, 5, 10, 20)))
print(f"    X_J(200) = {XJ(SM, TLONG):.6f} vs prediction T r^2 = {T * r_inf**2:.6f} (= Q12_G - Q12(0); marginals give 0)")
print("  -> with EVERY marginal at its equilibrium value, the missing cross-correlation alone drives a transient and a net")
print("     J transfer: CORRELATION / INTERACTION-ENERGY MISMATCH IS PHYSICALLY LOAD-BEARING.")

# =============================================================================================
hdr("B2-2 SAME MARGINALS, DIFFERENT S-B CORRELATIONS (S6 marginals fixed)")
def with_cross(Ts, Tb, alpha, beta=0.0, gamma=0.0):
    S = s6_product(Ts, Tb)
    S[0, 1:N] = S[1:N, 0] = alpha * Tb * Kinv[0, 1:N]                 # q1-q_B cross along the Gibbs direction
    S[N, N + 1] = S[N + 1, N] = beta * Tb                            # p1-p2 cross
    S[0, N + 1] = S[N + 1, 0] = gamma                                # q1-p2 (mixed) cross: J(0) = gamma
    return S
print(f"  admissibility of the q-cross family at S6 marginals (Ts = Tb = T): Schur complement T(1/2.3 - alpha^2 r^3) > 0 -> |alpha| < {np.sqrt((1/2.3) / r_inf**3):.4f}")
fwd = None
for (Ts, Tb) in [(1, 1), (2, 1)]:
    print(f"  marginals of the S6 member (Ts, Tb) = ({Ts}, {Tb}):")
    base = with_cross(Ts, Tb, 0.0)
    cases = [("1 zero cross (S6)", 0, 0, 0), ("2 weak q-cross alpha=0.3", 0.3, 0, 0), ("3 strong q-cross alpha=1.0 (Gibbs value)", 1.0, 0, 0),
             ("3' near-maximal alpha=1.4", 1.4, 0, 0), ("3'' p-cross beta=0.5", 0, 0.5, 0), ("3''' mixed q1-p2 gamma=0.3", 0, 0, 0.3)]
    for nm, a, b, g in cases:
        S = with_cross(Ts, Tb, a, b, g); me = mineig(S)
        if me <= 0: print(f"    {nm}: REJECTED (min eig {me:.2e})"); continue
        l0 = local(S, 0)[:3]; same = np.allclose(l0, local(base, 0)[:3])
        traj = np.array([local(S, t)[:3] for t in (2, 5, 10)]) - np.array([local(base, t)[:3] for t in (2, 5, 10)])
        xj = XJ(S, TLONG); xs = Dt(S, 0, Tb) - Dt(S, TLONG, Tb)
        Js = [local(S, t)[4] for t in (0.5, 1, 2)]
        print(f"    {nm}: min eig {me:.2e}; S1(0) same as S6: {same}; max |S1 - S1_S6| at t=2,5,10: {np.abs(traj).max(axis=1)};"
              f" J(0.5,1,2) = {np.round(Js, 4)}; X_J(inf) ~ {xj:+.5f}; X_sigma(inf) ~ {xs:.5f}")
    print(f"    closed form at fixed S6 marginals: X_J(inf) = (Ts-Tb) + Tb r^2/2 - Q12(0) = (Ts-Tb) + Tb r^2 (1/2 - alpha)")

# 4: cross block grafted from a time-evolved preparation, marginals restored to S6, scaled to admissibility
Ts, Tb = 2, 1; S0 = s6_product(Ts, Tb); tau = 3.0
Rq = rows_at(S0, tau, ["q1", "p1"])
for sc in (1.0, 0.7, 0.5, 0.3):
    S = S0.copy(); idx = [i for i in range(2 * N) if i not in (0, N)]
    for a, row in zip((0, N), Rq): S[a, idx] = S[idx, a] = sc * row[idx]
    me = mineig(S)
    if me > 0:
        xj = XJ(S, TLONG)
        print(f"  4 cross block grafted from the S6 (2,1) state evolved to t={tau}, marginals restored, scale {sc}: min eig {me:.2e};"
              f" X_J(inf) ~ {xj:+.5f} (S6: {XJ(S0, TLONG):+.5f}); max |S1 - S1_S6| at t = 2: {np.abs(np.array(local(S, 2)[:3]) - np.array(local(S0, 2)[:3])).max():.4f}")
        break
    print(f"  4 graft scale {sc}: rejected (min eig {me:.2e})")
print("  -> with every marginal fixed, correlation changes alter the transient, J(t) and the SIGN of X_J(inf); X_sigma(inf) = D(0)")
print("     is blind to them (it only sees S1(0) and S1(inf)). H_cross IS AN INDEPENDENT PREPARATION DATUM.")

# =============================================================================================
hdr("B2-3 is productness needed for return to equilibrium? correlated finite-rank / trace-class perturbations")
for nm, S in [("S6 product (2,1)", s6_product(2, 1)), ("q-cross alpha=1.4", with_cross(2, 1, 1.4)), ("mixed gamma=0.3 + p-cross 0.5", with_cross(2, 1, 0, 0.5, 0.3)),
              ("marginal-matched product", SM), ("graft (B2-2 item 4)", S)]:
    Tb = 1.0; dev = [np.abs(np.array(local(S, t)[:3]) - Tb * np.array([rN, 1, 0])).max() for t in (10, 40, 100, 200)]
    print(f"  {nm}: rank(Sigma0 - Sigma_G) = {np.linalg.matrix_rank(S - gibbs(Tb), tol=1e-9)}; |S1(t) - S_ref| at t=10,40,100,200: {np.array(dev)}")
print("  -> local return to equilibrium holds for correlated finite-rank deviations (numerical, finite N). The LS-1 argument")
print("     extends: see B2_HCORR_RESULT.md, Proposition B2-P1 (bridge proof sketch, not canonically verified).")

# =============================================================================================
hdr("B2-4 correlation reversal / anti-arrow (same D, Sigma, A)")
Ts, Tb = 2.0, 1.0; S0 = s6_product(Ts, Tb); tau = 8.0
Rm = np.diag(np.r_[np.ones(N), -np.ones(N)])
Eq, Ep = M_apply(tau, np.eye(N), np.zeros((N, N))); Fq, Fp = M_apply(tau, np.zeros((N, N)), np.eye(N))
Mt = np.block([[Eq, Fq], [Ep, Fp]])
Stau = Mt @ S0 @ Mt.T; Srev = Rm @ Stau @ Rm
print(f"  forward: S6 (2,1) evolved to tau={tau}; reversed state Sigma_R = R Sigma(tau) R (p -> -p) declared as new t=0; min eig {mineig(Srev):.2e}")
print(f"  entropy: global S(Sigma(tau)) - S(Sigma_R) = {entropy(Stau) - entropy(Srev):.1e}; retained-marginal entropies equal: "
      f"{abs(np.linalg.det(Stau[np.ix_([0,N],[0,N])]) - np.linalg.det(Srev[np.ix_([0,N],[0,N])])):.1e}")
for t in (0, 2, 4, 6, 8):
    print(f"   t={t}: D_reversed(t) = {Dt(Srev, t, Tb):.5f}   vs   D_forward(tau - t) = {Dt(S0, tau - t, Tb):.5f}")
print(f"  X_J over [0, tau]: forward {XJ(S0, tau):+.5f}, reversed {XJ(Srev, tau):+.5f}")
print("  -> the reversed preparation runs back to the less-equilibrated S6 state (D rises to D(0)); D + Sigma + A DO NOT SELECT")
print("     THE ARROW BOUNDARY. The two states have identical entropies (global and retained marginal) and differ only in the sign")
print("     of their q-p correlations.")

# =============================================================================================
hdr("B2-5 special epoch: the S6 product state at t = 0, evolved to +t and -t")
S0 = s6_product(2, 1)
for t in (0.5, 1, 2, 5, 10, 20):
    print(f"   t=+-{t}: C_SB(+t) = {CSB(S0, t):.5f}, C_SB(-t) = {CSB(S0, -t):.5f}; D(+t) = {Dt(S0, t, 1):.5f}, D(-t) = {Dt(S0, -t, 1):.5f};"
          f" J(+t) = {local(S0, t)[4]:+.5f}, J(-t) = {local(S0, -t)[4]:+.5f}")
print(f"  C_SB(0) = {CSB(S0, 0):.1e}; Gibbs cross norm (stationary) = {CSB(gibbs(1.0), 0):.5f}")
print("  time-reversal relation: Sigma0 = R Sigma0 R  =>  Sigma(-t) = R Sigma(t) R: C_SB and D even in t, J odd.")
print("  -> SPECIAL MOMENT SELECTED BY PREPARATION; ORIENTATION NOT SELECTED (Janus: relaxation toward both time directions).")

# =============================================================================================
hdr("B2-7 temperature parameters and the stationary-state class")
for (Ts, Tb) in [(2, 1), (10, 1), (0.5, 1), (0.1, 1), (1, 1), (1.05, 1), (1, 1.05)]:
    S = s6_product(Ts, Tb)
    print(f"  ({Ts},{Tb}): X_J(inf) ~ {XJ(S, TLONG):+.5f} (closed {(Ts - Tb) + 0.5 * Tb * r_inf ** 2:+.5f}); D(0) = {Dt(S, 0, Tb):.5f}; S1(inf) -> Tb diag(r,1)")
for T in (0.5, 1.0, 3.0):
    print(f"  Gibbs at T={T}: stationary (max drift {np.abs(np.array(local(gibbs(T), 7.0)[:5]) - np.array(local(gibbs(T), 0)[:5])).max():.1e}) -> every T is admissible")
F = (1 / lam) * (1 + 0.4 * np.cos(2.0 * lam))                          # generalized Gibbs ensemble Q = F(K), P = K F(K)
QG = V @ np.diag(F) @ V.T; SG = np.block([[QG, Z], [Z, K @ QG]])
print(f"  non-thermal GGE Q = F(K), P = K F(K): min eig {mineig(SG):.2e}; stationary drift {np.abs(np.array(local(SG, 7.0)[:5]) - np.array(local(SG, 0)[:5])).max():.1e};"
      f" retained marginal diag({QG[0,0]:.4f}, {(K@QG)[0,0]:.4f}) is not a Gibbs marginal at any single T")
print("  -> stationarity under D admits every T AND non-thermal GGEs (the harmonic chain is integrable): THERMAL STATE CLASS")
print("     CONDITIONALLY DEFINED (by a supplied Gibbs/KMS postulate); TEMPERATURE VALUE SUPPLIED.")

# =============================================================================================
hdr("B2-8 / B2-9 does the environment (bath state) carry H_corr? relaxation mechanism vs boundary selection")
Ts = 2.0
for nm, bathQ, bathP in [("bath Gibbs T_b=1", KBBinv, np.eye(N - 1)), ("bath Gibbs T_b=1.5", 1.5 * KBBinv, 1.5 * np.eye(N - 1))]:
    S = s6_product(Ts, 1.0); S[1:N, 1:N] = bathQ; S[N + 1:, N + 1:] = bathP
    print(f"  system Gibbs Ts=2 + {nm}: S1(200) = {np.array(local(S, TLONG)[:3])} -> the bath's state fixes the local limit")
QBg = QG[1:, 1:]; S = s6_product(Ts, 1.0); S[1:N, 1:N] = QBg; S[N + 1:, N + 1:] = (K @ QG)[1:, 1:]
print(f"  system Gibbs Ts=2 + bath in the restriction of the GGE: S1(200) = {np.array(local(S, TLONG)[:3])} vs GGE marginal diag({QG[0,0]:.4f}, {(K@QG)[0,0]:.4f})")
print("  same bath T_b = 1, different (H_sys, H_cross): the local limit is identical, the transients and X_J differ:")
for nm, S in [("S6 (2,1)", s6_product(2, 1)), ("S6 (2,1) + q-cross 1.4", with_cross(2, 1, 1.4)), ("S6 (0.1,1)", s6_product(0.1, 1)), ("marginal-matched", SM)]:
    print(f"    {nm}: S1(200) = {np.array(local(S, TLONG)[:3])}; X_J(inf) ~ {XJ(S, TLONG):+.5f}")
print("  -> the bath state determines the local asymptotic state CONDITIONAL on T_b (PARTIAL RELOCATION INTO ENVIRONMENT);")
print("     T_b, H_sys, H_cross and the epoch are not determined. D COMPRESSES MEMORY OF H_corr DOWNSTREAM (the local state")
print("     forgets it) but the integrated record X_J keeps it: X_J(inf) = E1(0) - E1_G + Q12_G - Q12(0).")

# =============================================================================================
hdr("B2-10 noise-layer firewall: overdamped stochastic environment x' = -K x + noise, Q = 2 diag(T_i) (L0-1e)")
from scipy.linalg import solve_continuous_lyapunov, expm
n = 40; Kn = 2.3 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
for nm, Tn in [("uniform T=1", np.ones(n)), ("non-uniform T_i", np.linspace(0.5, 1.5, n))]:
    Sst = solve_continuous_lyapunov(-Kn, -2 * np.diag(Tn))
    a = np.zeros((n, n)); a[0, 0] = Sst[0, 0]; a[1:, 1:] = Sst[1:, 1:]           # same marginals as stationary, zero cross
    b = Sst.copy()
    Mt = expm(-Kn * 3.0); Ma = Mt @ a @ Mt.T + (Sst - Mt @ Sst @ Mt.T); Mb = Mt @ b @ Mt.T + (Sst - Mt @ Sst @ Mt.T)
    print(f"  {nm}: stationary <x1 x2> = {Sst[0,1]:.5f} (cross-correlated; = T K^-1 for uniform: {Kinv[0,1] if nm.startswith('uniform') else float('nan'):.5f});"
          f" product-of-marginals prep vs stationary prep: |<x1x2> diff| at t=0 {abs(a[0,1]-b[0,1]):.4f}, t=3 {abs(Ma[0,1]-Mb[0,1]):.2e}")
print("  -> the noise layer fixes the STATIONARY correlation structure conditional on T_i (another environment state:")
print("     RELOCATION of the reference state into supplied T_i), and erases any initial H_cross at rate >= 2*0.3; it does not")
print("     select which initial correlation was prepared. Not identified with the S6 conservative bath.")

# =============================================================================================
hdr("B2-11 entropy firewall")
SA = s6_product(1, 1); SB = gibbs(1.0)
print(f"  equal-T factorized A vs Gibbs B: entropy S(A) - S(B) = {entropy(SA) - entropy(SB):+.5f}; A relaxes (X_J = +0.169), B is stationary")
print(f"  marginal-matched product (Gibbs marginals, zero cross) vs Gibbs: S(M) - S(B) = {entropy(SM) - entropy(SB):+.5f} (= the Gibbs S-B mutual"
      f" information); M is NOT stationary (X_J = +{XJ(SM, TLONG):.4f}), B is: the HIGHER-entropy state relaxes")
print("  forward-evolved vs reversed state (B2-4): identical global and marginal entropies, opposite local arrows.")
print("  -> neither 'low entropy' nor 'temperature gradient' labels the boundary: the load-bearing datum is the correlation")
print("     structure relative to the split (and, for reversal, the sign of q-p correlations).")
