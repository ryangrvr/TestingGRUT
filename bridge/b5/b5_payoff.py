"""B5 -- prediction / payoff firewall.  Every candidate relation is attacked by varying the still-supplied inputs and by a
standard (non-GRUT) comparator.  Independent re-implementation; canonical GRUT is read-only evidence.
Harmonic chain family: site 1 pin kappa, coupling g to the bath, bath = half-line chain with on-site pin c and unit hops
(S6: kappa = c = 2.3, g = 1).  Finite N = 600, times << recurrence; finite-N values are cross-checks."""
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(20261005)
def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
N = 600; TL = 200.0

def chain(kappa=2.3, g=1.0, c=2.3, n=N):
    K = c * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1); K[0, 0] = kappa; K[0, 1] = K[1, 0] = -g; return K
class Sys:
    def __init__(s, K):
        s.K = K; s.lam, s.V = np.linalg.eigh(K); s.om = np.sqrt(s.lam); s.Kinv = s.V @ np.diag(1 / s.lam) @ s.V.T; s.n = len(K)
    def M(s, t, vq, vp, tr=False):
        c = np.cos(s.om * t)[:, None]; sn = (np.sin(s.om * t) / s.om)[:, None]; ks = (s.om * np.sin(s.om * t))[:, None]
        aq, ap = s.V.T @ vq, s.V.T @ vp
        if not tr: return s.V @ (c * aq + sn * ap), s.V @ (-ks * aq + c * ap)
        return s.V @ (c * aq - ks * ap), s.V @ (sn * aq + c * ap)
    def rows(s, Sig, t, idx):
        n = s.n; E = np.zeros((2 * n, len(idx)))
        for k, i in enumerate(idx): E[i, k] = 1
        Rq, Rp = s.M(t, E[:n], E[n:], tr=True); W = Sig @ np.vstack([Rq, Rp]); Wq, Wp = s.M(t, W[:n], W[n:])
        return np.vstack([Wq, Wp]).T
    def local(s, Sig, t):
        R = s.rows(Sig, t, [0, s.n]); return R[0, 0], R[1, s.n], R[0, 1]          # Q11, P11, Q12
def gibbs(S, T): return np.block([[T * S.Kinv, np.zeros((S.n, S.n))], [np.zeros((S.n, S.n)), T * np.eye(S.n)]])
def product(S, Ts, Tb, alpha=0.0):
    """S6-type product: bare site-1 Gibbs (Ts/kappa, Ts), uncoupled-bath Gibbs at Tb; optional q1-qB cross alpha * (Gibbs cross at Tb)"""
    n = S.n; K = S.K; Q = np.zeros((n, n)); Q[0, 0] = Ts / K[0, 0]; Q[1:, 1:] = Tb * np.linalg.inv(K[1:, 1:])
    Q[0, 1:] = Q[1:, 0] = alpha * Tb * S.Kinv[0, 1:]
    P = Tb * np.eye(n); P[0, 0] = Ts
    return np.block([[Q, np.zeros((n, n))], [np.zeros((n, n)), P]])
def XJ(S, Sig, T):
    """bath self-energy gained = E1(0) - E1(T) + g(Q12(T) - Q12(0))  [energy bookkeeping, LS-2.1 generalized]"""
    k, g = S.K[0, 0], -S.K[0, 1]; a, b = S.local(Sig, 0.0), S.local(Sig, T)
    E1 = lambda l: 0.5 * k * l[0] + 0.5 * l[1]
    return E1(a) - E1(b) + g * (b[2] - a[2])
def bound_states(K, c):
    lam = np.linalg.eigvalsh(K); return [l for l in lam if l < c - 2 - 1e-3 or l > c + 2 + 1e-3]

# =============================================================================================
hdr("P1  R1 frame-overdetermination: frame(net) = frame(drift) = frame(noise)?")
n = 24; Kp = 0.3 * np.eye(n) + (2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)); Kp[0, 0] -= 1; Kp[-1, -1] -= 1
def drift_axes(O, beta=1.0):
    """axes of the quartic part V4 = beta * sum_i (o_i . x)^4, read from its Hessian at a random probe state (Jennrich-type
    recovery; observable as the nonlinear part of the force response)"""
    u = rng.normal(size=n); H4 = 12 * beta * O @ np.diag((O.T @ u) ** 2) @ O.T; return np.linalg.eigh(H4)[1]
def noise_axes(O, T):
    return np.linalg.eigh(O @ np.diag(T) @ O.T)[1]
def frame_overlap(A, B): return np.min(np.max(np.abs(A.T @ B), axis=1))          # 1 = same frame up to sign / permutation
I_ = np.eye(n); Ti = np.linspace(0.5, 1.5, n); Orot = np.linalg.qr(rng.normal(size=(n, n)))[0]
print(f"  A  common supplied site frame: overlap(net frame, drift axes) = {frame_overlap(I_, drift_axes(I_)):.6f};"
      f" overlap(net frame, noise axes) = {frame_overlap(I_, noise_axes(I_, Ti)):.6f}")
dA = drift_axes(Orot)
print(f"  B  drift written on-site in an independently rotated frame (admissible: V convex, gradient flow well posed):")
print(f"     overlap(net frame, drift axes) = {frame_overlap(I_, dA):.4f}; overlap(rotated frame, drift axes) = {frame_overlap(Orot, dA):.6f}")
x = rng.normal(size=n); x0n = np.linalg.norm(x); dt = 1e-3
for _ in range(20000): x = x - dt * (Kp @ x + 4 * Orot @ ((Orot.T @ x) ** 3))
print(f"     rotated-drift model relaxes: |x(20)| = {np.linalg.norm(x):.2e} (from |x(0)| = {x0n:.2f})")
print("     -> nothing earned forbids the independent rotation; agreement exists only because the declarations share a frame.")
print("  C  comparator: an ordinary nonlinear Langevin network x' = -K x - grad V4 + xi with site-local V4 and site-local noise")
print("     has the same frame agreement by construction (identical equations). -> CONSISTENCY-ONLY / STANDARD-STRUCTURE")

# =============================================================================================
hdr("P2  generalized integrated transfer X(inf) = E1(0) - E1_G + g Q12_G - g Q12(0): provenance and variation")
print("  provenance: d(E1 + E_int)/dt = -g q1 p2 (equations of motion) => X(T) = E1(0) - E1(T) + g(Q12(T) - Q12(0)) for ANY")
print("  quadratic chain and ANY state (energy accounting); the T -> inf value additionally needs local return to equilibrium.")
cases = [("S6 canonical (2.3,1,2.3)", 2.3, 1.0, 2.3), ("kappa=3.0", 3.0, 1.0, 2.3), ("g=0.6", 2.3, 0.6, 2.3),
         ("kappa=2.6,g=1.2,c=2.5", 2.6, 1.2, 2.5), ("bath pin c=2.8", 2.3, 1.0, 2.8)]
preps = [(2.0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 1.0, 1.0), (0.5, 1.5, 0.4)]
for nm, k, g, c in cases:
    S = Sys(chain(k, g, c)); bs = bound_states(S.K, c)
    out = []
    for (Ts, Tb, al) in preps:
        Sig = product(S, Ts, Tb, al); xnum = XJ(S, Sig, TL)
        Q11, P11, Q12 = S.local(Sig, 0)
        xth = (0.5 * k * Q11 + 0.5 * P11) - (0.5 * k * Tb * S.Kinv[0, 0] + 0.5 * Tb) + g * Tb * S.Kinv[0, 1] - g * Q12
        out.append(f"({Ts},{Tb},a={al}): {xnum:+.5f}/{xth:+.5f}")
    print(f"  {nm}: bound states {len(bs)}; X num/closed: " + "; ".join(out))
S = Sys(chain(1.0, 1.0, 2.3)); bs = bound_states(S.K, 2.3); Sig = product(S, 2.0, 1.0)
xs = [XJ(S, Sig, T) for T in np.linspace(150, 200, 11)]
print(f"  class boundary: kappa=1.0 (bound state(s) {np.round(bs, 4)}): X(T) over T in [150,200] spans [{min(xs):+.4f}, {max(xs):+.4f}]"
      f" -> no limit; the relation needs the a.c. relaxation class")
# comparator: arbitrary random harmonic network, multi-site system, finite time -- the identity is energy conservation
m = 30; W = np.triu(rng.uniform(0, 1, (m, m)) * (rng.uniform(size=(m, m)) < 0.25), 1); W = W + W.T
Kr = np.diag(W.sum(1) + rng.uniform(0.2, 0.6, m)) - W; Sset = [0, 1, 2]; Bset = list(range(3, m))
A = np.block([[np.zeros((m, m)), np.eye(m)], [-Kr, np.zeros((m, m))]])
Sig0 = np.block([[np.diag(rng.uniform(0.5, 2, m)), np.zeros((m, m))], [np.zeros((m, m)), np.diag(rng.uniform(0.5, 2, m))]])
def EB(Sg):
    Q, P = Sg[:m, :m], Sg[m:, m:]; return 0.5 * np.trace(P[np.ix_(Bset, Bset)]) + 0.5 * np.trace(Kr[np.ix_(Bset, Bset)] @ Q[np.ix_(Bset, Bset)])
def flux(Sg):  # dE_B/dt = -<p_B^T K_BS q_S>
    return -np.trace(Kr[np.ix_(Bset, Sset)] @ Sg[np.ix_(Sset, [m + b for b in Bset])])
ts = np.linspace(0, 6, 3001); Ms = [expm(A * t) for t in ts[::1]]
fl = np.array([flux(Mt @ Sig0 @ Mt.T) for Mt in Ms]); integ = (np.trapezoid if hasattr(np, 'trapezoid') else np.trapz)(fl, ts)
print(f"  comparator (random 30-node harmonic network, 3-site system): integral of bath flux over [0,6] = {integ:.6f};"
      f" E_B(6) - E_B(0) = {EB(Ms[-1] @ Sig0 @ Ms[-1].T) - EB(Sig0):.6f} -> the same identity holds in any harmonic network")

# =============================================================================================
hdr("P3  equal-temperature offset X(inf) = (1/2) T_b r^2 (S6): which parts are invariant?")
print("  general closed form (equal T, S6-type product with q1-qB cross alpha * Gibbs cross):")
print("     X(inf) = g Q12_G (1/2 - alpha) = -<E_int>_G (1/2 - alpha);  S6: g = 1, Q12_G = T r^2")
for nm, k, g, c in cases:
    S = Sys(chain(k, g, c)); T = 1.3
    for al in (0.0, 1.0):
        Sig = product(S, T, T, al); x = XJ(S, Sig, TL); EintG = -g * T * S.Kinv[0, 1]
        print(f"  {nm}, T={T}, alpha={al}: X = {x:+.5f}; -<E_int>_G (1/2 - alpha) = {-EintG * (0.5 - al):+.5f};"
              f" 'r^2' analogue (K^-1)_12 = {S.Kinv[0,1]:.5f}; ratio X/|<E_int>_G| = {x / abs(EintG):+.4f}")
print("  static identity in ANY quadratic network (multi-site system S, equal T, product of the two Gibbs marginals):")
for trial in range(3):
    m = 25; W = np.triu(rng.uniform(0, 1, (m, m)) * (rng.uniform(size=(m, m)) < 0.3), 1); W = W + W.T
    Kr = np.diag(W.sum(1) + rng.uniform(0.2, 0.6, m)) - W; Ss = list(range(4)); Bs = list(range(4, m)); T = rng.uniform(0.5, 2)
    Kinv = np.linalg.inv(Kr); ES0 = 0.5 * T * len(Ss) + 0.5 * T * len(Ss)                      # <E_S> in the product of Gibbs(K_SS), equipartition
    ESG = 0.5 * T * len(Ss) + 0.5 * T * np.trace(Kr[np.ix_(Ss, Ss)] @ Kinv[np.ix_(Ss, Ss)])
    EintG = T * np.trace(Kr[np.ix_(Ss, Bs)] @ Kinv[np.ix_(Bs, Ss)])
    print(f"   random network {trial}: E_S(prod) - E_S(Gibbs) = {ES0 - ESG:+.6f};  (1/2)<E_int>_G = {0.5 * EintG:+.6f}")
print("  -> the 1/2 is the quadratic-network equipartition identity (STANDARD); r is model-specific (kappa, g, c);")
print("     the sign flips for alpha > 1/2 (INPUT-DEPENDENT).")

# =============================================================================================
hdr("P6  correlation-boundary sign: sweep of all PSD-admissible q1-q2 correlations at fixed S6 marginals")
S = Sys(chain()); r = (2.3 - np.sqrt(1.29)) / 2
for (Ts, Tb) in [(1.0, 1.0), (2.0, 1.0), (10.0, 1.0), (0.5, 1.0)]:
    Q11 = Ts / 2.3; QBB = Tb * np.linalg.inv(S.K[1:, 1:]); Q22 = QBB[0, 0]
    qmax = np.sqrt(Q11 * Q22)                                          # Cauchy-Schwarz / PSD bound on Q12(0)
    vals = []
    for f in (-0.999, -0.5, 0.0, 0.5, 0.999):
        Sig = product(S, Ts, Tb); lamb = f * qmax / Q22; X = lamb * QBB[0, :]                  # extremal rank-1 cross along Q_BB e_2
        Sig[0, 1:S.n] = Sig[1:S.n, 0] = X; me = np.linalg.eigvalsh(Sig).min()
        vals.append((Sig[0, 1], XJ(S, Sig, TL), me))
    lo, hi = (Ts - Tb) + 0.5 * Tb * r ** 2 - qmax, (Ts - Tb) + 0.5 * Tb * r ** 2 + qmax
    print(f"  (Ts,Tb)=({Ts},{Tb}): admissible Q12(0) in [-{qmax:.4f}, {qmax:.4f}] -> X(inf) in [{lo:+.4f}, {hi:+.4f}];"
          f" sampled (Q12, X, min eig): " + "; ".join(f"({a:+.3f},{b:+.4f},{c_:.1e})" for a, b, c_ in vals))
print("  -> at equal T both signs occur: no parameter-free arrow sign (BOUNDARY-STATE-DEPENDENT). For large |Ts-Tb| the sign")
print("     is fixed, but only by a Cauchy-Schwarz / PSD bound with supplied marginals (generic positivity: not distinctive).")

# =============================================================================================
hdr("B5-9 dimensionless combinations suggested by the exact equations (variation firewall applied)")
for nm, k, g, c in cases[:3]:
    S = Sys(chain(k, g, c))
    for (Ts, Tb, al) in [(1, 1, 0.0), (1, 1, 0.7), (2, 1, 0.0)]:
        Sig = product(S, Ts, Tb, al); x = XJ(S, Sig, TL); q = g * Tb * S.Kinv[0, 1]
        rt = f"{x / (Ts - Tb):+.4f}" if Ts != Tb else "undef"
        print(f"  {nm} ({Ts},{Tb},a={al}): X/(g Q12_G) = {x / q:+.4f}; X/(Ts-Tb) = {rt}")
print("  -> X/(g Q12_G) = 1/2 only at alpha = 0 and Ts = Tb (killed by H_cross / H_marginal variation); X/(Ts-Tb) depends on K and alpha.")
