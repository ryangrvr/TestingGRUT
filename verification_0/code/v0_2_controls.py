#!/usr/bin/env python3
"""V0-2 / P-15 blind reproduction: controls C1-C3 plus hostile side-checks.

Exact symbolic (sympy) wherever possible; mpmath quadrature at 40 digits otherwise.
No Monte Carlo.  Run:  python3 v0_2_controls.py > v0_2_controls.log
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
x, y, p, u = sp.symbols('x y p u', real=True)
D, kap, lam = sp.symbols('D kappa lambda', positive=True)


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


def feller(name, b, sig2, a_half=sp.Rational(1, 2), sdens=None):
    """Feller classification at the left endpoint 0 (by symmetry the same at 1 for all examples here).
    s'(x) = exp(-int_{1/2}^x 2b/sigma^2), m(x) = 2/(s' sigma^2) (speed density),
    S = int_0^{1/4} s', M = int_0^{1/4} m, Sigma = int_0^{1/4} S(0,x] m(x) dx, N = int_0^{1/4} M(0,x] s'(x) dx.
    Tested numerically by truncating at eps -> 0 (mpmath, symbols set to 1) and watching for divergence;
    symbolic s', m printed for the record."""
    sp_dens = sdens if sdens is not None else sp.simplify(sp.exp(-sp.integrate(sp.simplify(2 * b / sig2).subs(x, y), (y, a_half, x))))
    if sdens is not None:  # verify the supplied s' solves s'' = -(2b/sigma^2) s'
        assert sp.simplify(sp.diff(sp_dens, x) + 2 * b / sig2 * sp_dens) == 0
    m = sp.simplify(2 / (sp_dens * sig2))
    print(f"[{name}] s'(x) = {sp_dens},   m(x) = {m}")
    subsd = {D: 1, kap: 1, lam: 2}
    sf = sp.lambdify(x, sp_dens.subs(subsd), 'mpmath'); mf = sp.lambdify(x, m.subs(subsd), 'mpmath')
    q = mp.mpf(1) / 4
    rows = []
    for k in [4, 8, 16]:
        e = mp.mpf(10) ** (-k)
        S = mp.quad(sf, [e, q]); M = mp.quad(mf, [e, q])
        # Sigma = int_e^q S(e,x] m(x) dx = int_e^q s'(y) M[y,q] dy (Fubini) ; N = int_e^q m(y) S[y,q] dy
        Sig = mp.quad(lambda yy: sf(yy) * mp.quad(mf, [yy, q]), [e, mp.mpf(10)**(-k//2) if k>4 else e*10, q])
        N = mp.quad(lambda yy: mf(yy) * mp.quad(sf, [yy, q]), [e, mp.mpf(10)**(-k//2) if k>4 else e*10, q])
        rows.append((k, S, M, Sig, N))
        print(f"   eps=1e-{k:<2d}: S={mp.nstr(S,8):12s} M={mp.nstr(M,8):12s} Sigma={mp.nstr(Sig,8):12s} N={mp.nstr(N,8)}")
    fin = lambda i: abs(rows[-1][i] - rows[-2][i]) < 1e-3 * (1 + abs(rows[-1][i]))
    Sf, Sgf, Nf = fin(1), fin(3), fin(4)
    if Sgf and Nf: kind = "regular"
    elif Sgf: kind = "exit (accessible, attracting; finite-time absorption)"
    elif Nf: kind = "entrance (inaccessible, non-attracting)"
    else: kind = "natural, " + ("ATTRACTING (s(0+) finite): approached only as t->inf" if Sf else "non-attracting (s(0+) infinite)")
    print(f"   => boundary 0 is {kind}")
    return sp_dens, m, kind


# ---------------------------------------------------------------------------
hdr("C1(a)  driftless Wright-Fisher: sigma^2 = 2 D p(1-p), b = 0")
sig2 = 2 * D * x * (1 - x)
feller("WF", sp.Integer(0), sig2)
# h: scale s(x)=x -> h(p)=p.  Check BVP: (1/2)sigma^2 h'' + b h' = 0 with h=x
h = x
print("ODE residual for h=x:", sp.simplify(sp.Rational(1, 2) * sig2 * sp.diff(h, x, 2)))
# Expected exit time: (1/2) sigma^2 T'' = -1, T(0)=T(1)=0
T = -(x * sp.log(x) + (1 - x) * sp.log(1 - x)) / D
res = sp.simplify(sp.Rational(1, 2) * sig2 * sp.diff(T, x, 2) + 1)
print("T(p) = -[p ln p + (1-p) ln(1-p)]/D ; residual of (1/2)sigma^2 T'' + 1 :", res)
print("T(0+) =", sp.limit(T, x, 0, '+'), "  T(1-) =", sp.limit(T, x, 1, '-'))
# Independent check: Green's function representation E_p tau = int_0^1 G(p,y) m(y) dy,
# natural scale s(x)=x on (0,1):  G(p,y) = 2*min(p,y)*(1-max(p,y))/(s(1)-s(0)) with m(y)=1/(sigma^2(y)) ...
# Use E_p tau = int_0^1 G(p,y) * 2/sigma^2(y) dy with G(p,y) = min(p,y)(1-max(p,y)).
Dv = mp.mpf(1)
for pv in [mp.mpf('0.1'), mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.9')]:
    s2 = lambda yy: 2 * Dv * yy * (1 - yy)
    g = lambda yy: (yy * (1 - pv) if yy < pv else pv * (1 - yy)) * 2 / s2(yy)
    q = mp.quad(g, [0, pv, 1])
    closed = -(pv * mp.log(pv) + (1 - pv) * mp.log(1 - pv)) / Dv
    print(f"   D=1, p={float(pv):5.2f}: Green quadrature E tau = {mp.nstr(q, 20)}   closed form = {mp.nstr(closed, 20)}   diff = {mp.nstr(q - closed, 3)}")

# ---------------------------------------------------------------------------
hdr("C1(b)  variant sigma = 4 sqrt(kappa) p(1-p), b = 0")
sig2v = 16 * kap * x ** 2 * (1 - x) ** 2
feller("QSD-variant", sp.Integer(0), sig2v)
print("Scale s(x)=x is finite at both ends => both endpoints ATTRACTING, but Sigma=inf => NOT exit:")
print("natural attracting boundaries: p_t -> {0,1} a.s. as t->inf, never reached in finite time.")
# exit-time equation has no finite solution:
Tv_pp = -2 / sig2v
Tv = sp.integrate(sp.integrate(Tv_pp.subs(x, y), (y, sp.Rational(1, 2), u)), (u, sp.Rational(1, 2), x))
print("Particular sol. of (1/2)sigma^2 T''=-1 (T(1/2)=T'(1/2)=0):", sp.simplify(Tv))
print("  limit at 0+ :", sp.limit(sp.simplify(Tv), x, 0, '+'), " -> no bounded solution, E tau = +inf (indeed tau = inf a.s.)")
# Independent logit check: y = ln(p/(1-p)).  Ito: dy = y'(p) dp + (1/2) y''(p) sigma^2 dt
Y = sp.log(x / (1 - x))
drift_y = sp.simplify(sp.Rational(1, 2) * sp.diff(Y, x, 2) * sig2v)
diff_y = sp.simplify(sp.diff(Y, x) * sp.sqrt(sig2v))
print("logit y: drift =", drift_y, "  diffusion coeff =", diff_y)
yy = sp.symbols('yy', real=True)
drift_in_y = sp.simplify(drift_y.subs(x, 1 / (1 + sp.exp(-yy))))
print("   drift in terms of y:", sp.simplify(drift_in_y.rewrite(sp.tanh)), " (= 8 kappa tanh(y/2))")
sy_dens = sp.exp(-sp.integrate(2 * 8 * kap * sp.tanh(u / 2) / (16 * kap), (u, 0, yy)))
sy_dens = sp.simplify(sy_dens)
sy = sp.simplify(sp.integrate(sy_dens.subs(yy, u), (u, 0, yy)))
print("   scale of y-process: s_y'(y) =", sy_dens, " s_y(y) =", sy)
print("   s_y(+-inf) =", sp.limit(sy, yy, sp.oo), sp.limit(sy, yy, -sp.oo), " (finite: y -> +-inf a.s.)")
hy = sp.simplify((sy - sp.limit(sy, yy, -sp.oo)) / (sp.limit(sy, yy, sp.oo) - sp.limit(sy, yy, -sp.oo)))
pp_ = sp.symbols('pp_', positive=True)
hyp = sp.simplify(hy.rewrite(sp.exp).subs(yy, sp.log(pp_ / (1 - pp_))))
print("   P(y -> +inf) =", sp.simplify(hy.rewrite(sp.exp)), "; with y = logit(p):", hyp)
assert sp.simplify(hyp - pp_) == 0
print("   => h(p) = p  CONFIRMED (independent of the martingale argument)")

# ---------------------------------------------------------------------------
hdr("C2  b = lambda p(1-p)(2p-1), sigma^2 = 2 D p(1-p); frozen lambda=2, D=1")
b2 = lam * x * (1 - x) * (2 * x - 1)
sig2 = 2 * D * x * (1 - x)
ratio = sp.simplify(2 * b2 / sig2)
print("2b/sigma^2 =", ratio)
sdens = sp.simplify(sp.exp(-sp.integrate(ratio.subs(x, y), (y, sp.Rational(1, 2), x))))
print("s'(x) =", sdens)
s_sym = sp.integrate(sdens.subs(x, y), (y, 0, x))
h_sym = sp.simplify(s_sym / s_sym.subs(x, 1))
print("h(p) = s(p)/s(1) =", h_sym)
print("ODE residual b h' + (1/2) sigma^2 h'':", sp.simplify(b2 * sp.diff(h_sym, x) + sp.Rational(1, 2) * sig2 * sp.diff(h_sym, x, 2)))
print("h(0) =", sp.simplify(h_sym.subs(x, 0)), " h(1) =", sp.simplify(h_sym.subs(x, 1)))
feller("C2", b2.subs({lam: 2, D: 1}), sig2.subs(D, 1))
hnum = sp.lambdify(x, h_sym.subs({lam: 2, D: 1}), 'mpmath')
c = mp.sqrt(2)
h_closed = lambda pv: (mp.erf(c * (pv - mp.mpf(1) / 2)) + mp.erf(1 / c)) / (2 * mp.erf(1 / c))
h_quad = lambda pv: mp.quad(lambda t: mp.exp(-2 * (t - mp.mpf(1) / 2) ** 2), [0, pv]) / mp.quad(lambda t: mp.exp(-2 * (t - mp.mpf(1) / 2) ** 2), [0, 1])
spec = {'0.1': 0.07793, '0.25': 0.21955, '0.4': 0.38390, '0.5': 0.50000, '0.6': 0.61610, '0.75': 0.78045, '0.9': 0.92207}
print("\n  p     h (closed form erf, 20 d)      h (quadrature)            |diff|    spec      spec-ok(5dp)  h-p")
for ps, sv in spec.items():
    pv = mp.mpf(ps)
    a = h_closed(pv); q = h_quad(pv); s3 = hnum(pv)
    assert abs(a - s3) < mp.mpf('1e-30')
    ok = round(float(a), 5) == sv
    print(f"  {ps:4s}  {mp.nstr(a, 20):24s}  {mp.nstr(q, 20):24s}  {mp.nstr(abs(a-q), 2):8s}  {sv:.5f}  {ok}        {mp.nstr(a - pv, 8)}")
print("  closed form: h(p) = [erf(sqrt2 (p-1/2)) + erf(1/sqrt2)] / [2 erf(1/sqrt2)]")
print("  symmetry: h(p)+h(1-p) =", mp.nstr(h_closed(mp.mpf('0.1')) + h_closed(mp.mpf('0.9')), 25))

# ---------------------------------------------------------------------------
hdr("C3  decomposition counterexample under C2")
h02, h04, h06 = h_closed(mp.mpf('0.2')), h_closed(mp.mpf('0.4')), h_closed(mp.mpf('0.6'))
mix = (h02 + h06) / 2
print(f"  h(0.2) = {mp.nstr(h02, 15)}   h(0.6) = {mp.nstr(h06, 15)}")
print(f"  pure p=0.4           : {mp.nstr(h04, 15)}   (spec 0.38390)")
print(f"  50/50 mix of 0.2,0.6 : {mp.nstr(mix, 15)}   (spec 0.39271)")
print(f"  difference           : {mp.nstr(mix - h04, 10)}")
print("  Under C1 (h=p): pure 0.4, mixture (0.2+0.6)/2 = 0.4 exactly.")
print("\n  NOTE: pure |psi> with p=0.4 and the 50/50 mixture of pure p=0.2, p=0.6 states have DIFFERENT")
print("  density matrices (off-diagonals).  Same-rho comparisons (Bloch vector r=(0,0,-0.2), rho=diag(0.4,0.6)):")
z = lambda pv: 2 * pv - 1
# decomposition A: |1>,|0> with weights 0.4, 0.6 (p=1 and p=0)
decA = mp.mpf('0.4') * 1 + mp.mpf('0.6') * 0
# decomposition B: p=0.4 with phases 0, pi (weights 1/2 each)
decB = h04
# decomposition C: p=0.2 phases 0,pi and p=0.6 phases 0,pi, weights 1/4 each
decC = (h02 + h06) / 2
for name, pts in [("A {p=1,p=0} w=.4,.6", [(mp.mpf('0.4'), 1, 0), (mp.mpf('0.6'), 0, 0)]),
                  ("B {p=.4 phi=0, p=.4 phi=pi}", [(mp.mpf('0.5'), mp.mpf('0.4'), 0), (mp.mpf('0.5'), mp.mpf('0.4'), mp.pi)]),
                  ("C {.2/0,.2/pi,.6/0,.6/pi}", [(mp.mpf('0.25'), mp.mpf('0.2'), 0), (mp.mpf('0.25'), mp.mpf('0.2'), mp.pi),
                                                 (mp.mpf('0.25'), mp.mpf('0.6'), 0), (mp.mpf('0.25'), mp.mpf('0.6'), mp.pi)])]:
    rho00 = sum(w * pp for w, pp, ph in pts)
    rho01 = sum(w * mp.sqrt(pp * (1 - pp)) * mp.exp(1j * ph) for w, pp, ph in pts)
    val = sum(w * (h_closed(pp) if 0 < pp < 1 else pp) for w, pp, ph in pts)
    print(f"   {name:32s} rho_11={mp.nstr(rho00, 6)} |rho_01|={mp.nstr(abs(rho01), 3)}  outcome-1 freq = {mp.nstr(val, 10)}")
print("  => same rho, three different frequencies under C2: rho-dependence also violated.")

# ---------------------------------------------------------------------------
hdr("HOSTILE side-checks")
print("(H1) interior zero of sigma, b = 0:  sigma(p) = 4 sqrt(kappa) p(1-p)|2p-1|")
sig2i = 16 * kap * x ** 2 * (1 - x) ** 2 * (2 * x - 1) ** 2
fI = sp.lambdify(x, (1 / sig2i).subs(kap, 1), 'mpmath')
for k in [2, 4, 8]:
    print(f"   int_(1/4)^(1/2-1e-{k}) dx/sigma^2 = {mp.nstr(mp.quad(fI, [mp.mpf(1)/4, mp.mpf(1)/2 - mp.mpf(10)**(-k)]), 10)}")
print("   -> diverges like 1/(4 (1/2-x)): 1/sigma^2 not locally integrable at 1/2 => 1/2 is an absorbing (natural, attracting) trap")
print("   started at p<1/2: p_t bounded martingale on [0,1/2], converges to 0 or 1/2,")
print("   P(->1/2) = 2p, P(->1) = 0.  So h(p) = 0 on (0,1/2], h(p) = 2p-1 on [1/2,1): martingale but NOT Born.")
print("   e.g. h(0.4) = 0 vs Born 0.4;  50/50 mix of {0,1}: 0.5 vs h(0.5) = 0  -> decomposition-dependent.")

print("\n(H2) both endpoints non-attracting (entrance): sigma^2=2p(1-p), b=(1/2-p)*c with c=3 (b set to 0 at endpoints)")
sd = (4 * x * (1 - x)) ** sp.Rational(-3, 2)
feller("H2", sp.Rational(3) * (sp.Rational(1, 2) - x), 2 * x * (1 - x), sdens=sd)
print("   s(0+)=-inf, s(1-)=+inf: recurrent in (0,1); p_t never converges => h = 0 on (0,1), h(1)=1.")
print("   h is constant on (0,1): interior decompositions give p-bar-independent value 0, yet h != p and b != 0.")
print("   With decompositions allowed to use the endpoints: 0.5 h(0)+0.5 h(1) = 0.5 != h(0.5) = 0, detected.")

print("\n(H3) b nonzero on a Lebesgue-null set N (e.g. b = 1 on rationals): occupation-time formula")
print("   int_0^t 1_N(p_s) ds = int 1_N(a) L_t^a / sigma^2(a) da = 0, so the drift integral vanishes:")
print("   law identical to b=0, h = p, p_t martingale, yet b != 0 pointwise.  Only b=0 a.e. is determined.")
