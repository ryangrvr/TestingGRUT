/* sk_sym.c -- skeptic of track pt_m1tilt (Claude review, issue #13).  Independent integrator for the homogeneous
 * classical X + phi system under Codex's GRUT-M1 potential (d4747cb), written from scratch:
 *   NOT the track's adaptive DOP853 (m1_dop.c).  A fixed-order symplectic composition (Kahan-Li 17-stage order 8,
 *   Yoshida 7-stage order 6 as a check) of Strang (drift-kick-drift) steps in the EXTENDED phase space (s is a
 *   coordinate with conjugate momentum p_s, H_ext = p^2/2 + V(q, s) + p_s), so the explicit time dependence M(s) is
 *   handled exactly and the extended energy E + p_s is the conserved check quantity: W = -(p_s - p_s0) is the work
 *   integral and E(s) - E(0) - W(s) = Delta H_ext.  The step h depends on s only (h = 2 pi / (nper Omega_max(M(s)))),
 *   so every step is a symplectic map.
 * Potential (conformal variables chi = a phi, M = a v):
 *   V = lamX/4 r^4 + g2/2 chi^2 r^2 + lamP/4 (chi^2 - M^2)^2 - epst M (chi - M)^3
 * Derived here (checked against finite differences of V in the Python driver):
 *   dV/dx_i  = (lamX r^2 + g2 chi^2) x_i
 *   dV/dchi  = g2 r^2 chi + lamP (chi^2 - M^2) chi - 3 epst M (chi - M)^2
 *   dV/dM    = -lamP M (chi^2 - M^2) - epst (chi - M)^3 + 3 epst M (chi - M)^2 ;  dV/ds = M'(s) dV/dM
 * M(s): growth 0: M0 exp(kappa s); growth 1: M0 (1 + kappa s).
 * State y[8] = {x1, x2, chi, p1, p2, pchi, s, ps}. */
#include <math.h>

typedef struct { double lamX, g2, lamP, epst, M0, kappa, om_fix, fm; int growth, method; } SkPar;

static const double KL8[17] = {
  0.13020248308889008087881763, 0.56116298177510838456196441, -0.38947496264484728640807860,
  0.15884190655515560089621075, -0.39590389413323757733623154, 0.18453964097831570709183254,
  0.25837438768632204729397911, 0.29501172360931029887096624, -0.60550853383003451169892108,
  0.29501172360931029887096624, 0.25837438768632204729397911, 0.18453964097831570709183254,
  -0.39590389413323757733623154, 0.15884190655515560089621075, -0.38947496264484728640807860,
  0.56116298177510838456196441, 0.13020248308889008087881763 };
/* Yoshida 1990 order 6, solution A (7 stages) */
static const double Y6W[3] = { -1.17767998417887, 0.235573213359357, 0.784513610477560 };

static inline void Mfun(const SkPar *P, double s, double *M, double *Mp) {
  if (P->growth == 0) { *M = P->M0 * exp(P->kappa * s); *Mp = P->kappa * (*M); }
  else { *M = P->M0 * (1.0 + P->kappa * s); *Mp = P->M0 * P->kappa; }
}

double sk_V(const SkPar *P, const double *q, double s) {
  double M, Mp; Mfun(P, s, &M, &Mp);
  double r2 = q[0] * q[0] + q[1] * q[1], c = q[2], u = c * c - M * M, d = c - M;
  return 0.25 * P->lamX * r2 * r2 + 0.5 * P->g2 * c * c * r2 + 0.25 * P->lamP * u * u - P->epst * M * d * d * d;
}

/* gradient wrt (x1, x2, chi) and explicit dV/ds */
static inline void sk_grad(const SkPar *P, const double *q, double s, double *G, double *dVds) {
  double M, Mp; Mfun(P, s, &M, &Mp);
  double r2 = q[0] * q[0] + q[1] * q[1], c = q[2], u = c * c - M * M, d = c - M;
  double cx = P->lamX * r2 + P->g2 * c * c;
  G[0] = cx * q[0]; G[1] = cx * q[1];
  G[2] = P->g2 * r2 * c + P->lamP * u * c - 3.0 * P->epst * M * d * d;
  *dVds = Mp * (-P->lamP * M * u - P->epst * d * d * d + 3.0 * P->epst * M * d * d);
}

void sk_grad_py(const SkPar *P, const double *q, double s, double *G, double *dVds) { sk_grad(P, q, s, G, dVds); }

double sk_E(const SkPar *P, const double *y) {
  return 0.5 * (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) + sk_V(P, y, y[6]);
}

static inline void drift(double *y, double tau) {
  y[0] += tau * y[3]; y[1] += tau * y[4]; y[2] += tau * y[5]; y[6] += tau;
}
static inline void kick(const SkPar *P, double *y, double tau) {
  double G[3], dVds; sk_grad(P, y, y[6], G, &dVds);
  y[3] -= tau * G[0]; y[4] -= tau * G[1]; y[5] -= tau * G[2]; y[7] -= tau * dVds;
}

/* one composition step of size h: merged drifts */
static void step(const SkPar *P, double *y, double h) {
  const double *g; int n;
  double Y6[7];
  if (P->method == 8) { g = KL8; n = 17; }
  else {
    double w0 = 1.0 - 2.0 * (Y6W[0] + Y6W[1] + Y6W[2]);
    Y6[0] = Y6W[2]; Y6[1] = Y6W[1]; Y6[2] = Y6W[0]; Y6[3] = w0; Y6[4] = Y6W[0]; Y6[5] = Y6W[1]; Y6[6] = Y6W[2];
    g = Y6; n = 7;
  }
  drift(y, 0.5 * g[0] * h);
  for (int i = 0; i < n; i++) {
    kick(P, y, g[i] * h);
    if (i + 1 < n) drift(y, 0.5 * (g[i] + g[i + 1]) * h);
  }
  drift(y, 0.5 * g[n - 1] * h);
}

static inline double hstep(const SkPar *P, double s, int nper) {
  double M, Mp; Mfun(P, s, &M, &Mp);
  double om = sqrt(P->om_fix * P->om_fix + P->fm * 2.0 * P->lamP * M * M);
  return 2.0 * M_PI / (nper * om);
}

/* integrate from y[6] to s1; return number of steps */
long sk_run(const SkPar *P, double *y, double s1, int nper) {
  long n = 0;
  while (y[6] < s1) {
    double h = hstep(P, y[6], nper);
    if (y[6] + h >= s1) { h = s1 - y[6]; step(P, y, h); y[6] = s1; n++; break; }
    step(P, y, h); n++;
  }
  return n;
}

/* integrate and sample: ns samples at s = s0 + k ds (k = 1..ns) into out[ns][8] */
long sk_sample(const SkPar *P, double *y, double ds, int ns, int nper, double *out) {
  long n = 0;
  for (int k = 0; k < ns; k++) {
    double target = y[6] + ds;
    n += sk_run(P, y, target, nper);
    for (int c = 0; c < 8; c++) out[k * 8 + c] = y[c];
  }
  return n;
}
