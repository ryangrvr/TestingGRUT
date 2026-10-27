/* hom.c -- independent homogeneous X(2 real) + phi integrator (Claude review track pt_matching).
 *
 * Fields (conformal, canonical): x1, x2 (charged X, r^2 = x1^2 + x2^2), phi; momenta p1, p2, pphi.
 * Potential (per unit comoving volume, conformal units, radiation era so a''/a = 0):
 *   V = 1/2 m2(s) r^2 + 1/4 lamX(s) r^4 + 1/2 g2 phi^2 r^2 + 1/4 lamP phi^4 - 1/2 lamP M(s)^2 phi^2
 * mode 0 (matching ramps, phi frozen at 0): m2 = (mout f(s/T))^2, lamX = lamX0 (1 - f) if lamramp else lamX0.
 * mode 1 (phi condensation): m2 = 0, lamX = lamX0, M(s) = M0 exp(kappa s) (growth 0) or M0 (1 + kappa s) (growth 1).
 * Integrator: position-Verlet base map in extended phase space (s advanced in the drift, p_s kicked by
 * -dV/ds), symmetrically composed to order 6 (Yoshida 1990, 7 stages) or order 8 (Kahan-Li 1997, 15 stages).
 * The step h(s) = h0 / (wb0 + wbM * M(s)) depends on s only, so every step is a symplectic map.
 * y = {x1, x2, phi, p1, p2, pphi, s, ps}; H + ps is conserved up to integrator error (work check).
 * Written independently; does not use or import any Codex code.
 */
#include <math.h>
#include <string.h>

typedef struct {
  int mode, lamramp, ramp_type, growth;
  double lamX0, g2, lamP, T, mout, M0, kappa;
} Par;

/* optional linear wave test oscillators: chi_j'' = -(k_j^2 + MW2 + g2 phi^2 [+ wslave * g2 M^2 instead]) chi_j,
   stored in y[8 + 2j], y[9 + 2j]; they do not act back on the homogeneous fields. */
static int NW = 0; static double WK2[256]; static double MW2 = 0; static int WSLAVE = 0;
void set_waves(int nw, const double *k2, double mw2, int wslave) {
  int j; NW = nw; MW2 = mw2; WSLAVE = wslave; for (j = 0; j < nw; j++) WK2[j] = k2[j];
}

static double psi_fn(double x) { return x > 0 ? exp(-1.0 / x) : 0.0; }

/* ramp profile f(u) on [0,1] and df/du */
static double ramp_f(double u, int type, double *fp) {
  if (u <= 0) { *fp = 0; return 0; }
  if (u >= 1) { *fp = 0; return 1; }
  if (type == 0) { /* cosine ramp: f = (1 - cos(pi u))/2 */
    *fp = 0.5 * M_PI * sin(M_PI * u);
    return 0.5 * (1 - cos(M_PI * u));
  } else { /* C-infinity ramp: f = psi(u)/(psi(u)+psi(1-u)) */
    double a = psi_fn(u), b = psi_fn(1 - u), d = a + b;
    double ap = a / (u * u), bp = b / ((1 - u) * (1 - u));
    *fp = (ap * b + a * bp) / (d * d);
    return a / d;
  }
}

static double Mfun(const Par *P, double s, double *Mp) {
  if (P->growth == 0) { double M = P->M0 * exp(P->kappa * s); *Mp = P->kappa * M; return M; }
  *Mp = P->M0 * P->kappa; return P->M0 * (1 + P->kappa * s);
}

/* parameters at time s: m2, lamX, M2 and their s-derivatives */
static void params(const Par *P, double s, double *m2, double *dm2, double *lx, double *dlx, double *M2, double *dM2) {
  if (P->mode == 0) {
    double fp, f = ramp_f(s / P->T, P->ramp_type, &fp);
    fp /= P->T;
    *m2 = P->mout * P->mout * f * f; *dm2 = 2 * P->mout * P->mout * f * fp;
    if (P->lamramp) { *lx = P->lamX0 * (1 - f); *dlx = -P->lamX0 * fp; } else { *lx = P->lamX0; *dlx = 0; }
    *M2 = 0; *dM2 = 0;
  } else if (P->mode == 2) { /* control: X mass g M(s) imposed, no phi dynamics, no backreaction */
    double Mp, M = Mfun(P, s, &Mp);
    *m2 = P->g2 * M * M; *dm2 = 2 * P->g2 * M * Mp; *lx = P->lamX0; *dlx = 0; *M2 = M * M; *dM2 = 2 * M * Mp;
  } else {
    double Mp, M = Mfun(P, s, &Mp);
    *m2 = 0; *dm2 = 0; *lx = P->lamX0; *dlx = 0; *M2 = M * M; *dM2 = 2 * M * Mp;
  }
}

/* mode 3 (control): phi integrated out at its instantaneous minimum (Born-Oppenheimer, infinitely fast phi):
   V(u=r^2) = g2 M^2 u/2 + (lamX - g2^2/lamP) u^2/4 - lamP M^4/4  for g2 u < lamP M^2;  lamX u^2/4 otherwise. */
static void bo_terms(const Par *P, double u, double M2, double dM2, double *V, double *dVdu, double *dVds) {
  if (P->g2 * u < P->lamP * M2) {
    double le = P->lamX0 - P->g2 * P->g2 / P->lamP;
    *V = 0.5 * P->g2 * M2 * u + 0.25 * le * u * u - 0.25 * P->lamP * M2 * M2;
    *dVdu = 0.5 * P->g2 * M2 + 0.5 * le * u;
    *dVds = 0.5 * P->g2 * dM2 * u - 0.5 * P->lamP * M2 * dM2;
  } else {
    *V = 0.25 * P->lamX0 * u * u; *dVdu = 0.5 * P->lamX0 * u; *dVds = 0;
  }
}

static void kick(const Par *P, double *y, double h) {
  double m2, dm2, lx, dlx, M2, dM2;
  params(P, y[6], &m2, &dm2, &lx, &dlx, &M2, &dM2);
  if (P->mode == 3) {
    double V, dVdu, dVds, u = y[0] * y[0] + y[1] * y[1];
    bo_terms(P, u, M2, dM2, &V, &dVdu, &dVds);
    y[3] -= h * 2 * dVdu * y[0]; y[4] -= h * 2 * dVdu * y[1]; y[7] -= h * dVds;
    return;
  }
  double r2 = y[0] * y[0] + y[1] * y[1], f2 = y[2] * y[2];
  double cx = m2 + lx * r2 + P->g2 * f2;
  double cf = P->g2 * r2 + P->lamP * f2 - P->lamP * M2;
  y[3] -= h * cx * y[0];
  y[4] -= h * cx * y[1];
  if (P->mode == 1) y[5] -= h * cf * y[2];
  if (NW) {
    int j; double wm = WSLAVE ? P->g2 * M2 : P->g2 * f2;
    for (j = 0; j < NW; j++) y[9 + 2 * j] -= h * (WK2[j] + MW2 + wm) * y[8 + 2 * j];
  }
  double dVds = 0.5 * dm2 * r2 + 0.25 * dlx * r2 * r2 - 0.5 * P->lamP * dM2 * f2;
  y[7] -= h * dVds;
}

static void drift(double *y, double h) {
  int j;
  y[0] += h * y[3]; y[1] += h * y[4]; y[2] += h * y[5]; y[6] += h;
  for (j = 0; j < NW; j++) y[8 + 2 * j] += h * y[9 + 2 * j];
}

static void base(const Par *P, double *y, double h) { drift(y, 0.5 * h); kick(P, y, h); drift(y, 0.5 * h); }

static const double Y6[7] = {0.784513610477560, 0.235573213359357, -1.17767998417887, 1.31518632068391,
                             -1.17767998417887, 0.235573213359357, 0.784513610477560};
static const double KL8w[8] = {0.74167036435061295344822780, -0.40910082580003159399730010,
                               0.19075471029623837995387626, -0.57386247111608226665638773,
                               0.29906418130365592384446354, 0.33462491824529818378495798,
                               0.31529309239676659663205666, -0.79688793935291635401978884};

static void step(const Par *P, double *y, double h, int order) {
  int i;
  if (order == 2) { base(P, y, h); return; }
  if (order == 6) { for (i = 0; i < 7; i++) base(P, y, Y6[i] * h); return; }
  for (i = 0; i < 7; i++) base(P, y, KL8w[i] * h);
  base(P, y, KL8w[7] * h);
  for (i = 6; i >= 0; i--) base(P, y, KL8w[i] * h);
}

double energy(const Par *P, const double *y) {
  double m2, dm2, lx, dlx, M2, dM2;
  params(P, y[6], &m2, &dm2, &lx, &dlx, &M2, &dM2);
  if (P->mode == 3) {
    double V, dVdu, dVds, u = y[0] * y[0] + y[1] * y[1];
    bo_terms(P, u, M2, dM2, &V, &dVdu, &dVds);
    return 0.5 * (y[3] * y[3] + y[4] * y[4]) + V;
  }
  if (P->mode == 2) {
    double r2 = y[0] * y[0] + y[1] * y[1];
    return 0.5 * (y[3] * y[3] + y[4] * y[4]) + 0.5 * m2 * r2 + 0.25 * lx * r2 * r2;
  }
  double r2 = y[0] * y[0] + y[1] * y[1], f2 = y[2] * y[2];
  return 0.5 * (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) + 0.5 * m2 * r2 + 0.25 * lx * r2 * r2 +
         0.5 * P->g2 * f2 * r2 + 0.25 * P->lamP * f2 * f2 - 0.5 * P->lamP * M2 * f2;
}

double Mvalue(const Par *P, double s) { double Mp; return P->mode >= 1 ? Mfun(P, s, &Mp) : 0.0; }

/* Integrate y from y[6] to s_end. Step h = min(h0/(wb0 + wbM*M(s)), remaining).
 * If out != NULL, store a sample of the 8-vector y every time s passes the next multiple of dsample
 * (dsample <= 0: every step), up to nmax samples. Returns number of samples written (or steps if out==NULL). */
long run(const Par *P, double *y, double s_end, double h0, double wb0, double wbM, int order,
         double dsample, double *out, long nmax) {
  long ns = 0, nsteps = 0;
  double next = y[6] + (dsample > 0 ? dsample : 0);
  double tol = 1e-12 * (fabs(s_end) > 1 ? fabs(s_end) : 1);
  while (y[6] < s_end - tol) {
    double M = Mvalue(P, y[6]);
    double h = h0 / (wb0 + wbM * M);
    int last = 0;
    if (y[6] + h >= s_end - tol) { h = s_end - y[6]; last = 1; }
    step(P, y, h, order);
    if (last) y[6] = s_end; /* remove round-off in the clock */
    nsteps++;
    if (out && ns < nmax && (dsample <= 0 || y[6] >= next)) {
      memcpy(out + 8 * ns, y, 8 * sizeof(double));
      ns++;
      if (dsample > 0) while (next <= y[6]) next += dsample;
    }
  }
  return out ? ns : nsteps;
}

/* Fixed-step integration with exactly n steps of size h (for periodic sampling). Stores y after each step. */
void run_fixed(const Par *P, double *y, double h, long n, int order, double *out) {
  long i;
  for (i = 0; i < n; i++) {
    step(P, y, h, order);
    if (out) memcpy(out + 8 * i, y, 8 * sizeof(double));
  }
}
