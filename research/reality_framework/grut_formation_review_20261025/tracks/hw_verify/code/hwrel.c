/* hw_verify track (Claude, own code; issue #13 review of ls_scaling / Codex 2908e91).
 *
 * Homogeneous X, Y, Sigma release in CONFORMAL time with CONFORMAL fields, couplings at run time.  This is a
 * different formulation and a different integrator from ls_scaling's scalrel.c (cosmic time, physical fields,
 * Dormand-Prince 5(4)).  Units H_i = 1.
 *
 *   chi_i = a z_i,  ' = d/deta,  a' = a^2 H.   Own derivation from z'' + 3 H z' + dV/dz* = 0:
 *   chi'' = (a''/a) chi - a^3 F(chi/a),   a''/a = a^2 rho_phi / 6   (0 in radiation).
 *   V = lx|X|^4 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + M^2|Y|^2 + 2Re(kap Y* X^3) + 2Re(mu S* Y^2) - cX |X|^2,
 *   cX = cH rho_phi/3 while the inflaton is present (background 1), 0 after an instant release (background 0).
 *
 * Backgrounds (fields do not back-react, as in every track):
 *   0: radiation from a = 1, a = eta, H = 1/a^2 (instant release and instant reheating);
 *   1: inflaton matter with rho_phi(1) = 3, decaying at width G into radiation; rho_r(1) = 0.
 *      Finite release duration ~ 1/G; G -> infinity reproduces background 0.
 * Initial state: chi(1) = z0, chi'(1) = z0 + v0 (v0 = physical dz/dt at release, default 0).
 *
 * Integrators (own):
 *   8: Yoshida (1990) 8th-order symmetric composition (solution A) of drift-kick-drift, extended phase space
 *      (eta and the background advance in the drift; background 1 by RK4 substeps).  Step c/Omega, Omega a
 *      Gershgorin bound at each species' quartic turning amplitude.  Gauss's law exact to rounding.
 *   5: Dormand-Prince 5(4) adaptive on the same conformal system (internal cross-check only).
 *
 * Diagnostics (all per comoving volume, H_i units):
 *   N_i = -2 Im(chi_i* chi_i');  E_X = |w_X|^2 + lx |chi_X|^4 with w = chi' - (a'/a) chi (= a^2 dz/dt);
 *   N_circ = 2 sqrt(2 lx) (E_X/(3 lx))^(3/4);  eps = N_X/N_circ;
 *   Y excitation about the forced response y_p (A y_p + B y_p* = -J, A = a^2 M^2, B = 2 a mu chi_S, J = kap chi_X^3)
 *   split into the two real Majorana modes in the frame co-rotating with arg(Sigma)/2 (velocity includes the
 *   -i theta'/2 term):  N_YH = E_H/omega_H, N_YL = E_L/omega_L, omega^2 = A +- |B| (conformal);
 *   Sigma radial quanta N_sig = a^3 [rdot^2 + ls (r^2 - h^2/2)^2]/m_sig, m_sig = sqrt(2 ls) h, r = |chi_S|/a;
 *   Gauss G = N_X + 3 N_Y + 6 N_S (should stay at its initial value, 0 for zero velocities).
 *
 * Usage:
 *  hwrel method cpar bg G cH lx ly ls kap M h mu a_end nsamp nwin w_1 ... w_nwin  x0r x0i y0r y0i s0r s0i
 *        [vx0r vx0i vy0r vy0i vs0r vs0i]
 * Output: S lines (samples), W lines (window means, trapezoid in a), F (final chi, chi'), END counters.
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double complex cplx;
static double lx, ly, ls, kap, MM, hh, mu, CH, GAM;
static int BG;

/* background state: eta, a, Rphi = rho_phi a^3, Rr = rho_r a^4 */
typedef struct { double eta, a, Rphi, Rr; } bgs;

static double Hof(const bgs *b) {
  if (BG == 0) return 1.0 / (b->a * b->a);
  double a = b->a;
  return sqrt((b->Rphi / (a * a * a) + b->Rr / (a * a * a * a)) / 3.0);
}
static double rhophi(const bgs *b) { return BG == 0 ? 0.0 : b->Rphi / (b->a * b->a * b->a); }
static double aprime(const bgs *b) { return b->a * b->a * Hof(b); }

static void bg_rhs(const bgs *b, double *d) { /* d/deta of (a, Rphi, Rr) */
  double a = b->a;
  d[0] = aprime(b); d[1] = -GAM * a * b->Rphi; d[2] = GAM * a * a * b->Rphi;
}
static void bg_advance(bgs *b, double h) {
  if (BG == 0) { b->eta += h; b->a = b->eta; return; }
  double k1[3], k2[3], k3[3], k4[3]; bgs t;
  bg_rhs(b, k1);
  t = *b; t.a += 0.5 * h * k1[0]; t.Rphi += 0.5 * h * k1[1]; t.Rr += 0.5 * h * k1[2]; bg_rhs(&t, k2);
  t = *b; t.a += 0.5 * h * k2[0]; t.Rphi += 0.5 * h * k2[1]; t.Rr += 0.5 * h * k2[2]; bg_rhs(&t, k3);
  t = *b; t.a += h * k3[0]; t.Rphi += h * k3[1]; t.Rr += h * k3[2]; bg_rhs(&t, k4);
  b->a += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
  b->Rphi += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
  b->Rr += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
  b->eta += h;
}

/* conformal force: chi'' = -Fc */
static void force(const cplx *c, const bgs *b, cplx *F) {
  double a = b->a;
  cplx x = c[0], y = c[1], s = c[2];
  double x2 = creal(x * conj(x)), y2 = creal(y * conj(y)), s2 = creal(s * conj(s));
  double A = a * a * MM * MM;
  double cX = CH * rhophi(b) / 3.0;              /* Hubble-induced X mass^2 (physical), background 1 only */
  double app = (BG == 0) ? 0.0 : a * a * rhophi(b) / 6.0;   /* a''/a */
  F[0] = 2 * lx * x2 * x + 3 * kap * y * conj(x) * conj(x) - a * a * cX * x - app * x;
  F[1] = A * y + 2 * ly * y2 * y + kap * x * x * x + 2 * a * mu * s * conj(y) - app * y;
  F[2] = 2 * ls * (s2 - a * a * hh * hh / 2) * s + a * mu * y * y - app * s;
}

static double omega2(const cplx *c, const cplx *p, const bgs *b) {
  double a = b->a;
  double x = pow(pow(cabs(c[0]), 4) + creal(p[0] * conj(p[0])) / lx, 0.25);
  double y = pow(pow(cabs(c[1]), 4) + creal(p[1] * conj(p[1])) / ly, 0.25);
  double s = pow(pow(cabs(c[2]), 4) + creal(p[2] * conj(p[2])) / ls, 0.25);
  double cX = CH * rhophi(b) / 3.0;
  double rx = 6 * lx * x * x + 6 * kap * x * y + 3 * kap * x * x + a * a * cX;
  double ry = a * a * MM * MM + 6 * ly * y * y + 2 * a * mu * s + 3 * kap * x * x + 2 * a * mu * y;
  double rs = fabs(2 * ls * (s * s - a * a * hh * hh / 2)) + 4 * ls * s * s + 2 * a * mu * y;
  double m = rx > ry ? rx : ry;
  m = m > rs ? m : rs;
  return m > 1e-30 ? m : 1e-30;
}

#define NQ 18
static const char *QN[NQ] = {"NX", "NY", "NS", "EX", "eps", "NYH", "NYL", "Nsig", "Gauss", "H", "EYexc", "Esig",
                             "absX", "absY", "absS", "rhophi_frac", "Ncirc", "EYfull"};

static void diag(const cplx *c, const cplx *p, const bgs *b, double *q) {
  double a = b->a, ap = aprime(b);
  double n[3]; cplx w[3];
  for (int i = 0; i < 3; i++) { n[i] = -2 * cimag(conj(c[i]) * p[i]); w[i] = p[i] - c[i] * ap / a; }
  double x2 = creal(c[0] * conj(c[0])), y2 = creal(c[1] * conj(c[1]));
  double EX = creal(w[0] * conj(w[0])) + lx * x2 * x2;
  double Nc = 2 * sqrt(2 * lx) * pow(EX / (3 * lx), 0.75);
  double A = a * a * MM * MM;
  cplx B = 2 * a * mu * c[2];
  double bb = cabs(B), D = A * A - bb * bb;
  double NYH = NAN, NYL = NAN, Eexc = NAN;
  if (D > 0 && bb > 0) {
    cplx J = kap * c[0] * c[0] * c[0];
    cplx yp = (-A * J + B * conj(J)) / D;
    cplx Jp = 3 * kap * c[0] * c[0] * p[0];
    cplx Bp = 2 * mu * (ap * c[2] + a * p[2]);
    double Ap = 2 * a * ap * MM * MM;
    cplx J2 = Jp + Ap * yp + Bp * conj(yp);
    cplx ypp = (-A * J2 + B * conj(J2)) / D;
    cplx d = c[1] - yp, dp = p[1] - ypp;
    cplx wd = dp - d * ap / a;                                 /* a^2 d(delta/a)/dt */
    double thp = cimag(conj(c[2]) * p[2]) / creal(c[2] * conj(c[2]));   /* d arg(Sigma)/deta */
    cplx rot = cexp(-I * carg(B) / 2);
    cplx u = rot * d, v = rot * (wd - I * 0.5 * thp * d);
    double u1 = creal(u), u2 = cimag(u), v1 = creal(v), v2 = cimag(v);
    double wH = sqrt(A + bb), wL = sqrt(A - bb);
    double EH = v1 * v1 + wH * wH * u1 * u1, EL = v2 * v2 + wL * wL * u2 * u2;
    NYH = EH / wH; NYL = EL / wL;
    Eexc = creal(wd * conj(wd)) + A * creal(d * conj(d)) + creal(conj(B) * d * d);
  }
  double r = cabs(c[2]) / a;
  double rdot = creal(conj(c[2]) * w[2]) / (a * a * cabs(c[2]));
  double msig = sqrt(2 * ls) * hh;
  double Es = rdot * rdot + ls * (r * r - hh * hh / 2) * (r * r - hh * hh / 2);
  double EYfull = creal(w[1] * conj(w[1])) + A * y2 + ly * y2 * y2;
  double scale = fabs(n[0]) + 3 * fabs(n[1]) + 6 * fabs(n[2]);
  q[0] = n[0]; q[1] = n[1]; q[2] = n[2]; q[3] = EX; q[4] = n[0] / Nc; q[5] = NYH; q[6] = NYL;
  q[7] = a * a * a * Es / msig; q[8] = (n[0] + 3 * n[1] + 6 * n[2]) / (scale > 1e-300 ? scale : 1.0);
  q[9] = Hof(b); q[10] = Eexc; q[11] = a * a * a * a * Es; q[12] = cabs(c[0]) / a; q[13] = cabs(c[1]) / a;
  q[14] = cabs(c[2]) / a; q[15] = (BG == 0) ? 0.0 : rhophi(b) / (3 * Hof(b) * Hof(b)); q[16] = Nc; q[17] = EYfull;
}

static const double YW[8] = {0, -1.61582374150097, -2.44699182370524, -0.716989419708120e-2, 2.44002732616735,
                             0.157739928123617, 1.82020630970714, 1.04242620869991};

int main(int argc, char **argv) {
  int k = 1;
  if (argc < 17) { fprintf(stderr, "usage: see header\n"); return 1; }
  int method = atoi(argv[k++]); double cpar = atof(argv[k++]);
  BG = atoi(argv[k++]); GAM = atof(argv[k++]); CH = atof(argv[k++]);
  lx = atof(argv[k++]); ly = atof(argv[k++]); ls = atof(argv[k++]); kap = atof(argv[k++]);
  MM = atof(argv[k++]); hh = atof(argv[k++]); mu = atof(argv[k++]);
  double aend = atof(argv[k++]); int ns = atoi(argv[k++]); int nw = atoi(argv[k++]);
  if (argc < k + nw + 6) { fprintf(stderr, "usage: see header (windows/state)\n"); return 1; }
  double *we = malloc(sizeof(double) * (nw + 1));
  for (int i = 0; i < nw; i++) we[i] = atof(argv[k++]);
  cplx c[3], p[3];
  for (int i = 0; i < 3; i++) { double re = atof(argv[k++]), im = atof(argv[k++]); c[i] = re + I * im; p[i] = c[i]; }
  if (argc >= k + 6) for (int i = 0; i < 3; i++) { double re = atof(argv[k++]), im = atof(argv[k++]); p[i] += re + I * im; }
  bgs b = {1.0, 1.0, (BG == 0) ? 0.0 : 3.0, (BG == 0) ? 3.0 : 0.0};
  /* chi'(1) = a'(1) z0 + v0 with a'(1) = a^2 H = 1 in both backgrounds */
  double w[8]; w[0] = 1.0;
  for (int i = 1; i < 8; i++) { w[i] = YW[i]; w[0] -= 2 * YW[i]; }
  double sub[15];
  for (int i = 0; i < 7; i++) sub[i] = w[7 - i];
  sub[7] = w[0];
  for (int i = 0; i < 7; i++) sub[8 + i] = w[i + 1];

  double q[NQ], qo[NQ];
  diag(c, p, &b, q);
  printf("S %.17g", b.a); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n");
  double sums[64][NQ], span[64]; memset(sums, 0, sizeof sums); memset(span, 0, sizeof span);
  int si = 1; double nexts = exp(log(aend) * si / ns);
  long steps = 0, rejects = 0; double maxG = 0, hstep = 1e-3;
  static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
    a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
    a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
    b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
    e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
  while (b.a < aend * (1 - 1e-15)) {
    /* step in eta; targets in a are converted with da ~ a' deta (then landed exactly by a short final step) */
    double atarget = nexts < aend ? nexts : aend;
    for (int i = 0; i < nw; i++) if (we[i] > b.a * (1 + 1e-15) && we[i] < atarget) atarget = we[i];
    double ap = aprime(&b);
    double hmax = (atarget - b.a) / ap;   /* first-order estimate; exact landing enforced below */
    double h; int accepted = 1;
    for (int j = 0; j < NQ; j++) qo[j] = q[j];
    double aold = b.a; bgs bold = b;
    if (method == 8) {
      h = cpar / sqrt(omega2(c, p, &b));
      if (h > 0.05) h = 0.05;
      if (h > hmax) h = hmax;
      cplx c0[3], p0[3];
      for (int i = 0; i < 3; i++) { c0[i] = c[i]; p0[i] = p[i]; }
      for (;;) {
        bgs bt = bold;
        for (int st = 0; st < 15; st++) {
          double hs = sub[st] * h;
          for (int i = 0; i < 3; i++) c[i] += 0.5 * hs * p[i];
          bg_advance(&bt, 0.5 * hs);
          cplx F[3]; force(c, &bt, F);
          for (int i = 0; i < 3; i++) p[i] -= hs * F[i];
          for (int i = 0; i < 3; i++) c[i] += 0.5 * hs * p[i];
          bg_advance(&bt, 0.5 * hs);
        }
        double hn = cpar / sqrt(omega2(c, p, &bt));
        if (!isfinite(hn)) hn = 0.25 * h;
        if (h > 1.5 * hn && h > 1e-13) {
          rejects++;
          for (int i = 0; i < 3; i++) { c[i] = c0[i]; p[i] = p0[i]; }
          h = hn; continue;
        }
        b = bt; break;
      }
      steps++;
    } else {
      double rtol = cpar, atol = cpar * 1e-2;
      h = hstep; if (h > hmax) h = hmax;
      double y0[12], k1[12], k2[12], k3[12], k4[12], k5[12], k6[12], k7[12], yt[12], yn[12];
      bgs bs[7];
#define PACK(yy) for (int i = 0; i < 3; i++) { yy[2*i] = creal(c[i]); yy[2*i+1] = cimag(c[i]); yy[6+2*i] = creal(p[i]); yy[7+2*i] = cimag(p[i]); }
#define RHS(bb_, yy, ff) { cplx cc_[3], F_[3]; for (int i = 0; i < 3; i++) cc_[i] = yy[2*i] + I*yy[2*i+1]; force(cc_, bb_, F_); \
        for (int i = 0; i < 3; i++) { ff[2*i] = yy[6+2*i]; ff[2*i+1] = yy[7+2*i]; ff[6+2*i] = -creal(F_[i]); ff[7+2*i] = -cimag(F_[i]); } }
      PACK(y0);
      bs[0] = bold; RHS(&bs[0], y0, k1);
      bs[1] = bold; bg_advance(&bs[1], c2 * h);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * a21 * k1[i]; RHS(&bs[1], yt, k2);
      bs[2] = bold; bg_advance(&bs[2], c3 * h);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a31 * k1[i] + a32 * k2[i]); RHS(&bs[2], yt, k3);
      bs[3] = bold; bg_advance(&bs[3], c4 * h);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]); RHS(&bs[3], yt, k4);
      bs[4] = bold; bg_advance(&bs[4], c5 * h);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]); RHS(&bs[4], yt, k5);
      bs[5] = bold; bg_advance(&bs[5], h);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]); RHS(&bs[5], yt, k6);
      for (int i = 0; i < 12; i++) yn[i] = y0[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
      RHS(&bs[5], yn, k7);
      double err = 0;
      for (int i = 0; i < 12; i++) {
        double e = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
        double sc = atol + rtol * fmax(fabs(y0[i]), fabs(yn[i]));
        err += (e / sc) * (e / sc);
      }
      err = sqrt(err / 12);
      double fac = err > 0 ? 0.9 * pow(err, -0.2) : 5; if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
      if (err <= 1) {
        for (int i = 0; i < 3; i++) { c[i] = yn[2*i] + I * yn[2*i+1]; p[i] = yn[6+2*i] + I * yn[7+2*i]; }
        b = bs[5]; steps++; hstep = (h < hmax) ? h * fac : hstep;
        if (h >= hmax) hstep = fmax(hstep, h);
      } else { rejects++; hstep = h * fac; accepted = 0; }
    }
    if (!accepted) continue;
    if (fabs(b.a - atarget) < 1e-9 * atarget) { b.a = atarget; if (BG == 0) b.eta = atarget; }   /* land on edges exactly */
    diag(c, p, &b, q);
    if (fabs(q[8]) > maxG) maxG = fabs(q[8]);
    if (!isfinite(q[3])) { fprintf(stderr, "nonfinite at a=%g\n", b.a); return 3; }
    double da = b.a - aold;
    for (int i = 0; i + 1 < nw; i++) if (aold >= we[i] * (1 - 1e-12) && b.a <= we[i + 1] * (1 + 1e-12)) {
      span[i] += da; for (int j = 0; j < NQ; j++) sums[i][j] += 0.5 * da * (qo[j] + q[j]);
    }
    if (b.a >= nexts * (1 - 1e-12) || b.a >= aend * (1 - 1e-15)) {
      printf("S %.17g", b.a); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n");
      si++; nexts = si <= ns ? exp(log(aend) * si / ns) : 2 * aend;
    }
    if (steps > 4000000000L) { fprintf(stderr, "step limit\n"); return 4; }
  }
  for (int i = 0; i + 1 < nw; i++) if (span[i] > 0) {
    printf("W %.17g %.17g %.17g", we[i], we[i + 1], span[i]);
    for (int j = 0; j < NQ; j++) printf(" %.17g", sums[i][j] / span[i]);
    printf("\n");
  }
  printf("F");
  for (int i = 0; i < 3; i++) printf(" %.17g %.17g", creal(c[i]), cimag(c[i]));
  for (int i = 0; i < 3; i++) printf(" %.17g %.17g", creal(p[i]), cimag(p[i]));
  printf(" %.17g %.17g %.17g %.17g\n", b.eta, b.a, b.Rphi, b.Rr);
  printf("END %ld %ld %.6g\n", steps, rejects, maxG);
  return 0;
}
