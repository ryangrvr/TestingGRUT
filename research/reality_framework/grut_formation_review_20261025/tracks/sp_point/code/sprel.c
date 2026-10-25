/* sp_point track (Claude, own code, issue #13 review of Codex 373d37e).
 *
 * Homogeneous X, Y, Sigma release after inflation into radiation, H_i = 1 units.
 * Conformal fields chi = a z / H_i, conformal time with a = H_i eta, a(1) = 1, a' = 1 (so H = H_i / a^2).
 * Derivation (own): z'' + 3 H z' + dV/dz* = 0  ==>  chi'' - (a''/a) chi + a^3 F(chi/a) = 0, a'' = 0 in radiation.
 *   U(chi, a) = lx|x|^4 + ly|y|^4 + ls(|s|^2 - a^2 h^2/2)^2 + a^2 M^2 |y|^2
 *               + 2 Re(k y* x^3) + 2 a mu Re(s* y^2),          chi'' = -dU/dchi*.
 * Zero physical velocity at release: chi'(1) = chi(1).
 *
 * Integrators (both own implementations):
 *   method 8: Yoshida (1990) 8th-order symmetric composition of the drift-kick-drift leapfrog in extended phase space
 *             (time advances with the drift). Step h = c / Omega(state, a), Omega a Gershgorin bound on the local
 *             linear frequency. Gauss law is preserved to rounding by construction (drift and gauge-covariant kick).
 *             The tangent (variational) map is the exact derivative of the discrete map.
 *   method 5: Dormand-Prince 5(4) adaptive (no tangent), used as an independent cross-check.
 *
 * Diagnostics: N_i = -2 Im(chi_i* chi_i'), E_X = |chi_X' - chi_X/a|^2 + lx|chi_X|^4,
 *   N_circ = 2 sqrt(2 lx) (E_X/(3 lx))^(3/4), eps = N_X/N_circ,
 *   J_Y (Codex definition: forced-solution subtraction, E_exc/(a M)) and J_Y,eig (excitation split into the two
 *   Majorana-split real eigenmodes, each divided by its own frequency sqrt(A +- |B|)).
 *   With tangent: d(eps)/d(param) and d(N_X)/d(param) for the supplied initial tangent vector.
 *
 * Usage:
 *   sprel method c_or_rtol M h mu a_end nsamp tangent(0/1) nwin w0 w1 ... w_nwin  x0r x0i y0r y0i s0r s0i
 *         [vx0r vx0i vy0r vy0i vs0r vs0i]   (initial physical velocities dz/dt in H_i^2 units; default 0)
 *         [tangent: 12 numbers dchi(6) dpi(6)]
 * Output lines: S (samples), W (window averages), F (final state), END (counters).
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double lx = 1e-3, ly = 1.0, ls = 1.0, kap = 1e-3;
static double MM, hh, mu;

typedef double complex cplx;

static void force(const cplx *c, double a, cplx *F) {
  cplx x = c[0], y = c[1], s = c[2];
  double x2 = creal(x * conj(x)), y2 = creal(y * conj(y)), s2 = creal(s * conj(s));
  double A = a * a * MM * MM;
  cplx B = 2 * a * mu * s;
  F[0] = 2 * lx * x2 * x + 3 * kap * y * conj(x) * conj(x);
  F[1] = A * y + 2 * ly * y2 * y + kap * x * x * x + B * conj(y);
  F[2] = 2 * ls * (s2 - a * a * hh * hh / 2) * s + a * mu * y * y;
}

static void dforce(const cplx *c, double a, const cplx *d, cplx *dF) {
  cplx x = c[0], y = c[1], s = c[2], dx = d[0], dy = d[1], ds = d[2];
  double x2 = creal(x * conj(x)), y2 = creal(y * conj(y)), s2 = creal(s * conj(s));
  double A = a * a * MM * MM;
  dF[0] = 2 * lx * (2 * x2 * dx + x * x * conj(dx)) + 3 * kap * (dy * conj(x) * conj(x) + 2 * y * conj(x) * conj(dx));
  dF[1] = A * dy + 2 * ly * (2 * y2 * dy + y * y * conj(dy)) + 3 * kap * x * x * dx
          + 2 * a * mu * (ds * conj(y) + s * conj(dy));
  dF[2] = 2 * ls * ((s2 - a * a * hh * hh / 2) * ds + (s * conj(ds) + conj(s) * ds) * s) + 2 * a * mu * y * dy;
}

static double omega2(const cplx *c, const cplx *p, double a) {
  /* Gershgorin bound on the linear frequency, evaluated at the quartic turning-point amplitude of each species
     (|chi|^4 + |chi'|^2/lambda)^(1/4), so that fast passages through small |chi| do not shrink the bound. */
  double x = pow(pow(cabs(c[0]), 4) + creal(p[0] * conj(p[0])) / lx, 0.25);
  double y = pow(pow(cabs(c[1]), 4) + creal(p[1] * conj(p[1])) / ly, 0.25);
  double s = pow(pow(cabs(c[2]), 4) + creal(p[2] * conj(p[2])) / ls, 0.25);
  double rx = 6 * lx * x * x + 6 * kap * x * y + 3 * kap * x * x;
  double ry = a * a * MM * MM + 6 * ly * y * y + 2 * a * mu * s + 3 * kap * x * x + 2 * a * mu * y;
  double rs = fabs(2 * ls * (s * s - a * a * hh * hh / 2)) + 4 * ls * s * s + 2 * a * mu * y;
  double m = rx > ry ? rx : ry;
  return m > rs ? m : rs;
}

/* diagnostics */
#define NQ 16
static const char *QN[NQ] = {"NX", "NY", "NS", "Neff", "EX", "epsX", "NcircX", "JY", "JYeig", "EYexc",
                             "EYtot", "Gauss", "deps", "dNX", "lnT", "absY"};
static int TANG = 0;
static double LNSCALE = 0; static int RESCALED = 0;

static void diag(double a, const cplx *c, const cplx *p, const cplx *dc, const cplx *dp, double *q) {
  double n[3];
  for (int i = 0; i < 3; i++) n[i] = -2 * cimag(conj(c[i]) * p[i]);
  cplx wx = p[0] - c[0] / a;
  double x2 = creal(c[0] * conj(c[0])), y2 = creal(c[1] * conj(c[1]));
  double EX = creal(wx * conj(wx)) + lx * x2 * x2;
  double Nc = 2 * sqrt(2 * lx) * pow(EX / (3 * lx), 0.75);
  double A = a * a * MM * MM;
  cplx B = 2 * a * mu * c[2];
  double D = A * A - creal(B * conj(B));
  double JY = -1, JYe = -1, Eexc = -1;
  if (D > 0) {
    cplx J = kap * c[0] * c[0] * c[0];
    cplx yp = (-A * J + B * conj(J)) / D;
    cplx JD = 3 * kap * c[0] * c[0] * p[0] + 2 * MM * MM * a * yp + 2 * mu * (c[2] + a * p[2]) * conj(yp);
    cplx ypd = (-A * JD + B * conj(JD)) / D;
    cplx dy = c[1] - yp, dw = p[1] - ypd - dy / a;
    Eexc = creal(dw * conj(dw)) + A * creal(dy * conj(dy)) + creal(conj(B) * dy * dy);
    JY = Eexc / (a * MM);
    /* eigenmode split: d = e^{i th/2}(u + i w), Re(B* d^2) = |B|(u^2 - w^2) */
    double th = carg(B), bb = cabs(B);
    cplx rot = cexp(-I * th / 2);
    cplx dd = rot * dy, vv = rot * dw;
    double u = creal(dd), w = cimag(dd), uv = creal(vv), wv = cimag(vv);
    double Eu = uv * uv + (A + bb) * u * u, Ew = wv * wv + (A - bb) * w * w;
    JYe = Eu / sqrt(A + bb) + Ew / sqrt(A - bb);
  }
  cplx wy = p[1] - c[1] / a;
  double EYt = creal(wy * conj(wy)) + A * y2 + ly * y2 * y2;
  q[0] = n[0]; q[1] = n[1]; q[2] = n[2]; q[3] = n[0] + 3 * n[1];
  q[4] = EX; q[5] = n[0] / Nc; q[6] = Nc; q[7] = JY; q[8] = JYe; q[9] = Eexc; q[10] = EYt;
  q[11] = n[0] + 3 * n[1] + 6 * n[2];
  q[12] = q[13] = q[14] = 0;
  if (TANG) {
    double dN = -2 * cimag(conj(dc[0]) * p[0] + conj(c[0]) * dp[0]);
    cplx dwx = dp[0] - dc[0] / a;
    double dE = 2 * creal(conj(wx) * dwx) + 4 * lx * x2 * creal(conj(c[0]) * dc[0]);
    double dNc = 0.75 * Nc * dE / EX;
    q[12] = dN / Nc - n[0] * dNc / (Nc * Nc);
    q[13] = dN;
    double t = 0;
    for (int i = 0; i < 3; i++) t += creal(dc[i] * conj(dc[i])) + creal(dp[i] * conj(dp[i]));
    q[14] = 0.5 * log(t) + LNSCALE;
    if (RESCALED) { q[12] = NAN; q[13] = NAN; }
  }
  q[15] = sqrt(y2);
}

/* Yoshida solution A */
static const double YW[8] = {0, -1.61582374150097, -2.44699182370524, -0.716989419708120e-2, 2.44002732616735,
                             0.157739928123617, 1.82020630970714, 1.04242620869991};

int main(int argc, char **argv) {
  int k = 1;
  if (argc < 12) { fprintf(stderr, "usage, see header\n"); return 1; }
  int method = atoi(argv[k++]);
  double cpar = atof(argv[k++]);
  MM = atof(argv[k++]); hh = atof(argv[k++]); mu = atof(argv[k++]);
  double aend = atof(argv[k++]);
  int ns = atoi(argv[k++]);
  TANG = atoi(argv[k++]);
  int nw = atoi(argv[k++]);
  double *we = malloc(sizeof(double) * (nw + 1));
  for (int i = 0; i < nw; i++) we[i] = atof(argv[k++]);
  cplx c[3], p[3], dc[3] = {0, 0, 0}, dp[3] = {0, 0, 0};
  for (int i = 0; i < 3; i++) { double re = atof(argv[k++]), im = atof(argv[k++]); c[i] = re + I * im; p[i] = c[i]; }
  if (argc >= k + 6) { /* initial physical velocities dz/dt; chi' = z + a^2 dz/dt at a = 1 */
    for (int i = 0; i < 3; i++) { double re = atof(argv[k++]), im = atof(argv[k++]); p[i] += re + I * im; }
  }
  if (TANG) {
    if (argc < k + 12) { fprintf(stderr, "tangent needs 12 numbers\n"); return 1; }
    for (int i = 0; i < 3; i++) { double re = atof(argv[k++]), im = atof(argv[k++]); dc[i] = re + I * im; }
    for (int i = 0; i < 3; i++) { double re = atof(argv[k++]), im = atof(argv[k++]); dp[i] = re + I * im; }
  }
  double w[8];
  w[0] = 1.0;
  for (int i = 1; i < 8; i++) { w[i] = YW[i]; w[0] -= 2 * YW[i]; }
  double sub[15];
  for (int i = 0; i < 7; i++) sub[i] = w[7 - i];
  sub[7] = w[0];
  for (int i = 0; i < 7; i++) sub[8 + i] = w[i + 1];

  double a = 1.0, q[NQ], qo[NQ];
  diag(a, c, p, dc, dp, q);
  double G0 = fabs(q[6]);
  printf("S %.17g", a); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n");
  /* windows */
  double sums[64][NQ], sq[64][3], span[64];
  memset(sums, 0, sizeof sums); memset(sq, 0, sizeof sq); memset(span, 0, sizeof span);
  int si = 1;
  double nexts = exp(log(aend) * si / ns);
  long steps = 0, rejects = 0;
  double maxG = 0;
  double hstep = 1e-3;
  /* DP5(4) tableau */
  static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
                      a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
                      a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
                      b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
                      e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525,
                      e7 = -1. / 40;
  while (a < aend * (1 - 1e-15)) {
    double target = nexts < aend ? nexts : aend;
    for (int i = 0; i < nw; i++) if (we[i] > a * (1 + 1e-15) && we[i] < target) target = we[i];
    double h;
    int accepted = 1;
    for (int j = 0; j < NQ; j++) qo[j] = q[j];
    double aold = a;
    if (method == 8) {
      h = cpar / sqrt(omega2(c, p, a));
      if (h > 0.05) h = 0.05;
      if (a + h > target) h = target - a;
      cplx c0[3], p0[3], dc0[3], dp0[3];
      for (int i = 0; i < 3; i++) { c0[i] = c[i]; p0[i] = p[i]; dc0[i] = dc[i]; dp0[i] = dp[i]; }
      for (;;) {
        /* composition: for each substep: drift sub/2, kick sub, drift sub/2 */
        double t = a;
        for (int st = 0; st < 15; st++) {
          double hs = sub[st] * h;
          for (int i = 0; i < 3; i++) { c[i] += 0.5 * hs * p[i]; if (TANG) dc[i] += 0.5 * hs * dp[i]; }
          t += 0.5 * hs;
          cplx F[3], dF[3];
          force(c, t, F);
          if (TANG) dforce(c, t, dc, dF);
          for (int i = 0; i < 3; i++) { p[i] -= hs * F[i]; if (TANG) dp[i] -= hs * dF[i]; }
          for (int i = 0; i < 3; i++) { c[i] += 0.5 * hs * p[i]; if (TANG) dc[i] += 0.5 * hs * dp[i]; }
          t += 0.5 * hs;
        }
        /* guard: if the local frequency at the end of the step demands a step < h/1.5, redo with that step */
        double hn = cpar / sqrt(omega2(c, p, aold + h));
        if (!isfinite(hn)) hn = 0.25 * h;
        if (h > 1.5 * hn && h > 1e-12) {
          rejects++;
          for (int i = 0; i < 3; i++) { c[i] = c0[i]; p[i] = p0[i]; dc[i] = dc0[i]; dp[i] = dp0[i]; }
          h = hn;
          continue;
        }
        break;
      }
      a = (h == target - aold) ? target : aold + h;
      steps++;
      if (TANG) { /* renormalise the tangent to avoid overflow in chaotic runs; keep the log scale */
        double tt = 0;
        for (int i = 0; i < 3; i++) tt += creal(dc[i] * conj(dc[i])) + creal(dp[i] * conj(dp[i]));
        if (tt > 1e200) {
          for (int i = 0; i < 3; i++) { dc[i] *= 1e-100; dp[i] *= 1e-100; }
          LNSCALE += 100 * log(10.0); RESCALED = 1;
        }
      }
    } else {
      /* DP5(4) on y = (chi, pi), 12 real components */
      double rtol = cpar, atol = cpar * 1e-2;
      double htrial = hstep;
      h = htrial;
      if (a + h > target) h = target - a;
      double y0[12], k1[12], k2[12], k3[12], k4[12], k5[12], k6[12], k7[12], yt[12], yn[12];
#define PACK(yy, cc, pp) for (int i = 0; i < 3; i++) { yy[2*i] = creal(cc[i]); yy[2*i+1] = cimag(cc[i]); yy[6+2*i] = creal(pp[i]); yy[7+2*i] = cimag(pp[i]); }
#define RHS(tt, yy, ff) { cplx cc_[3], F_[3]; for (int i = 0; i < 3; i++) cc_[i] = yy[2*i] + I*yy[2*i+1]; force(cc_, tt, F_); \
        for (int i = 0; i < 3; i++) { ff[2*i] = yy[6+2*i]; ff[2*i+1] = yy[7+2*i]; ff[6+2*i] = -creal(F_[i]); ff[7+2*i] = -cimag(F_[i]); } }
      PACK(y0, c, p);
      RHS(a, y0, k1);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * a21 * k1[i];
      RHS(a + c2 * h, yt, k2);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a31 * k1[i] + a32 * k2[i]);
      RHS(a + c3 * h, yt, k3);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      RHS(a + c4 * h, yt, k4);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      RHS(a + c5 * h, yt, k5);
      for (int i = 0; i < 12; i++) yt[i] = y0[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      RHS(a + h, yt, k6);
      for (int i = 0; i < 12; i++) yn[i] = y0[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
      RHS(a + h, yn, k7);
      double err = 0;
      for (int i = 0; i < 12; i++) {
        double e = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
        double sc = atol + rtol * fmax(fabs(y0[i]), fabs(yn[i]));
        err += (e / sc) * (e / sc);
      }
      err = sqrt(err / 12);
      double fac = err > 0 ? 0.9 * pow(err, -0.2) : 5;
      if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
      if (err <= 1) {
        for (int i = 0; i < 3; i++) { c[i] = yn[2*i] + I * yn[2*i+1]; p[i] = yn[6+2*i] + I * yn[7+2*i]; }
        a = (h == target - aold) ? target : aold + h;
        steps++;
        hstep = (h < htrial) ? htrial : h * fac;
      } else { rejects++; hstep = h * fac; accepted = 0; }
    }
    if (!accepted) continue;
    diag(a, c, p, dc, dp, q);
    double g = fabs(q[11]) / G0; if (g > maxG) maxG = g;
    if (!isfinite(q[4])) { fprintf(stderr, "nonfinite at a=%g\n", a); return 3; }
    double hh_ = a - aold;
    for (int i = 0; i + 1 < nw; i++) {
      if (aold >= we[i] * (1 - 1e-15) && a <= we[i + 1] * (1 + 1e-15)) {
        span[i] += hh_;
        for (int j = 0; j < NQ; j++) sums[i][j] += 0.5 * hh_ * (qo[j] + q[j]);
        int jj[3] = {0, 3, 5};
        for (int m = 0; m < 3; m++) sq[i][m] += 0.5 * hh_ * (qo[jj[m]] * qo[jj[m]] + q[jj[m]] * q[jj[m]]);
      }
    }
    if (a >= nexts * (1 - 1e-14) || a >= aend * (1 - 1e-15)) {
      printf("S %.17g", a); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n");
      si++; nexts = si <= ns ? exp(log(aend) * si / ns) : 2 * aend;
    }
    if (steps > 4000000000L) { fprintf(stderr, "step limit\n"); return 4; }
  }
  for (int i = 0; i + 1 < nw; i++) if (span[i] > 0) {
    printf("W %.17g %.17g %.17g", we[i], we[i + 1], span[i]);
    for (int j = 0; j < NQ; j++) printf(" %.17g", sums[i][j] / span[i]);
    for (int m = 0; m < 3; m++) printf(" %.17g", sq[i][m] / span[i]);
    printf("\n");
  }
  printf("F");
  for (int i = 0; i < 3; i++) printf(" %.17g %.17g", creal(c[i]), cimag(c[i]));
  for (int i = 0; i < 3; i++) printf(" %.17g %.17g", creal(p[i]), cimag(p[i]));
  printf("\n");
  printf("END %ld %ld %.6g %d\n", steps, rejects, maxG, RESCALED);
  return 0;
}
