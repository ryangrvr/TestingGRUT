/* hw_gates track (Claude, own code): homogeneous X, Y, Sigma release at run-time couplings, cosmic time, physical
 * fields, Dormand-Prince 5(4); instant radiation release (a = sqrt(2t), H = 1/(2t)), H_i = 1 units.
 * Lineage: our ls_scaling scalrel.c / sp_yledger yrelease.c (same session, our code). Additions here:
 *   - unwrapped phases theta_X, theta_Sigma and the gauge-invariant Goldstone angle theta_G (charges 1, 3, 6);
 *   - min |Sigma|/vev over the run; Higgs norm F^2 = 2 sum q_i^2 |z_i|^2 at samples;
 *   - supplied, tagged test terms that may persist into the release (environment HWG_*):
 *       V_P = gA(a) [ 2Re(kapP Y* X^3) + 2Re(muP Sigma* Y^2) + 2Re(c7 X^6 Sigma*) ] - gm(a) cX |X|^2
 *     pmode 0: gA = gm = 0 (instant removal); 1: gA = gm = exp(-(a-1)/Dl); 2: gA = a^-2, gm = a^-4 (terms that
 *     scale with H and H^2 in radiation; DRT-like persistence); 3: gA = gm = 1 (never removed).
 * V = lx|X|^4 + ly|Y|^4 + ls(|S|^2-h^2/2)^2 + M^2|Y|^2 + 2Re(kap Y* X^3) + 2Re(mu S* Y^2) (+ V_P).
 * Usage: hwgrel lx ly ls kap M h mu t_end rtol atol nsamp a_win0 x0r x0i y0r y0i s0r s0i
 * Output: S lines (samples), one W line (da-weighted window mean over a in [a_win0, a_end]), END counters.
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

static double lx, ly, ls, kap, MM, hh, mu;
static int pmode = 0; static double Dl = 1.0, cX = 0.0; static double complex muP = 0, kapP = 0, c7 = 0;

static void gfac(double a, double *gA, double *gm) {
  if (pmode == 1) { *gA = *gm = exp(-(a - 1) / Dl); }
  else if (pmode == 2) { *gA = 1 / (a * a); *gm = 1 / (a * a * a * a); }
  else if (pmode == 3) { *gA = *gm = 1.0; }
  else { *gA = *gm = 0.0; }
}

/* state: 0..5 z (re,im)x3, 6..11 z' (re,im)x3, 12 thX, 13 thS, 14 W (energy ledger), 15 ln a (unused) */
#define N 16

static void force(const double complex *z, double a, double complex *F) {
  double complex X = z[0], Y = z[1], S = z[2];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  F[0] = 2 * lx * ax2 * X + 3 * kap * Y * conj(X) * conj(X);
  F[1] = MM * MM * Y + 2 * ly * ay2 * Y + kap * X * X * X + 2 * mu * S * conj(Y);
  F[2] = 2 * ls * (as2 - hh * hh / 2) * S + mu * Y * Y;
  if (pmode) {
    double gA, gm; gfac(a, &gA, &gm);
    double complex Xc = conj(X), Xc2 = Xc * Xc;
    F[0] += gA * (3 * conj(kapP) * Y * Xc2 + 6 * conj(c7) * Xc2 * Xc2 * Xc * S) - gm * cX * X;
    F[1] += gA * (kapP * X * X * X + 2 * conj(muP) * S * conj(Y));
    F[2] += gA * (muP * Y * Y + c7 * X * X * X * X * X * X);
  }
}

static double potential(const double complex *z, double a) {
  double complex X = z[0], Y = z[1], S = z[2];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  double v = lx * ax2 * ax2 + ly * ay2 * ay2 + ls * (as2 - hh * hh / 2) * (as2 - hh * hh / 2) + MM * MM * ay2
       + 2 * creal(kap * conj(Y) * X * X * X) + 2 * creal(mu * conj(S) * Y * Y);
  if (pmode) {
    double gA, gm; gfac(a, &gA, &gm);
    v += gA * 2 * creal(kapP * conj(Y) * X * X * X + muP * conj(S) * Y * Y + c7 * X * X * X * X * X * X * conj(S)) - gm * cX * ax2;
  }
  return v;
}

static void rhs(double t, const double *y, double *dy) {
  double a = sqrt(2 * t), H = 1.0 / (2 * t);
  double complex z[3], zd[3], F[3];
  for (int i = 0; i < 3; i++) { z[i] = y[2 * i] + I * y[2 * i + 1]; zd[i] = y[6 + 2 * i] + I * y[7 + 2 * i]; }
  force(z, a, F);
  double K = 0;
  for (int i = 0; i < 3; i++) {
    dy[2 * i] = creal(zd[i]); dy[2 * i + 1] = cimag(zd[i]);
    double complex acc = -3 * H * zd[i] - F[i];
    dy[6 + 2 * i] = creal(acc); dy[7 + 2 * i] = cimag(acc);
    K += creal(zd[i] * conj(zd[i]));
  }
  double rX2 = creal(z[0] * conj(z[0])), rS2 = creal(z[2] * conj(z[2]));
  dy[12] = cimag(conj(z[0]) * zd[0]) / rX2;
  dy[13] = cimag(conj(z[2]) * zd[2]) / rS2;
  double Vv = potential(z, a);
  /* energy ledger d(a^3 rho)/dt = 3 H a^3 (V - K) + a^3 dV/dt|explicit; explicit term only for pmode != 0 (numerical) */
  double dVdt = 0;
  if (pmode) { double da = 1e-6 * a; dVdt = (potential(z, a + da) - potential(z, a - da)) / (2 * da) * H * a; }
  dy[14] = 3 * H * a * a * a * (Vv - K) + a * a * a * dVdt;
  dy[15] = H;
}

#define NO 26
static void diagnostics(double t, const double *y, double *o) {
  double a = sqrt(2 * t), H = 1.0 / (2 * t);
  double complex z[3], zd[3];
  double n[3], K = 0;
  for (int i = 0; i < 3; i++) {
    z[i] = y[2 * i] + I * y[2 * i + 1]; zd[i] = y[6 + 2 * i] + I * y[7 + 2 * i];
    n[i] = -2 * a * a * a * cimag(conj(z[i]) * zd[i]);
    K += creal(zd[i] * conj(zd[i]));
  }
  double a3 = a * a * a, a4 = a3 * a;
  double ax2 = creal(z[0] * conj(z[0])), ay2 = creal(z[1] * conj(z[1])), as2 = creal(z[2] * conj(z[2]));
  double EX = a4 * (creal(zd[0] * conj(zd[0])) + lx * ax2 * ax2);
  double Nc = 2 * sqrt(2 * lx) * pow(EX / (3 * lx), 0.75);
  /* forced Y solution (static) and the excitation split into the two Majorana-split real eigenmodes */
  double complex X = z[0], Xd = zd[0], S = z[2], Sd = zd[2];
  double M2 = MM * MM, as = cabs(S);
  double D = M2 * M2 - 4 * mu * mu * as * as;
  double complex J = kap * X * X * X, Jd = 3 * kap * X * X * Xd;
  double complex yp = (-M2 * J + 2 * mu * S * conj(J)) / D;
  double complex J2 = Jd + 2 * mu * Sd * conj(yp);
  double complex ypd = (-M2 * J2 + 2 * mu * S * conj(J2)) / D;
  double complex d = z[1] - yp, dd = zd[1] - ypd;
  double phid = cimag(conj(S) * Sd) / (as * as);
  double complex rot = cexp(-I * carg(S) / 2);
  double complex u = d * rot, ud = (dd - I * 0.5 * phid * d) * rot;
  double u1 = sqrt(2.) * creal(u), u2 = sqrt(2.) * cimag(u), u1d = sqrt(2.) * creal(ud), u2d = sqrt(2.) * cimag(ud);
  double w1 = sqrt(M2 + 2 * mu * as), w2 = sqrt(M2 - 2 * mu * as);
  double E1 = 0.5 * (u1d * u1d + w1 * w1 * u1 * u1), E2 = 0.5 * (u2d * u2d + w2 * w2 * u2 * u2);
  double Ecod = creal(dd * conj(dd)) + M2 * creal(d * conj(d)) + 2 * mu * creal(conj(S) * d * d);
  double r = as, rd = creal(conj(S) * Sd) / as;
  double Es = rd * rd + ls * (r * r - hh * hh / 2) * (r * r - hh * hh / 2);
  double msig = sqrt(2 * ls) * hh;
  double rho = K + potential(z, a);
  double scale = fabs(n[0]) + 3 * fabs(n[1]) + 6 * fabs(n[2]);
  /* Goldstone angle: theta_G = thS/6 + wX rX + wY rY, r_X = thX - thS/6, r_Y = principal(arg Y - thS/2) */
  double F2 = 2 * (ax2 + 9 * ay2 + 36 * as2);
  double wX = 2 * ax2 / F2, wY = 6 * ay2 / F2;   /* q_i 2|z_i|^2 / F^2 with q = 1, 3 */
  double thX = y[12], thS = y[13];
  double rY = carg(z[1]) - thS / 2; rY -= 2 * M_PI * floor((rY + M_PI) / (2 * M_PI));
  double thG = thS / 6 + wX * (thX - thS / 6) + wY * rY;
  o[0] = t; o[1] = a; o[2] = H; o[3] = n[0]; o[4] = n[1]; o[5] = n[2]; o[6] = EX; o[7] = Nc; o[8] = n[0] / Nc;
  o[9] = a3 * E1 / w1; o[10] = a3 * E2 / w2; o[11] = a3 * Ecod / MM; o[12] = a3 * Es / msig;
  o[13] = (n[0] + 3 * n[1] + 6 * n[2]) / (scale > 1e-300 ? scale : 1.0);
  o[14] = a3 * rho - y[14];
  o[15] = a4 * rho; o[16] = 2 * mu * as / M2; o[17] = a4 * (E1 + E2); o[18] = a4 * Es;
  o[19] = thX; o[20] = thS; o[21] = thG; o[22] = as / (hh / sqrt(2.0)); o[23] = sqrt(ax2); o[24] = sqrt(ay2); o[25] = F2;
}

static double genv(const char *k, double d) { const char *s = getenv(k); return s ? atof(s) : d; }

int main(int argc, char **argv) {
  if (argc != 19) { fprintf(stderr, "usage: hwgrel lx ly ls kap M h mu t_end rtol atol nsamp a_win0 x0r x0i y0r y0i s0r s0i\n"); return 1; }
  lx = atof(argv[1]); ly = atof(argv[2]); ls = atof(argv[3]); kap = atof(argv[4]);
  MM = atof(argv[5]); hh = atof(argv[6]); mu = atof(argv[7]);
  double t_end = atof(argv[8]), rtol = atof(argv[9]), atol = atof(argv[10]);
  int nsamp = atoi(argv[11]); double awin0 = atof(argv[12]);
  pmode = (int)genv("HWG_PMODE", 0); Dl = genv("HWG_DL", 1.0); cX = genv("HWG_CX", 0.0);
  muP = genv("HWG_MUPR", 0) + I * genv("HWG_MUPI", 0); kapP = genv("HWG_KPR", 0) + I * genv("HWG_KPI", 0);
  c7 = genv("HWG_C7R", 0) + I * genv("HWG_C7I", 0);
  double y[N] = {0};
  for (int i = 0; i < 6; i++) { y[i] = atof(argv[13 + i]); y[6 + i] = 0.0; }
  double t = 0.5;
  { double complex z[3]; for (int i = 0; i < 3; i++) z[i] = y[2 * i] + I * y[2 * i + 1];
    y[12] = carg(z[0]); y[13] = carg(z[2]); y[14] = potential(z, 1.0); }
  double W0 = fabs(y[14]) > 0 ? fabs(y[14]) : 1.0;
  double wsum[NO] = {0}, wspan = 0;
  static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
    a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
    a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
    b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
    e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
  double k1[N], k2[N], k3[N], k4[N], k5[N], k6[N], k7[N], yt[N], yn[N];
  double h = 1e-4; long steps = 0, rejects = 0;
  double t0 = t, lr = log(t_end / t0); int si = 0;
  double next_s = t0 * exp(lr * (si + 1) / nsamp);
  double o[NO], oo[NO]; double maxg = 0, maxl = 0, minS = 1e300;
  diagnostics(t, y, o);
  printf("S"); for (int j = 0; j < NO; j++) printf(" %.15g", o[j]); printf("\n");
  rhs(t, y, k1);
  while (t < t_end) {
    int land = 0;
    if (t + h >= next_s) { h = next_s - t; land = 1; }
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * a21 * k1[i];
    rhs(t + c2 * h, yt, k2);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
    rhs(t + c3 * h, yt, k3);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    rhs(t + c4 * h, yt, k4);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    rhs(t + c5 * h, yt, k5);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    rhs(t + h, yt, k6);
    for (int i = 0; i < N; i++) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    rhs(t + h, yn, k7);
    double err = 0;
    for (int i = 0; i < 14; i++) {
      double ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
      double sc = atol + rtol * fmax(fabs(y[i]), fabs(yn[i]));
      double rr = fabs(ei) / sc; if (rr > err) err = rr;
    }
    if (err <= 1.0) {
      for (int j = 0; j < NO; j++) oo[j] = o[j];
      t += h; steps++;
      for (int i = 0; i < N; i++) { y[i] = yn[i]; k1[i] = k7[i]; }
      diagnostics(t, y, o);
      if (fabs(o[13]) > maxg) maxg = fabs(o[13]);
      if (fabs(o[14]) / W0 > maxl) maxl = fabs(o[14]) / W0;
      if (o[22] < minS) minS = o[22];
      if (o[1] > awin0) { double da = o[1] - oo[1]; wspan += da; for (int j = 0; j < NO; j++) wsum[j] += 0.5 * da * (o[j] + oo[j]); }
      if (land) {
        printf("S"); for (int j = 0; j < NO; j++) printf(" %.15g", o[j]); printf("\n");
        si++; next_s = (si + 1 <= nsamp) ? t0 * exp(lr * (si + 1) / nsamp) : t_end * 2;
        if (si >= nsamp) break;
      }
      double fac = (err > 0) ? 0.9 * pow(err, -0.2) : 5.0;
      if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
      h *= fac;
    } else {
      rejects++;
      double fac = 0.9 * pow(err, -0.2); if (fac < 0.1) fac = 0.1;
      h *= fac;
    }
    if (!isfinite(y[0]) || h < 1e-15) { printf("FAIL %.17g\n", t); return 2; }
  }
  if (wspan > 0) { printf("W %.15g", wspan); for (int j = 0; j < NO; j++) printf(" %.15g", wsum[j] / wspan); printf("\n"); }
  printf("END %ld %ld %.6g %.6g %.6g\n", steps, rejects, maxg, maxl, minS);
  return 0;
}
