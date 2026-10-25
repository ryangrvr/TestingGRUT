/* Track sp_vector (Claude, own code): homogeneous X, Y, Sigma release with the eaten-phase bookkeeping.
 *
 * Conformal fields chi_i = a z_i / H_i, a = H_i eta, a(1) = 1, radiation era (Codex MATHEMATICS.md dictionary,
 * equations written independently):
 *   U = lx|cX|^4 + ly|cY|^4 + ls(|cS|^2 - a^2 h^2/2)^2 + a^2 M^2 |cY|^2 + 2 Re[kap cY* cX^3 + a mu cS* cY^2]
 *   chi_i'' = -dU/dchi_i*,  prime = d/da.   Zero physical velocity at a = 1: chi_i' = chi_i.
 * Extra ODE variables:
 *   theta_X, theta_S  unwrapped phases, theta' = Im(chi* chi')/|chi|^2
 *   window integrals over [aw0, a_end] of N_X, N_Y, N_S, E_X, N_circ(E_X), eps_X, theta_G
 * Gauge-invariant late Goldstone angle (charges q = 1, 3, 6):
 *   theta_G = theta_S/6 + w_X r_X + w_Y r_Y,  r_X = theta_X - theta_S/6 (unwrapped),
 *   r_Y = principal(arg cY - theta_S/2),  w_i = q_i 2|chi_i|^2 / sum_j q_j^2 2|chi_j|^2.
 *   theta_S(inf) is also estimated as theta_S(a_end) + tail, tail = -<N_S>/(h^2 a_end) (vacuum |Sigma|).
 * Integrator: Dormand-Prince 5(4), adaptive, own implementation.
 * Usage: vrel delta sfac a_end aw0 rtol HI [histfile|-] [xfac yfac]
 *   initial: X = +X0, Y = -Y0, Sigma = -S0*sfac*exp(i delta) (aligned at delta = 0), minimum from vmodel.py
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

static double lx = 1e-3, ly = 1.0, ls = 1.0, kap = 1e-3, mu, M, hh;
static double aw0;

#define NV 23
/* y: 0-5 chi (re,im) X,Y,S; 6-11 chi'; 12 thX; 13 thS; 14..20 window integrals NX NY NS EX Nc eps thG; 21,22 spare */

static void dU(double a, const double complex *c, double complex *F) {
  double complex X = c[0], Y = c[1], S = c[2];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  double complex Xc = conj(X), Yc = conj(Y);
  F[0] = 2 * lx * ax2 * X + 3 * kap * Y * Xc * Xc;
  F[1] = 2 * ly * ay2 * Y + a * a * M * M * Y + kap * X * X * X + 2 * a * mu * S * Yc;
  F[2] = 2 * ls * (as2 - a * a * hh * hh / 2) * S + a * mu * Y * Y;
}

static double princ(double x) { return x - 2 * M_PI * floor((x + M_PI) / (2 * M_PI)); }

static double thetaG(const double *y) {
  double complex X = y[0] + I * y[1], Y = y[2] + I * y[3], S = y[4] + I * y[5];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  double F2 = 2 * (ax2 + 9 * ay2 + 36 * as2);
  double wX = 2 * ax2 / F2, wY = 6 * ay2 / F2;
  double thS = y[13], thX = y[12];
  double rX = thX - thS / 6;
  double rY = (ay2 > 0) ? princ(carg(Y) - thS / 2) : 0.0;
  return thS / 6 + wX * rX + wY * rY;
}

static void charges(const double *y, double *N) {
  for (int i = 0; i < 3; i++) {
    double complex c = y[2 * i] + I * y[2 * i + 1], p = y[6 + 2 * i] + I * y[7 + 2 * i];
    N[i] = -2 * cimag(conj(c) * p);
  }
}

static double EX(double a, const double *y) {
  double complex c = y[0] + I * y[1], p = y[6] + I * y[7];
  double complex d = p - c / a;
  double r2 = creal(c * conj(c));
  return creal(d * conj(d)) + lx * r2 * r2;
}

static void rhs(double a, const double *y, double *dy) {
  double complex c[3], p[3], F[3];
  for (int i = 0; i < 3; i++) { c[i] = y[2 * i] + I * y[2 * i + 1]; p[i] = y[6 + 2 * i] + I * y[7 + 2 * i]; }
  dU(a, c, F);
  for (int i = 0; i < 3; i++) {
    dy[2 * i] = creal(p[i]); dy[2 * i + 1] = cimag(p[i]);
    dy[6 + 2 * i] = -creal(F[i]); dy[7 + 2 * i] = -cimag(F[i]);
  }
  double rX2 = creal(c[0] * conj(c[0])), rS2 = creal(c[2] * conj(c[2]));
  dy[12] = cimag(conj(c[0]) * p[0]) / rX2;
  dy[13] = cimag(conj(c[2]) * p[2]) / rS2;
  for (int j = 14; j < NV; j++) dy[j] = 0;
  if (a >= aw0) {
    double N[3]; charges(y, N);
    double E = EX(a, y);
    double Nc = 2 * sqrt(2 * lx) * pow(E / (3 * lx), 0.75);
    dy[14] = N[0]; dy[15] = N[1]; dy[16] = N[2]; dy[17] = E; dy[18] = Nc; dy[19] = N[0] / Nc; dy[20] = thetaG(y);
  }
}

/* Dormand-Prince 5(4) coefficients */
static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
  a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
  a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
  a71 = 35. / 384, a73 = 500. / 1113, a74 = 125. / 192, a75 = -2187. / 6784, a76 = 11. / 84;
static const double e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525,
  e7 = -1. / 40;

static long steps = 0, rejects = 0;

static int integrate(double *t, double tend, double *y, double *h, double rtol, double atol, FILE *hist,
                     double *minS, double *maxGauss, double Nc0, double *F2bins, double *wbins, int nb, double lbmin,
                     double lbmax) {
  double k1[NV], k2[NV], k3[NV], k4[NV], k5[NV], k6[NV], k7[NV], yt[NV], yn[NV];
  rhs(*t, y, k1);
  while (*t < tend) {
    if (*t + *h > tend) *h = tend - *t;
    double hs = *h;
    for (int i = 0; i < NV; i++) yt[i] = y[i] + hs * a21 * k1[i];
    rhs(*t + c2 * hs, yt, k2);
    for (int i = 0; i < NV; i++) yt[i] = y[i] + hs * (a31 * k1[i] + a32 * k2[i]);
    rhs(*t + c3 * hs, yt, k3);
    for (int i = 0; i < NV; i++) yt[i] = y[i] + hs * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    rhs(*t + c4 * hs, yt, k4);
    for (int i = 0; i < NV; i++) yt[i] = y[i] + hs * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    rhs(*t + c5 * hs, yt, k5);
    for (int i = 0; i < NV; i++)
      yt[i] = y[i] + hs * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    rhs(*t + hs, yt, k6);
    for (int i = 0; i < NV; i++)
      yn[i] = y[i] + hs * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    rhs(*t + hs, yn, k7);
    double err = 0;
    for (int i = 0; i < 14; i++) { /* error control on fields, momenta and phases */
      double sc = atol + rtol * fmax(fabs(y[i]), fabs(yn[i]));
      double ei = hs * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]) / sc;
      err += ei * ei;
    }
    err = sqrt(err / 14);
    if (!isfinite(err)) return 2;
    if (err <= 1) {
      *t += hs;
      for (int i = 0; i < NV; i++) { y[i] = yn[i]; k1[i] = k7[i]; }
      steps++;
      double a = *t;
      double rs = sqrt(y[4] * y[4] + y[5] * y[5]) / a;
      double sv = rs / (hh / sqrt(2.0));
      if (a > 1.0 && sv < *minS) *minS = sv;
      double N[3]; charges(y, N);
      double g = fabs(N[0] + 3 * N[1] + 6 * N[2]) / Nc0;
      if (g > *maxGauss) *maxGauss = g;
      if (F2bins) {
        double ax2 = (y[0] * y[0] + y[1] * y[1]), ay2 = (y[2] * y[2] + y[3] * y[3]), as2 = (y[4] * y[4] + y[5] * y[5]);
        double F2 = 2 * (ax2 + 9 * ay2 + 36 * as2) / (a * a);
        int b = (int)floor((log(a) - lbmin) / (lbmax - lbmin) * nb);
        if (b >= 0 && b < nb) { F2bins[b] += F2 * hs; wbins[b] += hs; }
        F2bins[nb + (b >= 0 && b < nb ? b : 0)] += (b >= 0 && b < nb) ? 2 * ax2 / (a * a) * hs : 0;
      }
    } else rejects++;
    double fac = 0.9 * pow(err > 1e-30 ? err : 1e-30, -0.2);
    if (fac > 5) fac = 5;
    if (fac < 0.2) fac = 0.2;
    *h = hs * fac;
    if (*h < 1e-14) return 3;
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 7) { fprintf(stderr, "usage: vrel delta sfac a_end aw0 rtol HI [histfile]\n"); return 1; }
  double delta = atof(argv[1]), sfac = atof(argv[2]), aend = atof(argv[3]);
  aw0 = atof(argv[4]);
  double rtol = atof(argv[5]), HI = atof(argv[6]);
  const char *hf = (argc > 7 && argv[7][0] != '-') ? argv[7] : NULL;
  double xfac = argc > 9 ? atof(argv[8]) : 1.0, yfac = argc > 9 ? atof(argv[9]) : 1.0;
  mu = 5e11 / HI; M = 1e12 / HI; hh = 1e11 / HI;
  /* inflationary minimum at the leading point (vmodel.py, |grad| = 2.9e-15) */
  double X0 = 23.826064298969072, Y0 = 1.8937619219396156, S0 = 0.43086809346988353;
  double y[NV];
  for (int i = 0; i < NV; i++) y[i] = 0;
  double complex cX = X0 * xfac, cY = -Y0 * yfac, cS = -S0 * sfac * cexp(I * delta);
  double complex c0[3] = {cX, cY, cS};
  for (int i = 0; i < 3; i++) {
    y[2 * i] = creal(c0[i]); y[2 * i + 1] = cimag(c0[i]);
    y[6 + 2 * i] = creal(c0[i]); y[7 + 2 * i] = cimag(c0[i]); /* chi' = chi at a = 1 */
  }
  y[12] = 0.0; y[13] = M_PI + delta;
  double thG0 = thetaG(y);
  double E0 = EX(1.0, y);
  double Nc0 = 2 * sqrt(2 * lx) * pow(E0 / (3 * lx), 0.75);
  double t = 1.0, h = 1e-4, minS = 1e300, maxG = 0;
  int nb = 40; double lbmin = 0, lbmax = log(aend);
  double F2b[80] = {0}, wb[40] = {0};
  FILE *fh = hf ? fopen(hf, "w") : NULL;
  /* history output points: log-spaced */
  int nh = hf ? 120 : 0;
  double tpts[200]; int np = 0;
  for (int j = 1; j <= nh; j++) tpts[np++] = exp(log(aend) * j / nh);
  tpts[np++] = aw0;
  /* sort */
  for (int i = 0; i < np; i++) for (int j = i + 1; j < np; j++) if (tpts[j] < tpts[i]) { double z = tpts[i]; tpts[i] = tpts[j]; tpts[j] = z; }
  int rc = 0;
  for (int j = 0; j < np && rc == 0; j++) {
    if (tpts[j] <= t) continue;
    double tend = tpts[j] > aend ? aend : tpts[j];
    rc = integrate(&t, tend, y, &h, rtol, rtol * 1e-3, fh, &minS, &maxG, Nc0, hf ? F2b : NULL, wb, nb, lbmin, lbmax);
    if (fh && rc == 0) {
      double N[3]; charges(y, N);
      double a = t, rs = sqrt(y[4] * y[4] + y[5] * y[5]) / a;
      fprintf(fh, "H %.10g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", a, y[13], thetaG(y), N[0], N[1], N[2],
              rs / (hh / sqrt(2.0)), EX(a, y));
    }
  }
  if (rc == 0 && t < aend) rc = integrate(&t, aend, y, &h, rtol, rtol * 1e-3, fh, &minS, &maxG, Nc0, hf ? F2b : NULL, wb, nb, lbmin, lbmax);
  if (rc) { printf("{\"delta\": %.17g, \"sfac\": %.17g, \"fail\": %d, \"a\": %.6g}\n", delta, sfac, rc, t); return rc; }
  double W = aend - aw0;
  double NXm = y[14] / W, NYm = y[15] / W, NSm = y[16] / W, EXm = y[17] / W, Ncm = y[18] / W, epsm = y[19] / W,
         thGm = y[20] / W;
  double N[3]; charges(y, N);
  double tail = -NSm / (hh * hh * aend);
  if (fh) {
    for (int b = 0; b < nb; b++) if (wb[b] > 0)
      fprintf(fh, "F %.10g %.12g %.12g\n", exp(lbmin + (b + 0.5) * (lbmax - lbmin) / nb), F2b[b] / wb[b], F2b[nb + b] / wb[b]);
    fclose(fh);
  }
  printf("{\"xfac\": %.17g, \"yfac\": %.17g, \"delta\": %.17g, \"sfac\": %.17g, \"a_end\": %.6g, \"rtol\": %.3g, \"NX_mean\": %.12g, \"NY_mean\": %.12g, "
         "\"NS_mean\": %.12g, \"EX_mean\": %.12g, \"eps_ratio\": %.12g, \"eps_mean\": %.12g, \"NX_end\": %.12g, "
         "\"NY_end\": %.12g, \"NS_end\": %.12g, \"thS_end\": %.12g, \"thX_end\": %.12g, \"thG_end\": %.12g, "
         "\"thG_window_mean\": %.12g, \"thG0\": %.15g, \"thS_tail\": %.12g, \"thS_inf\": %.12g, \"minS_over_vac\": %.6g, "
         "\"max_gauss\": %.3g, \"steps\": %ld, \"rejects\": %ld}\n",
         xfac, yfac, delta, sfac, aend, rtol, NXm, NYm, NSm, EXm, NXm / Ncm, epsm, N[0], N[1], N[2], y[13], y[12], thetaG(y), thGm,
         thG0, tail, y[13] + tail, minS, maxG, steps, rejects);
  return 0;
}
