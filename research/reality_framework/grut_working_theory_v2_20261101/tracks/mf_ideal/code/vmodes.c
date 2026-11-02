/* inf_chain track (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): longitudinal hidden-vector modes
 * (Stueckelberg / Goldstone-gauge angle theta_k) in cosmic time on a tabulated background.
 *
 *   S = 1/2 int dt [ p theta_dot^2 - q theta^2 ],  p = a^3 F^2 k^2/(k^2 + a^2 e^2 F^2),  q = a F^2 k^2
 *   (the single-Higgs-direction Proca action of dq_heat / its skeptic's sk_modes.py, A_0 integrated out).
 *   Variables (theta, pi = p theta_dot), complex; theta' = pi/p, pi' = -q theta.
 * Background: t < 0 de Sitter, H = 1, a = e^t, F = FS (the prepared state at the slice); t >= 0 the table
 *   (t, ln a, H, F^2), linear interpolation (binary search), constant extrapolation is an error.
 * Start: mode "D": exact massless de Sitter Goldstone u = a F theta = e^{-ik eta}(1 - i/(k eta))/sqrt(2k), eta = -1/a, at
 *   k/aH = KS (default 50);  mode "F": frozen angle theta = 1, pi = 0 at t = TF (per unit amplitude).
 * Read-out: adiabatic number  n = (|pi|^2/p + q|theta|^2)/(2 Omega) [- 1/2 for "D"], Omega^2 = q/p, at the first time
 *   t > TR with Omega/H > AD, and again at Omega/H > 2 AD (convergence check).
 * Usage: vmodes TABLE FS TR AD MODE TF k1 e1 [k2 e2 ...]     prints "M k e n1 n2 t1 t2 steps" per mode.
 * Units: H_i = 1 at the slice, a = 1 at the slice.  Labels: [C].
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int NT; static double *TT, *LA, *HH, *F2;
static double FS = 276.0, KS = 50.0;

static void bgat(double t, double *a, double *H, double *f2) {
  if (t < TT[0]) { *a = exp(t); *H = 1.0; *f2 = FS * FS; return; }
  if (t >= TT[NT - 1]) { fprintf(stderr, "t=%g beyond table end %g\n", t, TT[NT - 1]); exit(3); }
  int lo = 0, hi = NT - 1;
  while (hi - lo > 1) { int m = (lo + hi) / 2; if (TT[m] <= t) lo = m; else hi = m; }
  double w = (t - TT[lo]) / (TT[hi] - TT[lo]);
  *a = exp(LA[lo] + w * (LA[hi] - LA[lo])); *H = HH[lo] + w * (HH[hi] - HH[lo]); *f2 = F2[lo] + w * (F2[hi] - F2[lo]);
}

static double K, E2;
static void pq(double t, double *p, double *q, double *Om, double *H) {
  double a, f2; bgat(t, &a, H, &f2);
  double k2 = K * K;
  *p = a * a * a * f2 * k2 / (k2 + a * a * E2 * f2);
  *q = a * f2 * k2;
  *Om = sqrt(*q / *p);
}

static void rhs(double t, const double *y, double *d) {
  double p, q, Om, H; pq(t, &p, &q, &Om, &H);
  d[0] = y[2] / p; d[1] = y[3] / p; d[2] = -q * y[0]; d[3] = -q * y[1];
}

static double number(double t, const double *y, int sub) {
  double p, q, Om, H; pq(t, &p, &q, &Om, &H);
  double n = ((y[2] * y[2] + y[3] * y[3]) / p + q * (y[0] * y[0] + y[1] * y[1])) / (2 * Om);
  return sub ? n - 0.5 : n;
}

int main(int argc, char **argv) {
  if (argc < 9) { fprintf(stderr, "usage: vmodes TABLE FS TR AD MODE TF k e ...\n"); return 1; }
  FILE *f = fopen(argv[1], "r"); if (!f) { perror("table"); return 1; }
  if (fscanf(f, "%d", &NT) != 1) return 1;
  TT = malloc(NT * sizeof(double)); LA = malloc(NT * sizeof(double)); HH = malloc(NT * sizeof(double)); F2 = malloc(NT * sizeof(double));
  for (int i = 0; i < NT; i++) if (fscanf(f, "%lf %lf %lf %lf", &TT[i], &LA[i], &HH[i], &F2[i]) != 4) { fprintf(stderr, "bad table row %d\n", i); return 1; }
  fclose(f);
  FS = atof(argv[2]); double TR = atof(argv[3]), AD = atof(argv[4]); char MODE = argv[5][0]; double TF = atof(argv[6]);
  if (getenv("VM_KS")) KS = atof(getenv("VM_KS"));
  double rtol = getenv("VM_RTOL") ? atof(getenv("VM_RTOL")) : 1e-9;
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
    a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
    a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
    b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
    e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
  for (int m = 7; m + 1 < argc; m += 2) {
    K = atof(argv[m]); double e = atof(argv[m + 1]); E2 = e * e;
    double y[4], t;
    int sub = 0;
    if (MODE == 'D') {
      t = log(K / KS); if (t > -1e-12) { fprintf(stderr, "k too large for dS start\n"); return 2; }
      double a = exp(t), eta = -1.0 / a, ke = K * eta;
      /* u = e^{-ik eta}(1 - i/(k eta))/sqrt(2k); du/deta = e^{-ik eta}[-ik(1 - i/(k eta)) + i/(k eta^2)]/sqrt(2k) */
      double c = cos(ke), s = -sin(ke), n2 = 1.0 / sqrt(2 * K);
      double ur = (c * 1 - s * (-1.0 / ke)) * n2, ui = (c * (-1.0 / ke) + s * 1) * n2;
      double wr = -K * (1.0 / ke) * 0 + 0;  (void)wr;
      /* g = -ik(1 - i/(k eta)) + i/(k eta^2) = -ik - 1/eta + i/(k eta^2) */
      double gr = -1.0 / eta, gi = -K + 1.0 / (K * eta * eta);
      double dur = (c * gr - s * gi) * n2, dui = (c * gi + s * gr) * n2;
      double F = FS, H = 1.0;
      /* theta = u/(aF); theta_dot = (u_dot - H u)/(aF), u_dot = (du/deta)/a */
      double thr = ur / (a * F), thi = ui / (a * F);
      double tdr = (dur / a - H * ur) / (a * F), tdi = (dui / a - H * ui) / (a * F);
      double p, q, Om, Hh; pq(t, &p, &q, &Om, &Hh);
      y[0] = thr; y[1] = thi; y[2] = p * tdr; y[3] = p * tdi; sub = 1;
    } else { t = TF; y[0] = 1; y[1] = 0; y[2] = 0; y[3] = 0; }
    double h = 1e-3 * fmin(1.0, 1.0 / K), n1 = NAN, n2 = NAN, t1 = NAN, t2 = NAN; long steps = 0;
    double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4], yt[4], yn[4];
    while (1) {
      double p, q, Om, H; pq(t, &p, &q, &Om, &H);
      if (t > TR && Om / H > AD && isnan(n1)) { n1 = number(t, y, sub); t1 = t; }
      if (t > TR && Om / H > 2 * AD) { n2 = number(t, y, sub); t2 = t; break; }
      double hmax = 0.1 / Om; if (h > hmax) h = hmax;
      rhs(t, y, k1);
      for (int i = 0; i < 4; i++) yt[i] = y[i] + h * a21 * k1[i]; rhs(t + h / 5, yt, k2);
      for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]); rhs(t + 3 * h / 10, yt, k3);
      for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]); rhs(t + 4 * h / 5, yt, k4);
      for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]); rhs(t + 8 * h / 9, yt, k5);
      for (int i = 0; i < 4; i++) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]); rhs(t + h, yt, k6);
      for (int i = 0; i < 4; i++) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
      rhs(t + h, yn, k7);
      double err = 0, scale = 0;
      for (int i = 0; i < 4; i++) scale += yn[i] * yn[i];
      /* error control relative to the state norm in each pair (theta, pi) separately */
      double sth = sqrt(y[0] * y[0] + y[1] * y[1] + yn[0] * yn[0] + yn[1] * yn[1]) + 1e-300;
      double spi = sqrt(y[2] * y[2] + y[3] * y[3] + yn[2] * yn[2] + yn[3] * yn[3]);
      double spi_ref = sqrt(p * q) * sth;                    /* natural pi scale */
      if (spi < spi_ref) spi = spi_ref;
      for (int i = 0; i < 4; i++) {
        double ee = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
        double sc = rtol * (i < 2 ? sth : spi);
        err += (ee / sc) * (ee / sc);
      }
      err = sqrt(err / 4); (void)scale;
      double fac = err > 0 ? 0.9 * pow(err, -0.2) : 5; if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
      if (err <= 1) { for (int i = 0; i < 4; i++) y[i] = yn[i]; t += h; steps++; }
      h *= fac;
      if (steps > 200000000L) { fprintf(stderr, "step limit k=%g\n", K); break; }
    }
    printf("M %.10g %.6g %.12g %.12g %.8g %.8g %ld\n", K, e, n1, n2, t1, t2, steps);
    fflush(stdout);
  }
  return 0;
}
