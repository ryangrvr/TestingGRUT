/* dq_entropy SKEPTIC (Claude, own code, written fresh, 2 November 2026): an independent late-time integrator for the
 * X <-> Sigma exchange of the M1r release, to test dq_entropy item 1 (frozen N_X, E_X) with a different method.
 *
 * Differences from rpdq.c (a DP5(4) adaptive integrator in physical time with the inflaton, nonminimal gravity and
 * Friedmann back-reaction):
 *   - conformal fields chi_k = a z_k in conformal time eta, exact radiation background a(eta) = a0 + K (eta - eta0)
 *     (a' = a^2 H = K; inflaton decayed, hidden back-reaction ~1e-11 dropped; xi_X R |X|^2 dropped since R ~ eta2 ~ 1e-15);
 *   - equations chi'' = -a^3 F(chi/a), F = dV/dz* of the late potential (u = mu, no portal: du = 0 when u = mu);
 *   - 4th-order Yoshida composition of the velocity-Verlet (symplectic in the extended phase space; it preserves the
 *     U(1) charge N_X + 3 N_Y + 6 N_Sigma to round-off), step h = c_h /(omega_max a) with omega_max = 40;
 *   - diagnostics: N_k = -2 Im(conj(chi_k) chi_k'), the exact conformal X energy |chi_X'|^2 + lX |chi_X|^4, and the
 *     engine's definition E_X = |chi_X' - (K/a) chi_X|^2 + lX |chi_X|^4, all averaged over conformal windows.
 * Input (stdin, one line): a0 eta0 K lx ly ls kap M h mu  Xr Xi Yr Yi Sr Si  Vxr Vxi Vyr Vyi Vsr Vsi
 *   (fields z_k and physical-time velocities zdot_k from rpdq's 'Y' line, H at that time folded into K = a0^2 H).
 * Args: eta_end c_h win1 win2 [eta_s deta_r]   (win1: wide window width; win2: fine bin width; both in eta, aligned to
 *   multiples).  Optional dephasing toy: kappa -> kappa cos^2(pi/2 (eta - eta_s)/deta_r) for eta_s < eta < eta_s + deta_r and
 *   0 afterwards (the only X-number-violating coupling is switched off, so N_X is frozen at the end of the ramp).
 * Output lines: "W eta_lo eta_hi a_mid <NX> <NY> <NS> <Econf> <Eeng> std(NX over fine bins)" and "F ..." fine bins.
 */
#define _DEFAULT_SOURCE
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
typedef double complex cplx;
static double lx, ly, ls, kap, kap0, MM, hh, mu, a0, eta0, K, eta_s = 1e300, deta_r = 1.0;

static void force(const cplx *c, double a, cplx *acc) {
  /* acc_k = -a^3 F_k(chi/a) */
  cplx x = c[0] / a, y = c[1] / a, s = c[2] / a;
  double x2 = creal(x * conj(x)), y2 = creal(y * conj(y)), s2 = creal(s * conj(s));
  cplx F0 = 2 * lx * x2 * x + 3 * kap * y * conj(x) * conj(x);
  cplx F1 = MM * MM * y + 2 * ly * y2 * y + kap * x * x * x + 2 * mu * s * conj(y);
  cplx F2 = 2 * ls * (s2 - hh * hh / 2) * s + mu * y * y;
  double a3 = a * a * a;
  acc[0] = -a3 * F0; acc[1] = -a3 * F1; acc[2] = -a3 * F2;
}

static void vv(cplx *c, cplx *w, double *eta, double h) {
  for (int k = 0; k < 3; k++) c[k] += 0.5 * h * w[k];
  *eta += 0.5 * h;
  double a = a0 + K * (*eta - eta0);
  if (*eta > eta_s) { double u = (*eta - eta_s) / deta_r; kap = (u >= 1) ? 0.0 : kap0 * pow(cos(0.5 * M_PI * u), 2); }
  cplx acc[3]; force(c, a, acc);
  for (int k = 0; k < 3; k++) w[k] += h * acc[k];
  for (int k = 0; k < 3; k++) c[k] += 0.5 * h * w[k];
  *eta += 0.5 * h;
}

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "args: eta_end c_h win1 win2 < state\n"); return 1; }
  double eta_end = atof(argv[1]), ch = atof(argv[2]), win1 = atof(argv[3]), win2 = atof(argv[4]);
  double zr[6], vr[6];
  if (argc >= 7) { eta_s = atof(argv[5]); deta_r = atof(argv[6]); }
  if (scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &a0, &eta0, &K, &lx, &ly, &ls, &kap, &MM, &hh, &mu) != 10) return 2;
  for (int i = 0; i < 6; i++) if (scanf("%lf", &zr[i]) != 1) return 2;
  for (int i = 0; i < 6; i++) if (scanf("%lf", &vr[i]) != 1) return 2;
  kap0 = kap;
  cplx c[3], w[3];
  double H0 = K / (a0 * a0);
  for (int k = 0; k < 3; k++) {
    cplx z = zr[2 * k] + I * zr[2 * k + 1], zd = vr[2 * k] + I * vr[2 * k + 1];
    c[k] = a0 * z; w[k] = a0 * a0 * (zd + H0 * z);
  }
  const double cr = cbrt(2.0), w1 = 1.0 / (2.0 - cr), w0 = -cr / (2.0 - cr);
  double eta = eta0;
  /* window accumulators (trapezoid in eta) */
  double s1[6] = {0}, s2[6] = {0}, d1 = 0, d2 = 0, ssq = 0; int nf = 0;
  double lo1 = floor(eta0 / win1 + 1e-12) * win1 + win1, lo2 = floor(eta0 / win2 + 1e-12) * win2 + win2;
  int on1 = 0, on2 = 0;
  double q[6], qo[6], gmax = 0, gscale = 0; long steps = 0;
  #define DIAG(QQ) do { double a_ = a0 + K * (eta - eta0); \
      double n0 = -2 * cimag(conj(c[0]) * w[0]), n1 = -2 * cimag(conj(c[1]) * w[1]), n2 = -2 * cimag(conj(c[2]) * w[2]); \
      double ec = creal(w[0] * conj(w[0])) + lx * pow(cabs(c[0]), 4); \
      cplx we = w[0] - (K / a_) * c[0]; double ee = creal(we * conj(we)) + lx * pow(cabs(c[0]), 4); \
      QQ[0] = n0; QQ[1] = n1; QQ[2] = n2; QQ[3] = ec; QQ[4] = ee; QQ[5] = a_; } while (0)
  DIAG(q);
  double g0 = q[0] + 3 * q[1] + 6 * q[2]; gscale = fabs(q[0]) + 3 * fabs(q[1]) + 6 * fabs(q[2]);
  printf("I a0 %.17g eta0 %.17g K %.17g NX0 %.17g Econf0 %.17g Eeng0 %.17g gauss0 %.6e\n", a0, eta0, K, q[0], q[3], q[4], g0 / gscale);
  while (eta < eta_end) {
    double a = a0 + K * (eta - eta0);
    double h = ch / (40.0 * a);
    /* do not straddle a fine-bin edge: shorten the step to land on it */
    double nxt = lo2; if (eta + h > nxt) h = nxt - eta; if (h < 1e-15) h = 1e-15;
    for (int j = 0; j < 6; j++) qo[j] = q[j];
    double et0 = eta;
    vv(c, w, &eta, w1 * h); vv(c, w, &eta, w0 * h); vv(c, w, &eta, w1 * h);
    eta = et0 + h;                                   /* remove round-off drift of the composed time */
    steps++;
    DIAG(q);
    double g = fabs((q[0] + 3 * q[1] + 6 * q[2]) - g0) / gscale; if (g > gmax) gmax = g;
    double dt = h;
    if (et0 >= lo1 - win1 - 1e-12) { on1 = 1; }
    if (et0 >= lo2 - win2 - 1e-12) { on2 = 1; }
    if (on1) { d1 += dt; for (int j = 0; j < 6; j++) s1[j] += 0.5 * dt * (qo[j] + q[j]); }
    if (on2) { d2 += dt; for (int j = 0; j < 6; j++) s2[j] += 0.5 * dt * (qo[j] + q[j]); }
    if (eta >= lo2 - 1e-12) {
      if (on2 && d2 > 0.999 * win2) {
        printf("F %.10g %.10g", lo2 - win2, lo2); for (int j = 0; j < 6; j++) printf(" %.17g", s2[j] / d2); printf("\n");
        double m = s2[0] / d2; if (on1) { ssq += m * m; nf++; }
      }
      for (int j = 0; j < 6; j++) s2[j] = 0; d2 = 0; lo2 += win2;
      if (eta >= lo1 - 1e-12) {
        if (on1 && d1 > 0.999 * win1) {
          printf("W %.10g %.10g", lo1 - win1, lo1); for (int j = 0; j < 6; j++) printf(" %.17g", s1[j] / d1);
          double m = s1[0] / d1; printf(" %.6g\n", nf > 1 ? sqrt(fmax(ssq / nf - m * m, 0)) : 0.0);
        }
        for (int j = 0; j < 6; j++) s1[j] = 0; d1 = 0; ssq = 0; nf = 0; lo1 += win1;
      }
    }
  }
  printf("END steps %ld eta %.17g a %.17g gauss_drift_max %.3e\n", steps, eta, a0 + K * (eta - eta0), gmax);
  printf("Z"); for (int k = 0; k < 3; k++) printf(" %.17g %.17g %.17g %.17g", creal(c[k]), cimag(c[k]), creal(w[k]), cimag(w[k])); printf("\n");
  return 0;
}
