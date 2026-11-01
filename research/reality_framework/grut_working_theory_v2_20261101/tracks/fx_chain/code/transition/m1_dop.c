/* m1_dop.c -- track pt_m1tilt (Claude review, issue #13): the phi-condensation transition under Codex's GRUT-M1
 * potential (commit d4747cb).  Copied from the pt_matching skeptic's sk_dop.c (own adaptive DOP853, Hairer
 * coefficients from scipy, dop853_coef.h) and extended by one potential mode.
 * Homogeneous conformal X (two real fields x1, x2) + phi, conformal variables chi = a phi, M = a v:
 *   mode 1 (skeptic model):  V = lamX/4 r^4 + g2/2 chi^2 r^2 + lamP/4 chi^4 - lamP/2 M^2 chi^2 - beta lamP M^3 chi
 *   mode 2:                  X mass imposed as g2 M^2 r^2/2 (phi absent)
 *   mode 3 (GRUT-M1):        V = lamX/4 r^4 + g2/2 chi^2 r^2 + lamP/4 (chi^2 - M^2)^2 - epst M (chi - M)^3
 *     force on chi:  -g2 r^2 chi - lamP (chi^2 - M^2) chi + 3 epst M (chi - M)^2
 *     explicit dV/ds = M' [ -lamP M (chi^2 - M^2) - epst (chi - M)^3 + 3 epst M (chi - M)^2 ]
 *     epst = 0 is the skeptic's beta = 0 model up to the constant lamP/4 M^4 (carried consistently in E and W).
 *   M(s) = M0 exp(kappa s) (growth 0) or M0 (1 + kappa s) (growth 1).
 * State y = {x1, x2, chi, p1, p2, pchi, W} with W' = dV/ds|_explicit, so E(s) - E(0) - W(s) = 0 checks the run.
 * No Codex or pt_matching code is included or called. */
#include <math.h>
#include <string.h>
#include "dop853_coef.h"

#define NY 7
typedef struct { int mode, growth; double lamX, g2, lamP, M0, kappa, beta, epst; } SPar;

static void Mf(const SPar *P, double s, double *M, double *Mp) {
  if (P->growth == 0) { *M = P->M0 * exp(P->kappa * s); *Mp = P->kappa * *M; }
  else { *M = P->M0 * (1 + P->kappa * s); *Mp = P->M0 * P->kappa; }
}

static void rhs(const SPar *P, double s, const double *y, double *f) {
  double M, Mp; Mf(P, s, &M, &Mp);
  double r2 = y[0] * y[0] + y[1] * y[1], ph = y[2];
  f[0] = y[3]; f[1] = y[4]; f[2] = y[5];
  if (P->mode == 2) {
    double cx = P->lamX * r2 + P->g2 * M * M;
    f[3] = -cx * y[0]; f[4] = -cx * y[1]; f[5] = 0;
    f[6] = P->g2 * M * Mp * r2;
    return;
  }
  double cx = P->lamX * r2 + P->g2 * ph * ph;
  f[3] = -cx * y[0]; f[4] = -cx * y[1];
  if (P->mode == 3) {
    double d = ph - M, u = ph * ph - M * M;
    f[5] = -P->g2 * r2 * ph - P->lamP * u * ph + 3 * P->epst * M * d * d;
    f[6] = Mp * (-P->lamP * M * u - P->epst * d * d * d + 3 * P->epst * M * d * d);
    return;
  }
  f[5] = -(P->g2 * r2 + P->lamP * ph * ph - P->lamP * M * M) * ph + P->beta * P->lamP * M * M * M;
  f[6] = -P->lamP * M * Mp * ph * ph - 3 * P->beta * P->lamP * M * M * Mp * ph;
}

double m1_energy(const SPar *P, double s, const double *y) {
  double M, Mp; Mf(P, s, &M, &Mp);
  double r2 = y[0] * y[0] + y[1] * y[1], ph = y[2];
  if (P->mode == 2) return 0.5 * (y[3] * y[3] + y[4] * y[4]) + 0.25 * P->lamX * r2 * r2 + 0.5 * P->g2 * M * M * r2;
  double K = 0.5 * (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) + 0.25 * P->lamX * r2 * r2 + 0.5 * P->g2 * ph * ph * r2;
  if (P->mode == 3) {
    double d = ph - M, u = ph * ph - M * M;
    return K + 0.25 * P->lamP * u * u - P->epst * M * d * d * d;
  }
  return K + 0.25 * P->lamP * ph * ph * ph * ph - 0.5 * P->lamP * M * M * ph * ph - P->beta * P->lamP * M * M * M * ph;
}

/* Integrate from *s to s1 with adaptive DOP853.  *h carries the step guess between calls.
 * Returns number of accepted steps; *nrej counts rejected steps. */
long m1_run(const SPar *P, double *y, double *s, double s1, double rtol, double atol, double *h, long *nrej) {
  double K[13][NY], yt[NY], yn[NY], e5[NY], e3[NY];
  long nacc = 0;
  int i, j, c;
  if (*s >= s1) return 0;
  rhs(P, *s, y, K[0]);
  while (*s < s1) {
    double hh = *h;
    int last = 0;
    if (*s + hh >= s1) { hh = s1 - *s; last = 1; }
    for (i = 1; i < NST; i++) {
      for (c = 0; c < NY; c++) {
        double acc = 0;
        for (j = 0; j < i; j++) acc += DA[i * 13 + j] * K[j][c];
        yt[c] = y[c] + hh * acc;
      }
      rhs(P, *s + DC[i] * hh, yt, K[i]);
    }
    for (c = 0; c < NY; c++) {
      double acc = 0;
      for (j = 0; j < NST; j++) acc += DB[j] * K[j][c];
      yn[c] = y[c] + hh * acc;
    }
    rhs(P, *s + hh, yn, K[12]);
    double n5 = 0, n3 = 0;
    for (c = 0; c < NY; c++) {
      double a5 = 0, a3 = 0;
      for (j = 0; j < 13; j++) { a5 += DE5[j] * K[j][c]; a3 += DE3[j] * K[j][c]; }
      double sc = atol + rtol * fmax(fabs(y[c]), fabs(yn[c]));
      if (c == 6) sc = atol + rtol * fmax(1.0, fmax(fabs(y[c]), fabs(yn[c]))) * 1e6; /* work: not step-controlling */
      e5[c] = a5 / sc; e3[c] = a3 / sc;
      n5 += e5[c] * e5[c]; n3 += e3[c] * e3[c];
    }
    double err = 0;
    if (n5 > 0 || n3 > 0) err = fabs(hh) * n5 / sqrt((n5 + 0.01 * n3) * NY);
    if (err < 1.0) {
      *s = last ? s1 : *s + hh;
      memcpy(y, yn, sizeof yn);
      memcpy(K[0], K[12], sizeof K[0]);
      nacc++;
      double fac = err == 0 ? 10.0 : fmin(10.0, 0.9 * pow(err, -1.0 / 8.0));
      if (!last) *h = hh * fac; else *h = fmax(*h, hh);
    } else {
      double fac = fmax(0.2, 0.9 * pow(err, -1.0 / 8.0));
      *h = hh * fac;
      (*nrej)++;
    }
  }
  return nacc;
}
