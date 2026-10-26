/* rp_reheating track (Claude, own code; issue #13 review of Codex e160356 handoff item 1).
 *
 * Jordan-frame homogeneous X, Y, Sigma (charges 1,3,6) with the nonminimal term (1/2)(M_r^2 - 2 xi |X|^2) R,
 * a canonical quadratic inflaton I decaying into radiation through a Gamma dI/dt friction term, flat FLRW.
 * Own derivation (NOTES.md section 1):  F G_mn = T_mn + nabla_m nabla_n F - g_mn box F  with F = M_r^2 - 2 xi |X|^2,
 *   Friedmann      3 F H^2 + 3 H Fdot = rho,
 *   Raychaudhuri   Hdot = -(rho + p + Fddot - H Fdot) / (2 F),
 *   trace          (F + 12 xi^2 |X|^2) R = rho - 3p + 12 xi |Xdot|^2 - 12 xi Re(X^* dV/dX^*)   (Xddot eliminated),
 *   fields         zddot + 3 H zdot + dV/dz^* + xi R X delta_zX = 0,   Iddot + (3H + Gamma) Idot + m_I^2 I = 0,
 *   radiation      rhodot_r = Gamma Idot^2 - 4 H rho_r.
 * DIFFERENCE from Codex's engine: H is RECONSTRUCTED from the Friedmann constraint at every evaluation (quadratic in H),
 * not evolved; the Raychaudhuri equation is integrated into an auxiliary variable H_aux and |H_aux - H| is the
 * geometry consistency check (the mirror image of Codex's evolved-H / Friedmann-residual check).
 * Units: H_i = 1, tau = H_i t, fields in H_i, I in M_r, eta2 = (H_i/M_r)^2, rho and rho_r in M_r^2 H_i^2.
 *
 * Modes:  0 physical oscillator inflaton, curvature mass xi R (the Codex model);
 *         1 fluid inflaton (w = 0, decays at Gamma), curvature mass xi R with the fluid trace (R = rho_I/F: the
 *           matter AVERAGE of mode 0, no oscillation);
 *         2 fluid inflaton, curvature mass -cH rho_I/3 (the hw_verify finite-release model, -cH H^2 in matter era);
 *         3 radiation only from tau = 0 (eta2 = 0, rho_r = 3: the instant release / RD reference).
 * State y: 0 ln a | 1 i=I/M_r | 2 j=i' | 3..8 X,Y,S re/im | 9..14 velocities | 15 r=rho_r | 16 H_aux | 17 W_hid |
 *          18 rho_I (fluid modes only).
 * Integrator: Dormand-Prince 5(4) adaptive (coefficients as in our hw_verify/hwrel.c, same author), error per step on
 * all state variables, step clipped to land on sample times.
 * Usage: rpreheat mode b gamma xi eta2 cH lx ly ls kap M h mu end rtol atol nsamp Xre Xim Yre Yim Sre Sim
 * Output: S lines (samples), W lines (late windows), A lines (a-window 30-40 for mode 3), END counters.
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double complex cplx;
#define ND 19
#define NQ 26
static int MODE;
static double b, gam, xi, eta2, cH, lx, ly, ls, kap, MM, hh, mu;

static void potential(const cplx *z, double *V, cplx *F) {
  cplx x = z[0], y = z[1], s = z[2];
  double x2 = creal(x * conj(x)), y2 = creal(y * conj(y)), s2 = creal(s * conj(s));
  *V = lx * x2 * x2 + MM * MM * y2 + ly * y2 * y2 + ls * (s2 - hh * hh / 2) * (s2 - hh * hh / 2)
       + 2 * creal(kap * conj(y) * x * x * x + mu * conj(s) * y * y);
  F[0] = 2 * lx * x2 * x + 3 * kap * y * conj(x) * conj(x);
  F[1] = MM * MM * y + 2 * ly * y2 * y + kap * x * x * x + 2 * mu * s * conj(y);
  F[2] = 2 * ls * (s2 - hh * hh / 2) * s + mu * y * y;
}

/* geometry at a state: returns H (constraint), R, rho, f, fd, curvature mass^2 m2c, and the X acceleration pieces */
typedef struct { double H, R, rho, f, fd, m2c, rhoI, pI, K, V, w; cplx z[3], v[3], F[3]; } geo;

static void geometry(const double *y, geo *g) {
  for (int k = 0; k < 3; k++) { g->z[k] = y[3 + 2 * k] + I * y[4 + 2 * k]; g->v[k] = y[9 + 2 * k] + I * y[10 + 2 * k]; }
  potential(g->z, &g->V, g->F);
  g->K = 0; for (int k = 0; k < 3; k++) g->K += creal(g->v[k] * conj(g->v[k]));
  double x2 = creal(g->z[0] * conj(g->z[0]));
  g->f = 1 - 2 * xi * eta2 * x2;
  g->fd = -4 * xi * eta2 * creal(conj(g->z[0]) * g->v[0]);
  double trI;   /* rho_I - 3 p_I */
  if (MODE == 0) { g->rhoI = 0.5 * y[2] * y[2] + 0.5 * b * b * y[1] * y[1]; g->pI = 0.5 * y[2] * y[2] - 0.5 * b * b * y[1] * y[1]; trI = g->rhoI - 3 * g->pI; }
  else if (MODE == 3) { g->rhoI = 0; g->pI = 0; trI = 0; }
  else { g->rhoI = y[18]; g->pI = 0; trI = g->rhoI; }
  g->rho = g->rhoI + eta2 * (g->K + g->V) + y[15];
  /* Friedmann constraint solved for H (positive root) */
  g->H = (-3 * g->fd + sqrt(9 * g->fd * g->fd + 12 * g->f * g->rho)) / (6 * g->f);
  double vx2 = creal(g->v[0] * conj(g->v[0]));
  g->R = (trI + eta2 * (4 * g->V - 2 * g->K + 12 * xi * vx2 - 12 * xi * creal(conj(g->z[0]) * g->F[0])))
         / (g->f + 12 * xi * xi * eta2 * x2);
  g->m2c = (MODE == 2) ? -cH * g->rhoI / 3.0 : xi * g->R;
  g->w = (g->rhoI > 1e-300) ? g->pI / g->rhoI : 0;
}

static void rhs(const double *y, double *d) {
  geo g; geometry(y, &g);
  double H = g.H;
  d[0] = H;
  if (MODE == 0) { d[1] = y[2]; d[2] = -(3 * H + gam) * y[2] - b * b * y[1]; d[18] = 0; d[15] = gam * y[2] * y[2] - 4 * H * y[15]; }
  else if (MODE == 3) { d[1] = d[2] = d[18] = 0; d[15] = -4 * H * y[15]; }
  else { d[1] = d[2] = 0; d[18] = -(3 * H + gam) * y[18]; d[15] = gam * y[18] - 4 * H * y[15]; }
  cplx acc[3];
  for (int k = 0; k < 3; k++) acc[k] = -3 * H * g.v[k] - g.F[k];
  acc[0] -= g.m2c * g.z[0];
  for (int k = 0; k < 3; k++) { d[3 + 2 * k] = creal(g.v[k]); d[4 + 2 * k] = cimag(g.v[k]); d[9 + 2 * k] = creal(acc[k]); d[10 + 2 * k] = cimag(acc[k]); }
  double vx2 = creal(g.v[0] * conj(g.v[0]));
  double fdd = -4 * xi * eta2 * (vx2 + creal(conj(g.z[0]) * acc[0]));
  double rhoppI = (MODE == 0) ? y[2] * y[2] : g.rhoI;              /* rho_I + p_I */
  if (MODE == 3) rhoppI = 0;
  d[16] = -(rhoppI + 2 * eta2 * g.K + 4 * y[15] / 3 + fdd - H * g.fd) / (2 * g.f);   /* Raychaudhuri for H_aux */
  d[17] = -6 * H * g.K - 2 * g.m2c * creal(conj(g.z[0]) * g.v[0]);               /* d(K+V)/dt, own derivation */
}

static const char *QN[NQ] = {"tau", "a", "H", "R", "NX", "NY", "NS", "Gauss", "EX", "eps", "NYH", "NYL", "Nsig", "absX", "absY",
                             "S_over_vev", "rho_r", "rho_I", "radfrac", "Hcheck", "Wcheck", "w_I", "m2curv", "KV", "Haux", "rho"};
static double KV0 = 0;

static void diag(double t, const double *y, double *q) {
  geo g; geometry(y, &g);
  double a = exp(y[0]), a3 = a * a * a, a4 = a3 * a;
  double n[3]; for (int k = 0; k < 3; k++) n[k] = -2 * cimag(conj(g.z[k]) * g.v[k]) * a3;
  double EX = a4 * (creal(g.v[0] * conj(g.v[0])) + lx * pow(cabs(g.z[0]), 4));
  double Nc = 2 * sqrt(2 * lx) * pow(EX / (3 * lx), 0.75);
  /* free Y quanta about the forced response with instantaneous Sigma (same diagnostic definition as the parent
     packets and our hw tracks; own implementation) */
  double A = MM * MM; cplx B = 2 * mu * g.z[2], J = kap * g.z[0] * g.z[0] * g.z[0];
  double bb = cabs(B), D = A * A - bb * bb, NYH = NAN, NYL = NAN;
  if (D > 0 && bb > 1e-20) {
    cplx yp = (-A * J + B * conj(J)) / D;
    cplx Jd = 3 * kap * g.z[0] * g.z[0] * g.v[0] + 2 * mu * g.v[2] * conj(yp);
    cplx ypd = (-A * Jd + B * conj(Jd)) / D;
    double th = cimag(conj(g.z[2]) * g.v[2]) / creal(g.z[2] * conj(g.z[2]));
    cplx rot = cexp(-0.5 * I * carg(g.z[2]));
    cplx dy = (g.z[1] - yp) * rot, dv = (g.v[1] - ypd - 0.5 * I * th * (g.z[1] - yp)) * rot;
    double wH = sqrt(A + bb), wL = sqrt(A - bb);
    NYH = a3 * (creal(dv) * creal(dv) + wH * wH * creal(dy) * creal(dy)) / wH;
    NYL = a3 * (cimag(dv) * cimag(dv) + wL * wL * cimag(dy) * cimag(dy)) / wL;
  }
  double sr = cabs(g.z[2]), sd = sr > 1e-20 ? creal(conj(g.z[2]) * g.v[2]) / sr : 0;
  double Nsig = a3 * (sd * sd + ls * pow(sr * sr - hh * hh / 2, 2)) / (sqrt(2 * ls) * hh);
  double scale = fabs(n[0]) + 3 * fabs(n[1]) + 6 * fabs(n[2]); if (scale < 1) scale = 1;
  q[0] = t; q[1] = a; q[2] = g.H; q[3] = g.R; q[4] = n[0]; q[5] = n[1]; q[6] = n[2]; q[7] = (n[0] + 3 * n[1] + 6 * n[2]) / scale;
  q[8] = EX; q[9] = n[0] / Nc; q[10] = NYH; q[11] = NYL; q[12] = Nsig; q[13] = cabs(g.z[0]); q[14] = cabs(g.z[1]);
  q[15] = sr / (hh / sqrt(2)); q[16] = y[15]; q[17] = g.rhoI; q[18] = y[15] / g.rho;
  q[19] = (y[16] - g.H) / g.H; q[20] = (g.K + g.V - KV0 - y[17]) / KV0; q[21] = g.w; q[22] = g.m2c; q[23] = g.K + g.V;
  q[24] = y[16]; q[25] = g.rho;
}

int main(int argc, char **argv) {
  if (argc != 24) { fprintf(stderr, "usage: see header (%d args given)\n", argc); return 1; }
  int k = 1;
  MODE = atoi(argv[k++]); b = atof(argv[k++]); gam = atof(argv[k++]); xi = atof(argv[k++]); eta2 = atof(argv[k++]); cH = atof(argv[k++]);
  lx = atof(argv[k++]); ly = atof(argv[k++]); ls = atof(argv[k++]); kap = atof(argv[k++]); MM = atof(argv[k++]); hh = atof(argv[k++]); mu = atof(argv[k++]);
  double end = atof(argv[k++]), rtol = atof(argv[k++]), atol = atof(argv[k++]); int ns = atoi(argv[k++]);
  double y[ND] = {0};
  for (int i = 0; i < 6; i++) y[3 + i] = atof(argv[k++]);
  /* initial state: a = 1, zero velocities, H = 1 fixed by the inflaton amplitude (or fluid density / radiation) */
  { cplx z[3]; double V; cplx F[3]; for (int i = 0; i < 3; i++) z[i] = y[3 + 2 * i] + I * y[4 + 2 * i]; potential(z, &V, F);
    double f = 1 - 2 * xi * eta2 * creal(z[0] * conj(z[0]));
    double need = 3 * f - eta2 * V;                 /* rho_I (or rho_r) such that H = 1 with fd = 0 */
    if (MODE == 0) y[1] = sqrt(2 * need) / b;
    else if (MODE == 3) y[15] = need;
    else y[18] = need;
    y[16] = 1.0; y[17] = 0.0; KV0 = 0; double q0[NQ]; geo g; geometry(y, &g); KV0 = g.K + g.V; diag(0, y, q0); }
  static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
    a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
    a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
    b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
    e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
  (void)c2; (void)c3; (void)c4; (void)c5;
  double t = 0, h = 1e-4, q[NQ], qo[NQ];
  diag(0, y, q);
  printf("S"); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n");
  long steps = 0, rejects = 0; double maxG = 0, maxH = 0, maxW = 0;
  /* windows: tau-weighted means over [0.8,1] and [0.9,1] of end; a-window [30,40] (mode 3) */
  double s80[NQ] = {0}, s90[NQ] = {0}, sa[NQ] = {0}, w80 = 0, w90 = 0, wa = 0;
  int si = 1; double hmax_cap = 0.05 / fmax(1.0, b / 2.0);
  while (t < end * (1 - 1e-14)) {
    double target = end * si / ns; if (target > end) target = end;
    double hstep = h; if (hstep > hmax_cap) hstep = hmax_cap; if (t + hstep > target) hstep = target - t;
    double k1[ND], k2[ND], k3[ND], k4[ND], k5[ND], k6[ND], k7[ND], yt[ND], yn[ND];
    rhs(y, k1);
    for (int i = 0; i < ND; i++) yt[i] = y[i] + hstep * a21 * k1[i]; rhs(yt, k2);
    for (int i = 0; i < ND; i++) yt[i] = y[i] + hstep * (a31 * k1[i] + a32 * k2[i]); rhs(yt, k3);
    for (int i = 0; i < ND; i++) yt[i] = y[i] + hstep * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]); rhs(yt, k4);
    for (int i = 0; i < ND; i++) yt[i] = y[i] + hstep * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]); rhs(yt, k5);
    for (int i = 0; i < ND; i++) yt[i] = y[i] + hstep * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]); rhs(yt, k6);
    for (int i = 0; i < ND; i++) yn[i] = y[i] + hstep * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    rhs(yn, k7);
    double err = 0;
    for (int i = 0; i < ND; i++) {
      double e = hstep * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
      double sc = atol + rtol * fmax(fabs(y[i]), fabs(yn[i]));
      err += (e / sc) * (e / sc);
    }
    err = sqrt(err / ND);
    if (!isfinite(err)) { fprintf(stderr, "nonfinite at tau=%g\n", t); return 3; }
    double fac = err > 0 ? 0.9 * pow(err, -0.2) : 5; if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
    if (err <= 1) {
      double told = t; for (int j = 0; j < NQ; j++) qo[j] = q[j];
      memcpy(y, yn, sizeof y); t += hstep; steps++;
      diag(t, y, q);
      if (fabs(q[7]) > maxG) maxG = fabs(q[7]);
      if (fabs(q[19]) > maxH) maxH = fabs(q[19]);
      if (fabs(q[20]) > maxW) maxW = fabs(q[20]);
      double dt = t - told;
      if (told >= 0.8 * end * (1 - 1e-12)) { w80 += dt; for (int j = 0; j < NQ; j++) s80[j] += 0.5 * dt * (qo[j] + q[j]); }
      if (told >= 0.9 * end * (1 - 1e-12)) { w90 += dt; for (int j = 0; j < NQ; j++) s90[j] += 0.5 * dt * (qo[j] + q[j]); }
      if (qo[1] >= 30 && q[1] <= 40 + 1e-9) { double da = q[1] - qo[1]; wa += da; for (int j = 0; j < NQ; j++) sa[j] += 0.5 * da * (qo[j] + q[j]); }
      if (fabs(t - target) < 1e-12 * end) {
        printf("S"); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n"); si++;
      }
      if (hstep >= h * (1 - 1e-12)) h = hstep * fac; else h = fmax(h, hstep * fac);
    } else { rejects++; h = hstep * fac; }
    if (steps > 200000000L) { fprintf(stderr, "step limit\n"); return 4; }
  }
  if (w80 > 0) { printf("W80 %.17g", w80); for (int j = 0; j < NQ; j++) printf(" %.17g", s80[j] / w80); printf("\n"); }
  if (w90 > 0) { printf("W90 %.17g", w90); for (int j = 0; j < NQ; j++) printf(" %.17g", s90[j] / w90); printf("\n"); }
  if (wa > 0) { printf("A3040 %.17g", wa); for (int j = 0; j < NQ; j++) printf(" %.17g", sa[j] / wa); printf("\n"); }
  printf("Y"); for (int i = 0; i < ND; i++) printf(" %.17g", y[i]); printf("\n");
  printf("END %ld %ld %.6g %.6g %.6g\n", steps, rejects, maxG, maxH, maxW);
  (void)QN;
  return 0;
}
