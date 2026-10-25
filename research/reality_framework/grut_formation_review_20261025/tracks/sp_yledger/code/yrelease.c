/* sp_yledger track (Claude, own code): homogeneous X, Y, Sigma release in COSMIC time with PHYSICAL
 * fields (Codex 373d37e uses conformal fields; this is a separate formulation), with three backgrounds:
 *   bg 0: instant radiation release, H = 1/(2t), a = sqrt(2t), t0 = 1/2 (H(t0) = 1, a(t0) = 1).
 *   bg 1: inflaton matter era with perturbative decay Gamma_phi (H_i units); X Hubble mass removed at t0.
 *   bg 2: as bg 1 but the X Hubble mass persists as -cH (rho_phi/3)|X|^2 (inflaton-induced).
 * Units: H_i = 1, rho in units of M_red^2 H_i^2 for the background, fields in H_i units.
 * V = lx|X|^4 + ly|Y|^4 + ls(|S|^2-h^2/2)^2 + M^2|Y|^2 + 2Re(kap Y* X^3) + 2Re(mu S* Y^2)
 * EOM: z'' + 3H z' + dV/dz* = 0 (dot = d/dt).
 * Integrator: Dormand-Prince 5(4) with adaptive steps (own implementation, reused from our gp_phase track).
 * Diagnostics (all comoving, H_i^3 units): N_i = -2 a^3 Im(z* z'), E_X = a^4(|X'|^2 + lx|X|^4),
 *   eps = N_X / Ncirc(E_X); Y excitation about the static forced solution, projected onto the two real
 *   Majorana mass eigenmodes (adiabatic numbers N_YH, N_YL), plus the Codex-style single-mass number J_Y;
 *   Sigma radial excitation number N_sig = a^3 (r'^2 + ls (r^2-h^2/2)^2)/m_sig.
 * Energy control: W integrates d(a^3 rho)/dt = 3 H a^3 (V - K) - a^3 dV/dt|explicit; compared with a^3 rho.
 * Usage: yrelease bg Gphi cH M h mu t_end rtol atol nsamp x0r x0i y0r y0i s0r s0i
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

static int BG; static double GPHI, CH;
static const double lx = 1e-3, ly = 1.0, ls = 1.0, kap = 1e-3;
static double MM, hh, mu;

/* state: 0..5 z (re,im)x3, 6..11 z' (re,im)x3, 12 ln a, 13 Rphi = rho_phi a^3, 14 Rr = rho_r a^4, 15 W */
#define N 16

static void background(double t, const double *y, double *a, double *H, double *rphi) {
  if (BG == 0) { *a = sqrt(2 * t); *H = 1.0 / (2 * t); *rphi = 0; return; }
  double A = exp(y[12]);
  double rp = y[13] / (A * A * A), rr = y[14] / (A * A * A * A);
  *a = A; *H = sqrt((rp + rr) / 3.0); *rphi = rp;
}

static void force(const double complex *z, double complex *F) {
  double complex X = z[0], Y = z[1], S = z[2];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  F[0] = 2 * lx * ax2 * X + 3 * kap * Y * conj(X) * conj(X);
  F[1] = MM * MM * Y + 2 * ly * ay2 * Y + kap * X * X * X + 2 * mu * S * conj(Y);
  F[2] = 2 * ls * (as2 - hh * hh / 2) * S + mu * Y * Y;
}

static double potential(const double complex *z) {
  double complex X = z[0], Y = z[1], S = z[2];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  return lx * ax2 * ax2 + ly * ay2 * ay2 + ls * (as2 - hh * hh / 2) * (as2 - hh * hh / 2) + MM * MM * ay2
       + 2 * creal(kap * conj(Y) * X * X * X) + 2 * creal(mu * conj(S) * Y * Y);
}

static void rhs(double t, const double *y, double *dy) {
  double a, H, rphi; background(t, y, &a, &H, &rphi);
  double complex z[3], zd[3], F[3];
  for (int i = 0; i < 3; i++) { z[i] = y[2 * i] + I * y[2 * i + 1]; zd[i] = y[6 + 2 * i] + I * y[7 + 2 * i]; }
  force(z, F);
  double hub2 = (BG == 2) ? CH * rphi / 3.0 : 0.0;
  F[0] -= hub2 * z[0];
  double K = 0;
  for (int i = 0; i < 3; i++) {
    dy[2 * i] = creal(zd[i]); dy[2 * i + 1] = cimag(zd[i]);
    double complex acc = -3 * H * zd[i] - F[i];
    dy[6 + 2 * i] = creal(acc); dy[7 + 2 * i] = cimag(acc);
    K += creal(zd[i] * conj(zd[i]));
  }
  double V = potential(z) - hub2 * creal(z[0] * conj(z[0]));
  if (BG == 0) { dy[12] = H; dy[13] = 0; dy[14] = 0; }
  else { double A = exp(y[12]); dy[12] = H; dy[13] = -GPHI * y[13]; dy[14] = GPHI * y[13] * A; }
  /* explicit time dependence of V in bg 2: d(-cH rho_phi/3 |X|^2)/dt with d rho_phi/dt = -(3H+G) rho_phi */
  double dVdt = (BG == 2) ? CH * (3 * H + GPHI) * rphi / 3.0 * creal(z[0] * conj(z[0])) : 0.0;
  dy[15] = 3 * H * a * a * a * (V - K) + a * a * a * dVdt;
}

/* out: 0 t,1 a,2 H,3 NX,4 NY,5 NS,6 EX,7 Ncirc,8 eps,9 NYH,10 NYL,11 JYcodex,12 Nsig,13 gauss,14 ledger,
 *      15 a4rho_total,16 |B|/A (2mu|S|/M^2),17 EY_exc a^4 (sum modes),18 Esig a^4,19 rho_phi/(rho_phi+rho_r),
 *      20 a^4 rho_Y (full Y kinetic+quadratic+quartic), 21 Neff=NX+3NY */
#define NO 22
static void diagnostics(double t, const double *y, double *o) {
  double a, H, rphi; background(t, y, &a, &H, &rphi);
  double complex z[3], zd[3];
  double n[3], K = 0;
  for (int i = 0; i < 3; i++) {
    z[i] = y[2 * i] + I * y[2 * i + 1]; zd[i] = y[6 + 2 * i] + I * y[7 + 2 * i];
    n[i] = -2 * a * a * a * cimag(conj(z[i]) * zd[i]);
    K += creal(zd[i] * conj(zd[i]));
  }
  double a3 = a * a * a, a4 = a3 * a;
  double ax2 = creal(z[0] * conj(z[0]));
  double EX = a4 * (creal(zd[0] * conj(zd[0])) + lx * ax2 * ax2);
  double Nc = 2 * sqrt(2 * lx) * pow(EX / (3 * lx), 0.75);
  /* forced Y solution */
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
  /* Sigma radial */
  double r = as, rd = creal(conj(S) * Sd) / as;
  double Es = rd * rd + ls * (r * r - hh * hh / 2) * (r * r - hh * hh / 2);
  double msig = sqrt(2 * ls) * hh;
  double ay2 = creal(z[1] * conj(z[1]));
  double rhoY = creal(zd[1] * conj(zd[1])) + M2 * ay2 + ly * ay2 * ay2;
  double hub2 = (BG == 2) ? CH * rphi / 3.0 : 0.0;
  double rho = K + potential(z) - hub2 * ax2;
  double scale = fabs(n[0]) + 3 * fabs(n[1]) + 6 * fabs(n[2]);
  o[0] = t; o[1] = a; o[2] = H; o[3] = n[0]; o[4] = n[1]; o[5] = n[2]; o[6] = EX; o[7] = Nc; o[8] = n[0] / Nc;
  o[9] = a3 * E1 / w1; o[10] = a3 * E2 / w2; o[11] = a3 * Ecod / MM; o[12] = a3 * Es / msig;
  o[13] = (n[0] + 3 * n[1] + 6 * n[2]) / (scale > 1e-300 ? scale : 1.0);
  o[14] = a3 * rho - y[15];
  o[15] = a4 * rho; o[16] = 2 * mu * as / M2; o[17] = a4 * (E1 + E2); o[18] = a4 * Es;
  o[19] = (BG == 0) ? 0.0 : (y[13] / a3) / (y[13] / a3 + y[14] / a4);
  o[20] = a4 * rhoY; o[21] = n[0] + 3 * n[1];
}

int main(int argc, char **argv) {
  if (argc != 17) { fprintf(stderr, "usage: yrelease bg Gphi cH M h mu t_end rtol atol nsamp x0r x0i y0r y0i s0r s0i\n"); return 1; }
  BG = atoi(argv[1]); GPHI = atof(argv[2]); CH = atof(argv[3]);
  MM = atof(argv[4]); hh = atof(argv[5]); mu = atof(argv[6]);
  double t_end = atof(argv[7]), rtol = atof(argv[8]), atol = atof(argv[9]);
  int nsamp = atoi(argv[10]);
  double y[N] = {0};
  for (int i = 0; i < 6; i++) { y[i] = atof(argv[11 + i]); y[6 + i] = 0.0; }
  double t = (BG == 0) ? 0.5 : 2.0 / 3.0;
  y[12] = 0; y[13] = 3.0; y[14] = 0.0;
  { double complex z[3]; for (int i = 0; i < 3; i++) z[i] = y[2 * i] + I * y[2 * i + 1];
    double ax2 = creal(z[0] * conj(z[0]));
    y[15] = potential(z) - ((BG == 2) ? CH * ax2 : 0.0); }
  double W0 = fabs(y[15]) > 0 ? fabs(y[15]) : 1.0;
  /* windows by RD-equivalent a_eq = H^(-1/2) */
  const double edges[] = {100, 300, 600, 1000, 2000, 3000, 6000, 10000};
  const int NW = 7;
  double wsum[7][NO] = {{0}}, wspan[7] = {0};
  static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
    a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
    a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
    b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
    e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
  double k1[N], k2[N], k3[N], k4[N], k5[N], k6[N], k7[N], yt[N], yn[N];
  double h = 1e-3; long steps = 0, rejects = 0;
  double t0 = t, lr = log(t_end / t0); int si = 0;
  double next_s = t0 * exp(lr * (si + 1) / nsamp);
  double o[NO], oo[NO]; double maxg = 0, maxl = 0;
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
    for (int i = 0; i < 15; i++) {
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
      double aeq = pow(o[2], -0.5);
      for (int w = 0; w < NW; w++) if (aeq > edges[w] && aeq <= edges[w + 1]) {
        double da = o[1] - oo[1]; /* weight by da (conformal time in RD), as in Codex's window means */
        wspan[w] += da; for (int j = 0; j < NO; j++) wsum[w][j] += 0.5 * da * (o[j] + oo[j]); }
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
    if (!isfinite(y[0]) || h < 1e-14) { printf("FAIL %.17g\n", t); return 2; }
  }
  for (int w = 0; w < NW; w++) if (wspan[w] > 0) {
    printf("W %g %g %.15g", edges[w], edges[w + 1], wspan[w]);
    for (int j = 0; j < NO; j++) printf(" %.15g", wsum[w][j] / wspan[w]);
    printf("\n");
  }
  printf("F"); for (int i = 0; i < N; i++) printf(" %.17g", y[i]); printf("\n");
  printf("END %ld %ld %.6g %.6g\n", steps, rejects, maxg, maxl);
  return 0;
}
