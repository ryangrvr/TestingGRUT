/* gp_phase track (Claude, own code): homogeneous X, Y, Sigma release at the 7737487 gauge point,
 * with one optional extra gauge-invariant complex operator.
 *
 * Physical potential (H_inf = 1 units, kappa and mu real positive by rephasing):
 *   V = lx|X|^4 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + M^2|Y|^2
 *       + 2 Re(kappa Y* X^3) + 2 Re(mu S* Y^2) + 2 Re(C O_op)
 *   op 0: none
 *   op 1: O = S* X^6          (dimension 7, rephasing class kappa^2 mu)
 *   op 2: O = S* Y X^3        (dimension 5, class kappa mu)
 *   op 3: O = S* Y^2 |X|^2    (dimension 5, class mu)
 *   op 4: O = Y* X^3 |X|^2    (dimension 6, class kappa)
 * EOM z'' + 3H z' + dV/dz* = 0. Radiation era after release, a = tau (conformal), a(1) = 1, H(1) = 1.
 * Conformal chi = a z obeys chi'' = -a^3 F(chi/a), F = dV/dz*.
 * Extra ODE variable G integrates d(a^4 rho)/dtau = (a^4/tau)(4 rho - 6 K) for an energy-ledger control.
 * Integrator: Dormand-Prince 5(4), adaptive, own implementation.
 * Usage: release op Cre Cim lx ly ls kappa mu M h tau_end rtol atol nsamp  x0r x0i y0r y0i s0r s0i  (initial z, zero velocity)
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

static int OP; static double complex C;
static double lx, ly, ls, kap, mu, M, hh;

static void force(const double complex *z, double complex *F) {
  double complex X = z[0], Y = z[1], S = z[2];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  double complex X2 = X * X, X3 = X2 * X, Xc = conj(X), Yc = conj(Y), Sc = conj(S);
  F[0] = 2 * lx * ax2 * X + 3 * kap * Y * Xc * Xc;
  F[1] = M * M * Y + 2 * ly * ay2 * Y + kap * X3 + 2 * mu * S * Yc;
  F[2] = 2 * ls * (as2 - hh * hh / 2) * S + mu * Y * Y;
  double complex Cc = conj(C);
  switch (OP) {
    case 1: { double complex Xc5 = Xc * Xc * Xc * Xc * Xc;
      F[0] += 6 * Cc * S * Xc5; F[2] += C * X3 * X3; break; }
    case 2:
      F[0] += 3 * Cc * S * Yc * Xc * Xc; F[1] += Cc * S * Xc * Xc * Xc; F[2] += C * Y * X3; break;
    case 3:
      F[0] += 2 * creal(C * Sc * Y * Y) * X; F[1] += 2 * Cc * S * Yc * ax2; F[2] += C * Y * Y * ax2; break;
    case 4:
      F[0] += C * Yc * X3 * X + 4 * Cc * Y * X * Xc * Xc * Xc; F[1] += C * X3 * ax2; break;
    default: break;
  }
}

static double potential(const double complex *z) {
  double complex X = z[0], Y = z[1], S = z[2];
  double ax2 = creal(X * conj(X)), ay2 = creal(Y * conj(Y)), as2 = creal(S * conj(S));
  double complex X3 = X * X * X;
  double v = lx * ax2 * ax2 + ly * ay2 * ay2 + ls * (as2 - hh * hh / 2) * (as2 - hh * hh / 2) + M * M * ay2
           + 2 * creal(kap * conj(Y) * X3) + 2 * creal(mu * conj(S) * Y * Y);
  switch (OP) {
    case 1: v += 2 * creal(C * conj(S) * X3 * X3); break;
    case 2: v += 2 * creal(C * conj(S) * Y * X3); break;
    case 3: v += 2 * creal(C * conj(S) * Y * Y) * ax2; break;
    case 4: v += 2 * creal(C * conj(Y) * X3) * ax2; break;
    default: break;
  }
  return v;
}

/* state: y[0..5] = chi (re,im) x3, y[6..11] = chi' (re,im) x3, y[12] = G */
#define N 13
static void rhs(double tau, const double *y, double *dy) {
  double a = tau;
  double complex chi[3], pi[3], z[3], zd[3], F[3];
  for (int i = 0; i < 3; i++) {
    chi[i] = y[2 * i] + I * y[2 * i + 1];
    pi[i] = y[6 + 2 * i] + I * y[7 + 2 * i];
    z[i] = chi[i] / a;
    zd[i] = (pi[i] - chi[i] / tau) / (a * a); /* dz/dt with a'=1 */
  }
  force(z, F);
  double K = 0;
  for (int i = 0; i < 3; i++) {
    dy[2 * i] = y[6 + 2 * i]; dy[2 * i + 1] = y[7 + 2 * i];
    double complex acc = -a * a * a * F[i];
    dy[6 + 2 * i] = creal(acc); dy[7 + 2 * i] = cimag(acc);
    K += creal(zd[i] * conj(zd[i]));
  }
  double rho = K + potential(z);
  dy[12] = (a * a * a * a / tau) * (4 * rho - 6 * K);
}

static void diagnostics(double tau, const double *y, double *out) {
  /* out: a3n_X, a3n_Y, a3n_S, a4rho_X, eps_X, gauss_rel, a4rho_total, |X|a, |Y|a, |S|a */
  double a = tau; double complex chi[3], pi[3], z[3], zd[3];
  double n[3], K = 0;
  for (int i = 0; i < 3; i++) {
    chi[i] = y[2 * i] + I * y[2 * i + 1];
    pi[i] = y[6 + 2 * i] + I * y[7 + 2 * i];
    z[i] = chi[i] / a; zd[i] = (pi[i] - chi[i] / tau) / (a * a);
    n[i] = -2 * cimag(conj(chi[i]) * pi[i]);
    K += creal(zd[i] * conj(zd[i]));
  }
  double complex vX = pi[0] - chi[0] / tau;
  double ax2 = creal(chi[0] * conj(chi[0]));
  double a4rhoX = creal(vX * conj(vX)) + lx * ax2 * ax2;
  double nmax = 2 * sqrt(2 * lx) * pow(a4rhoX / (3 * lx), 0.75);
  double scale = fabs(n[0]) + 3 * fabs(n[1]) + 6 * fabs(n[2]);
  out[0] = n[0]; out[1] = n[1]; out[2] = n[2]; out[3] = a4rhoX; out[4] = n[0] / nmax;
  out[5] = (n[0] + 3 * n[1] + 6 * n[2]) / (scale > 1e-300 ? scale : 1.0);
  out[6] = a * a * a * a * (K + potential(z));
  out[7] = cabs(chi[0]); out[8] = cabs(chi[1]); out[9] = cabs(chi[2]);
}

int main(int argc, char **argv) {
  if (argc < 21) { fprintf(stderr, "bad args\n"); return 1; }
  OP = atoi(argv[1]); C = atof(argv[2]) + I * atof(argv[3]);
  lx = atof(argv[4]); ly = atof(argv[5]); ls = atof(argv[6]); kap = atof(argv[7]); mu = atof(argv[8]);
  M = atof(argv[9]); hh = atof(argv[10]);
  double tau_end = atof(argv[11]), rtol = atof(argv[12]), atol = atof(argv[13]);
  int nsamp = atoi(argv[14]);
  double y[N];
  for (int i = 0; i < 6; i++) y[i] = atof(argv[15 + i]);
  /* zero physical velocity at a=1, a'=1: chi' = z */
  for (int i = 0; i < 6; i++) y[6 + i] = y[i];
  { double a = 1.0; double complex z[3]; double K = 0;
    for (int i = 0; i < 3; i++) z[i] = y[2 * i] + I * y[2 * i + 1];
    y[12] = a * a * a * a * (K + potential(z)); }
  double E0 = y[12];
  double X0abs = hypot(y[0], y[1]), Y0abs = hypot(y[2], y[3]), S0abs = hypot(y[4], y[5]);
  /* Dormand-Prince 5(4) tableau */
  static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
    a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
    a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
    b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
    e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
  double k1[N], k2[N], k3[N], k4[N], k5[N], k6[N], k7[N], yt[N], yn[N];
  double tau = 1.0, h = 1e-4; long steps = 0, rejects = 0;
  double next_s; int si = 0;
  double lr = log(tau_end);
  next_s = exp(lr * (si + 1) / nsamp);
  double out[10]; double max_gauss = 0, max_ledger = 0;
  diagnostics(tau, y, out);
  printf("S %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", tau, out[0], out[1], out[2], out[3], out[4], out[5], 0.0, out[6], out[7], out[8], out[9]);
  rhs(tau, y, k1);
  while (tau < tau_end) {
    int land = 0;
    if (tau + h >= next_s) { h = next_s - tau; land = 1; }
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * a21 * k1[i];
    rhs(tau + c2 * h, yt, k2);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
    rhs(tau + c3 * h, yt, k3);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    rhs(tau + c4 * h, yt, k4);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    rhs(tau + c5 * h, yt, k5);
    for (int i = 0; i < N; i++) yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    rhs(tau + h, yt, k6);
    for (int i = 0; i < N; i++) yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    rhs(tau + h, yn, k7);
    double err = 0;
    for (int i = 0; i < 12; i++) {
      double ei = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
      double sc = atol + rtol * fmax(fabs(y[i]), fabs(yn[i]));
      double r = fabs(ei) / sc; if (r > err) err = r;
    }
    if (err <= 1.0) {
      tau += h; steps++;
      for (int i = 0; i < N; i++) { y[i] = yn[i]; k1[i] = k7[i]; }
      diagnostics(tau, y, out);
      double ledger = (out[6] - y[12]) / E0;
      if (fabs(out[5]) > max_gauss) max_gauss = fabs(out[5]);
      if (fabs(ledger) > max_ledger) max_ledger = fabs(ledger);
      if (land) {
        printf("S %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", tau, out[0], out[1], out[2], out[3], out[4], out[5], ledger, out[6], out[7], out[8], out[9]);
        si++; next_s = (si + 1 <= nsamp) ? exp(lr * (si + 1) / nsamp) : tau_end * 2;
        if (si >= nsamp) break;
      }
      double fac = (err > 0) ? 0.9 * pow(err, -0.2) : 5.0;
      if (fac > 5) fac = 5;
      if (fac < 0.2) fac = 0.2;
      h *= fac;
      if (land && h > 0) { /* keep h */ }
    } else {
      rejects++;
      double fac = 0.9 * pow(err, -0.2); if (fac < 0.1) fac = 0.1;
      h *= fac;
    }
    if (!isfinite(y[0]) || h < 1e-14) { printf("FAIL %.17g\n", tau); return 2; }
    { double ax = hypot(y[0], y[1]) / tau, ay = hypot(y[2], y[3]) / tau, as = hypot(y[4], y[5]) / tau;
      if (ax > 3 * X0abs || ay > 30 * Y0abs || as > 30 * S0abs) { printf("RUNAWAY %.17g %.6g %.6g %.6g\n", tau, ax, ay, as); return 3; } }
  }
  printf("END steps %ld rejects %ld max_gauss %.6g max_ledger %.6g\n", steps, rejects, max_gauss, max_ledger);
  return 0;
}
