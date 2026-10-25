/* kg_radial.c -- independent radial solver for the spectrum track (Claude, 25 Oct 2026).
 *
 * Fields: U = r psi (complex, components x,y), Z = r (S - 1) (real).
 * Equations (dimensionless, as in the Codex MATHEMATICS.md files, re-implemented here):
 *   U_tt = U_rr - (S^2 + 0.1 |psi|^2) U,
 *   Z_tt = Z_rr - r S [3 (S^2 - 1) + |psi|^2],      S = 1 + Z/r, |psi|^2 = |U|^2 / r^2.
 * Linear mode (nonlinear=0): U_tt = U_rr - U,  Z_tt = Z_rr - 6 Z.
 *
 * Method (chosen to differ from Codex ff49e5d and from the independent_solver track):
 *   space: 6th-order centred finite differences, odd continuation at r=0 and at the hard wall r=L;
 *   time : 4th-order Yoshida composition of velocity-Verlet (symplectic; the U(1) Noether
 *          charge 4 pi h sum (y Vx - x Vy) is conserved exactly by every drift and every kick);
 *   boundary: a hard wall far enough away that nothing reflected can reach the extraction radii
 *          inside the run window (checked in the driver), so no absorber is used at all;
 *   cost  : only nodes inside r < r_act0 + t + margin are advanced (field is ~0 outside).
 *
 * Usage: kg_radial INIT OUTPREFIX dt tmax nonlinear sample_dt Rext(csv) snaps(csv) hist_dt r_act0 Rin(csv)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define G 3
static long N;           /* nodes 0..N, wall at N */
static double h;
static int nonlinear;
static double *qx, *qy, *qz, *px, *py, *pz, *ax, *ay, *az;
static long iact;        /* last active interior node */

static void fill_ghosts(double *q) {
  for (int j = 1; j <= G; j++) q[G - j] = -q[G + j];
  if (iact >= N - G - 1)
    for (int j = 1; j <= G; j++) q[G + N + j] = -q[G + N - j];
}

static inline double lap6(const double *q, long i) {
  const double *c = q + G + i;
  return ((1.0 / 90) * (c[-3] + c[3]) - (3.0 / 20) * (c[-2] + c[2]) + 1.5 * (c[-1] + c[1]) - (49.0 / 18) * c[0]) / (h * h);
}
static inline double d1_6(const double *q, long i) {
  const double *c = q + G + i;
  return ((3.0 / 4) * (c[1] - c[-1]) - (3.0 / 20) * (c[2] - c[-2]) + (1.0 / 60) * (c[3] - c[-3])) / h;
}

static void accel(void) {
  fill_ghosts(qx); fill_ghosts(qy); fill_ghosts(qz);
  long top = iact < N - 1 ? iact : N - 1;
  for (long i = 1; i <= top; i++) {
    double r = i * h, x = qx[G + i], y = qy[G + i], z = qz[G + i];
    double fx, fy, fz;
    if (nonlinear) {
      double dens = (x * x + y * y) / (r * r), S = 1 + z / r, m2 = S * S + 0.1 * dens;
      fx = m2 * x; fy = m2 * y; fz = r * S * (3 * (S * S - 1) + dens);
    } else { fx = x; fy = y; fz = 6 * z; }
    ax[G + i] = lap6(qx, i) - fx;
    ay[G + i] = lap6(qy, i) - fy;
    az[G + i] = lap6(qz, i) - fz;
  }
}
static void drift(double c) {
  long top = iact < N - 1 ? iact : N - 1;
  for (long i = 1; i <= top; i++) { qx[G + i] += c * px[G + i]; qy[G + i] += c * py[G + i]; qz[G + i] += c * pz[G + i]; }
}
static void kick(double c) {
  accel();
  long top = iact < N - 1 ? iact : N - 1;
  for (long i = 1; i <= top; i++) { px[G + i] += c * ax[G + i]; py[G + i] += c * ay[G + i]; pz[G + i] += c * az[G + i]; }
}
static void step(double dt) {
  const double w1 = 1.0 / (2.0 - cbrt(2.0)), w0 = -cbrt(2.0) * w1;
  const double c1 = w1 / 2, c2 = (w0 + w1) / 2;
  drift(c1 * dt); kick(w1 * dt); drift(c2 * dt); kick(w0 * dt); drift(c2 * dt); kick(w1 * dt); drift(c1 * dt);
}

static double node_energy(long i) {
  double r = i * h, x = qx[G + i], y = qy[G + i], z = qz[G + i];
  double e = 0.5 * (px[G + i] * px[G + i] + py[G + i] * py[G + i] + pz[G + i] * pz[G + i]);
  e -= 0.5 * (x * lap6(qx, i) + y * lap6(qy, i) + z * lap6(qz, i));
  if (nonlinear) {
    double dens = (x * x + y * y) / (r * r), S = 1 + z / r, s2 = S * S - 1;
    e += r * r * (0.75 * s2 * s2 + 0.5 * S * S * dens + 0.025 * dens * dens);
  } else e += 0.5 * (x * x + y * y) + 3 * z * z;
  return e;
}

static int parse_csv(const char *s, double *out, int maxn) {
  int n = 0; char buf[4096]; strncpy(buf, s, sizeof buf - 1); buf[sizeof buf - 1] = 0;
  if (!strcmp(buf, "none")) return 0;
  for (char *t = strtok(buf, ","); t && n < maxn; t = strtok(NULL, ",")) out[n++] = atof(t);
  return n;
}

int main(int argc, char **argv) {
  if (argc != 12) { fprintf(stderr, "usage: INIT OUT dt tmax nonlinear sample_dt Rext snaps hist_dt r_act0 Rin\n"); return 1; }
  FILE *fi = fopen(argv[1], "rb"); if (!fi) { perror("init"); return 1; }
  long long NN; double hh;
  if (fread(&NN, sizeof NN, 1, fi) != 1 || fread(&hh, sizeof hh, 1, fi) != 1) return 1;
  N = NN; h = hh;
  double **arr[6] = {&qx, &qy, &qz, &px, &py, &pz};
  for (int k = 0; k < 6; k++) {
    *arr[k] = calloc(N + 1 + 2 * G, sizeof(double));
    if (fread(*arr[k] + G, sizeof(double), N + 1, fi) != (size_t)(N + 1)) { fprintf(stderr, "short init\n"); return 1; }
    (*arr[k])[G] = 0; (*arr[k])[G + N] = 0;
  }
  fclose(fi);
  ax = calloc(N + 1 + 2 * G, sizeof(double)); ay = calloc(N + 1 + 2 * G, sizeof(double)); az = calloc(N + 1 + 2 * G, sizeof(double));
  double dt = atof(argv[3]), tmax = atof(argv[4]); nonlinear = atoi(argv[5]);
  double sdt = atof(argv[6]); double Rext[64], snaps[64], Rin[64];
  int ne = parse_csv(argv[7], Rext, 64), ns = parse_csv(argv[8], snaps, 64);
  double hdt = atof(argv[9]), ract0 = atof(argv[10]); int nin = parse_csv(argv[11], Rin, 64);
  long iext[64]; for (int k = 0; k < ne; k++) iext[k] = lround(Rext[k] / h);
  long iin[64]; for (int k = 0; k < nin; k++) iin[k] = lround(Rin[k] / h);
  long nsteps = lround(tmax / dt), severy = lround(sdt / dt), hevery = lround(hdt / dt);
  char fn[4096];
  snprintf(fn, sizeof fn, "%s_ext.bin", argv[2]); FILE *fe = fopen(fn, "wb");
  snprintf(fn, sizeof fn, "%s_hist.bin", argv[2]); FILE *fh = fopen(fn, "wb");
  const double PI = 3.14159265358979323846, W = 4 * PI * h;
  int snapdone[64] = {0};
  for (long n = 0; n <= nsteps; n++) {
    double t = n * dt;
    iact = lround((ract0 + t + 400.0) / h); if (iact > N - 1) iact = N - 1;
    fill_ghosts(qx); fill_ghosts(qy); fill_ghosts(qz);
    if (n % severy == 0) {
      fwrite(&t, sizeof t, 1, fe);
      for (int k = 0; k < ne; k++) {
        long i = iext[k];
        double rec[9] = {qx[G + i], qy[G + i], qz[G + i], px[G + i], py[G + i], pz[G + i], d1_6(qx, i), d1_6(qy, i), d1_6(qz, i)};
        fwrite(rec, sizeof(double), 9, fe);
      }
    }
    if (n % hevery == 0 || n == nsteps) {
      double E = 0, Q = 0, Qin[64] = {0}, Ein[64] = {0}, maxdm2_far = 0;
      long top = iact < N - 1 ? iact : N - 1;
      for (long i = 1; i <= top; i++) {
        double e = node_energy(i), q = qy[G + i] * px[G + i] - qx[G + i] * py[G + i];
        E += W * e; Q += W * q;
        for (int k = 0; k < nin; k++) {
          if (i < iin[k]) { Qin[k] += W * q; Ein[k] += W * e; }
          else if (i == iin[k]) { Qin[k] += 0.5 * W * q; Ein[k] += 0.5 * W * e; }
        }
      }
      fwrite(&t, sizeof t, 1, fh); fwrite(&E, sizeof E, 1, fh); fwrite(&Q, sizeof Q, 1, fh);
      for (int k = 0; k < nin; k++) { fwrite(&Qin[k], sizeof(double), 1, fh); fwrite(&Ein[k], sizeof(double), 1, fh); }
      (void)maxdm2_far;
      if (n % (hevery * 1000) == 0) { fprintf(stderr, "t=%.1f E=%.12g Q=%.15g\n", t, E, Q); fflush(fe); fflush(fh); }
    }
    for (int k = 0; k < ns; k++)
      if (!snapdone[k] && fabs(t - snaps[k]) < 0.5 * dt) {
        snprintf(fn, sizeof fn, "%s_snap%.0f.bin", argv[2], snaps[k]); FILE *fs = fopen(fn, "wb");
        long long nn = N; fwrite(&nn, sizeof nn, 1, fs); fwrite(&h, sizeof h, 1, fs); fwrite(&t, sizeof t, 1, fs);
        for (int a = 0; a < 6; a++) fwrite(*arr[a] + G, sizeof(double), N + 1, fs);
        fclose(fs); snapdone[k] = 1;
      }
    if (n < nsteps) step(dt);
  }
  fclose(fe); fclose(fh);
  return 0;
}
