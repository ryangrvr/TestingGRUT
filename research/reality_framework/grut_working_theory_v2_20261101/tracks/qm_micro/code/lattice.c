/* qm_micro item 1b (Claude, 2 Nov 2026): classical lattice evolution of the full M1 X-phi field theory on an
 * Nx x Ny grid (x: zero-gradient ends, y: periodic), rescaled units m_X = 1:
 *     psi_tt = lap psi - (S^2 + a |psi|^2) psi,
 *     S_tt   = lap S   - [b (S^2 - 1) S - 3 c (S - 1)^2] - |psi|^2 S,      a = 0.1, b = 3, c = 0.2.
 * Second-order Laplacian (dx, dy), velocity-Verlet in time (dt).  Own code.
 * Usage: lattice in.bin out.bin Nx Ny dx dy dt nsteps
 * File layout (doubles): pr, pi, S, vr, vi, vS, each Nx*Ny, index i*Ny + j (x-major).
 * Prints: energy and charge at start and end (energy E = sum [1/2|v|^2 + 1/2|grad|^2 + U] dx dy; charge
 * Q = sum (pi vr - pr vi) dx dy, the sign that gives psi = F e^{-i Om t} positive charge).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static const double A = 0.1, B = 3.0, C = 0.2;
static int Nx, Ny;
static double dx, dy;

static inline int id(int i, int j) { return i * Ny + j; }

static void accel(const double *pr, const double *pi, const double *S, double *ar, double *ai, double *aS) {
    double ix2 = 1.0 / (dx * dx), iy2 = 1.0 / (dy * dy);
    for (int i = 0; i < Nx; i++) {
        int im = i > 0 ? i - 1 : 1;          /* zero-gradient (mirror) ends */
        int ip = i < Nx - 1 ? i + 1 : Nx - 2;
        for (int j = 0; j < Ny; j++) {
            int jm = (j + Ny - 1) % Ny, jp = (j + 1) % Ny;
            int k = id(i, j);
            double lr = (pr[id(ip, j)] + pr[id(im, j)] - 2 * pr[k]) * ix2;
            double li = (pi[id(ip, j)] + pi[id(im, j)] - 2 * pi[k]) * ix2;
            double ls = (S[id(ip, j)] + S[id(im, j)] - 2 * S[k]) * ix2;
            if (Ny > 1) {
                lr += (pr[id(i, jp)] + pr[id(i, jm)] - 2 * pr[k]) * iy2;
                li += (pi[id(i, jp)] + pi[id(i, jm)] - 2 * pi[k]) * iy2;
                ls += (S[id(i, jp)] + S[id(i, jm)] - 2 * S[k]) * iy2;
            }
            double s = S[k], r2 = pr[k] * pr[k] + pi[k] * pi[k];
            double m = s * s + A * r2;
            ar[k] = lr - m * pr[k];
            ai[k] = li - m * pi[k];
            aS[k] = ls - (B * (s * s - 1) * s - 3 * C * (s - 1) * (s - 1)) - r2 * s;
        }
    }
}

static void measures(const double *pr, const double *pi, const double *S, const double *vr, const double *vi,
                     const double *vS, double *E, double *Q) {
    double e = 0, q = 0;
    for (int i = 0; i < Nx; i++) {
        int ip = i < Nx - 1 ? i + 1 : i;     /* forward differences, last bond dropped (zero gradient) */
        for (int j = 0; j < Ny; j++) {
            int jp = (j + 1) % Ny, k = id(i, j);
            double gx = 0, gy = 0;
            gx = ((pr[id(ip, j)] - pr[k]) * (pr[id(ip, j)] - pr[k]) + (pi[id(ip, j)] - pi[k]) * (pi[id(ip, j)] - pi[k]) +
                  (S[id(ip, j)] - S[k]) * (S[id(ip, j)] - S[k])) / (dx * dx);
            if (Ny > 1)
                gy = ((pr[id(i, jp)] - pr[k]) * (pr[id(i, jp)] - pr[k]) + (pi[id(i, jp)] - pi[k]) * (pi[id(i, jp)] - pi[k]) +
                      (S[id(i, jp)] - S[k]) * (S[id(i, jp)] - S[k])) / (dy * dy);
            double s = S[k], r2 = pr[k] * pr[k] + pi[k] * pi[k];
            double U = B / 4 * (s * s - 1) * (s * s - 1) - C * (s - 1) * (s - 1) * (s - 1) + 0.5 * s * s * r2 + A / 4 * r2 * r2;
            double wt = (i == 0 || i == Nx - 1) ? 0.5 : 1.0;   /* trapezoid end weights (mirror ends) */
            e += wt * (0.5 * (vr[k] * vr[k] + vi[k] * vi[k] + vS[k] * vS[k]) + U) + 0.5 * gx + 0.5 * wt * gy;
            q += wt * (pi[k] * vr[k] - pr[k] * vi[k]);
        }
    }
    *E = e * dx * dy; *Q = q * dx * dy;
}

int main(int argc, char **argv) {
    if (argc < 9) { fprintf(stderr, "usage: lattice in out Nx Ny dx dy dt nsteps\n"); return 1; }
    Nx = atoi(argv[3]); Ny = atoi(argv[4]); dx = atof(argv[5]); dy = atof(argv[6]);
    double dt = atof(argv[7]); long ns = atol(argv[8]);
    size_t n = (size_t)Nx * Ny;
    double *buf = malloc(6 * n * sizeof(double));
    double *acc = malloc(3 * n * sizeof(double));
    FILE *f = fopen(argv[1], "rb");
    if (!f || fread(buf, sizeof(double), 6 * n, f) != 6 * n) { fprintf(stderr, "read error\n"); return 2; }
    fclose(f);
    double *pr = buf, *pi = buf + n, *S = buf + 2 * n, *vr = buf + 3 * n, *vi = buf + 4 * n, *vS = buf + 5 * n;
    double *ar = acc, *ai = acc + n, *aS = acc + 2 * n;
    double E0, Q0, E1, Q1;
    measures(pr, pi, S, vr, vi, vS, &E0, &Q0);
    accel(pr, pi, S, ar, ai, aS);
    for (long t = 0; t < ns; t++) {
        for (size_t k = 0; k < n; k++) {
            vr[k] += 0.5 * dt * ar[k]; vi[k] += 0.5 * dt * ai[k]; vS[k] += 0.5 * dt * aS[k];
            pr[k] += dt * vr[k]; pi[k] += dt * vi[k]; S[k] += dt * vS[k];
        }
        accel(pr, pi, S, ar, ai, aS);
        for (size_t k = 0; k < n; k++) {
            vr[k] += 0.5 * dt * ar[k]; vi[k] += 0.5 * dt * ai[k]; vS[k] += 0.5 * dt * aS[k];
        }
    }
    measures(pr, pi, S, vr, vi, vS, &E1, &Q1);
    printf("E0 %.15e Q0 %.15e E1 %.15e Q1 %.15e\n", E0, Q0, E1, Q1);
    f = fopen(argv[2], "wb");
    fwrite(buf, sizeof(double), 6 * n, f);
    fclose(f);
    free(buf); free(acc);
    return 0;
}
