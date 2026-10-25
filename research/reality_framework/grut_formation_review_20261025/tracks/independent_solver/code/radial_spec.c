/* Fast C port of spectral_radial.py (same algorithm; this track's own code).
 * Fourier-sine pseudo-spectral radial solver, exact free Klein-Gordon rotation per sine mode,
 * nonlinear position kick, Strang (scheme 2) or Yoshida triple-jump (scheme 4) composition,
 * optional smooth velocity sponge applied exactly inside the kick.  Fields: rows 0,1 = U = r psi
 * (complex, mass 1), row 2 = Z = r(S-1) (mass sqrt 6).  a = 0.1, b = 3.
 *
 * Build:  gcc -O3 -march=native -o radial_spec radial_spec.c /lib/x86_64-linux-gnu/libfftw3.so.3 -lm
 * (No FFTW header is installed; the few prototypes used are declared below.)
 *
 * Usage: radial_spec init.bin out_prefix h L dt tmax scheme diag Rs D gmax power snaptimes
 *   init.bin: 6*N doubles, physical q rows (x,y,Z) then p rows (px,py,pz) at nodes r_j = j h.
 *   Rs < 0 disables the sponge.  snaptimes: comma list or "-".
 * Output: out_prefix_history.bin (double rows, NCOL columns), out_prefix_snapT.bin (6N doubles).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct fftw_plan_s *fftw_plan;
typedef int fftw_r2r_kind;
extern fftw_plan fftw_plan_many_r2r(int rank, const int *n, int howmany, double *in, const int *inembed,
                                    int istride, int idist, double *out, const int *onembed, int ostride,
                                    int odist, const fftw_r2r_kind *kind, unsigned flags);
extern void fftw_execute_r2r(const fftw_plan p, double *in, double *out);
#define FFTW_REDFT00 3
#define FFTW_RODFT00 7
#define FFTW_MEASURE 0u
#define FFTW_ESTIMATE (1u << 6)

static int N; static double h, L, *r, *kk, *om, *gam = NULL;
static double *Qh, *Ph, *q, *p, *F, *tmp, *dcoef, *dout;
static fftw_plan pdst, pdct;
static double dstnorm;
static double *amalloc(size_t n) { void *ptr = NULL; if (posix_memalign(&ptr, 64, n * sizeof(double))) { fprintf(stderr, "alloc\n"); exit(3); } return (double *)ptr; }
static double Erec = 0, Qrec = 0;
static const double A = 0.1;

static void dst3(double *in, double *out) { /* orthonormal DST-I of 3 rows */
    fftw_execute_r2r(pdst, in, out);
    for (int i = 0; i < 3 * N; i++) out[i] *= dstnorm;
}

typedef struct { double s; double *c, *so, *wo; } Rot;
static Rot rots[8]; static int nrot = 0;
static Rot *getrot(double s) {
    for (int i = 0; i < nrot; i++) if (fabs(rots[i].s - s) < 1e-15) return &rots[i];
    Rot *R = &rots[nrot++]; R->s = s;
    R->c = amalloc(3 * N); R->so = amalloc(3 * N); R->wo = amalloc(3 * N);
    for (int i = 0; i < 3 * N; i++) { double w = om[i], c = cos(w * s), sn = sin(w * s); R->c[i] = c; R->so[i] = sn / w; R->wo[i] = -sn * w; }
    return R;
}
static void drift(double s) {
    Rot *R = getrot(s);
    for (int i = 0; i < 3 * N; i++) { double a = Qh[i], b = Ph[i]; Qh[i] = R->c[i] * a + R->so[i] * b; Ph[i] = R->wo[i] * a + R->c[i] * b; }
}
static void force(const double *qq, double *FF) {
    const double *x = qq, *y = qq + N, *z = qq + 2 * N;
    for (int j = 0; j < N; j++) {
        double rr = r[j], u2 = x[j] * x[j] + y[j] * y[j], S = 1.0 + z[j] / rr;
        double g = (S * S - 1.0) + A * u2 / (rr * rr);
        FF[j] = -g * x[j]; FF[N + j] = -g * y[j];
        FF[2 * N + j] = -(9.0 * z[j] * z[j] / rr + 3.0 * z[j] * z[j] * z[j] / (rr * rr) + u2 * S / rr);
    }
}
static void kick(double s) {
    dst3(Qh, q); force(q, F);
    if (!gam) {
        dst3(F, tmp);
        for (int i = 0; i < 3 * N; i++) Ph[i] += s * tmp[i];
        return;
    }
    dst3(Ph, p);
    double c4 = 4 * M_PI * h, de = 0, dq = 0;
    for (int j = 0; j < N; j++) {
        double d = exp(-gam[j] * s);
        double p1[3], p2[3];
        for (int a = 0; a < 3; a++) { p1[a] = p[a * N + j] + 0.5 * s * F[a * N + j]; p2[a] = p1[a] * d; }
        if (gam[j] > 0) {
            for (int a = 0; a < 3; a++) de += 0.5 * (p1[a] * p1[a] - p2[a] * p2[a]);
            double x = q[j], y = q[N + j];
            dq += (y * p1[0] - x * p1[1]) - (y * p2[0] - x * p2[1]);
        }
        for (int a = 0; a < 3; a++) p[a * N + j] = p2[a] + 0.5 * s * F[a * N + j];
    }
    Erec += c4 * de; Qrec += c4 * dq;
    dst3(p, Ph);
}
/* derivative of 3 sine series at nodes j=1..N (returns into out[3N]) */
static void deriv(const double *Xh, double *out) {
    int M = N + 2;
    for (int a = 0; a < 3; a++) {
        dcoef[a * M] = 0; dcoef[a * M + M - 1] = 0;
        for (int j = 0; j < N; j++) dcoef[a * M + 1 + j] = Xh[a * N + j] * kk[j];
    }
    fftw_execute_r2r(pdct, dcoef, dout);
    double c = sqrt(2.0 / (N + 1)) * 0.5;
    for (int a = 0; a < 3; a++) for (int j = 0; j < N; j++) out[a * N + j] = c * dout[a * M + 1 + j];
}

#define NR 5
static const double RAD[NR] = {30, 60, 120, 240, 480};
#define NCOL (7 + 6 * NR)
static double *dq_, *dp_;
static void diagnostics(double tau, double *row) {
    dst3(Qh, q); dst3(Ph, p); deriv(Qh, dq_); deriv(Ph, dp_);
    double c4 = 4 * M_PI * h;
    double Ek = 0, Eg = 0, Ep = 0, Q = 0, El = 0;
    for (int i = 0; i < 3 * N; i++) { Ek += 0.5 * Ph[i] * Ph[i]; Eg += 0.5 * kk[i % N] * kk[i % N] * Qh[i] * Qh[i]; }
    double *x = q, *y = q + N, *z = q + 2 * N, *px = p, *py = p + N, *pz = p + 2 * N;
    double accQ[NR] = {0}, accE[NR] = {0}, Ech[NR] = {0}, Ene[NR] = {0}, Qtrap[NR] = {0}, Qcor[NR] = {0};
    int node[NR]; for (int k = 0; k < NR; k++) node[k] = (int)lround(RAD[k] / h) - 1;
    for (int j = 0; j < N; j++) {
        double rr = r[j], u2 = x[j] * x[j] + y[j] * y[j], S = 1 + z[j] / rr;
        double potU = 0.5 * S * S * u2 + 0.25 * A * u2 * u2 / (rr * rr);
        double potZ = 0.75 * z[j] * z[j] * (2 + z[j] / rr) * (2 + z[j] / rr);
        Ep += potU + potZ;
        double gU = pow(dq_[j] - x[j] / rr, 2) + pow(dq_[N + j] - y[j] / rr, 2), gZ = pow(dq_[2 * N + j] - z[j] / rr, 2);
        double ech = 0.5 * (px[j] * px[j] + py[j] * py[j]) + 0.5 * gU + 0.5 * u2;
        double eint = 0.5 * (S * S - 1) * u2 + 0.25 * A * u2 * u2 / (rr * rr);
        double ene = 0.5 * pz[j] * pz[j] + 0.5 * gZ + potZ;
        double qd = y[j] * px[j] - x[j] * py[j];
        Q += qd; El += ech + eint + ene;
        for (int k = 0; k < NR; k++) {
            double w = (j < node[k]) ? 1.0 : (j == node[k] ? 0.5 : 0.0);
            accQ[k] += w * qd; accE[k] += w * (ech + eint + ene);
            double wo = (j > node[k]) ? 1.0 : (j == node[k] ? 0.5 : 0.0);
            Ech[k] += wo * ech; Ene[k] += wo * ene;
            if (j == node[k]) {
                double dqd = dq_[N + j] * px[j] + y[j] * dp_[j] - dq_[j] * py[j] - x[j] * dp_[N + j];
                Qcor[k] = -4 * M_PI * h * h / 12.0 * dqd;
            }
        }
    }
    int c = 0;
    row[c++] = tau; row[c++] = c4 * (Ek + Eg + Ep); row[c++] = c4 * Q; row[c++] = Erec; row[c++] = Qrec;
    row[c++] = c4 * El; row[c++] = 1 + z[0] / r[0];
    for (int k = 0; k < NR; k++) {
        Qtrap[k] = c4 * accQ[k];
        row[c++] = Qtrap[k] + Qcor[k]; row[c++] = Qtrap[k]; row[c++] = c4 * accE[k];
        row[c++] = c4 * Ech[k]; row[c++] = c4 * Ene[k]; row[c++] = 0.0;
    }
}

int main(int argc, char **argv) {
    if (argc < 15) { fprintf(stderr, "usage\n"); return 1; }
    const char *init = argv[1], *pre = argv[2];
    h = atof(argv[3]); L = atof(argv[4]); double dt = atof(argv[5]), tmax = atof(argv[6]);
    int scheme = atoi(argv[7]); double diag = atof(argv[8]);
    double Rs = atof(argv[9]), D = atof(argv[10]), gmax = atof(argv[11]), pw = atof(argv[12]);
    const char *snaps = argv[13]; double ckpt = atof(argv[14]);
    N = (int)lround(L / h) - 1;
    r = amalloc(N); kk = amalloc(N); om = amalloc(3 * N);
    for (int j = 0; j < N; j++) { r[j] = h * (j + 1); kk[j] = M_PI * (j + 1) / L; }
    double m2[3] = {1, 1, 6};
    for (int a = 0; a < 3; a++) for (int j = 0; j < N; j++) om[a * N + j] = sqrt(kk[j] * kk[j] + m2[a]);
    if (Rs >= 0) { gam = amalloc(N); for (int j = 0; j < N; j++) { double xi = (r[j] - Rs) / D; xi = xi < 0 ? 0 : (xi > 1 ? 1 : xi); gam[j] = gmax * pow(xi, pw); } }
    Qh = amalloc(3 * N); Ph = amalloc(3 * N); q = amalloc(3 * N);
    p = amalloc(3 * N); F = amalloc(3 * N); tmp = amalloc(3 * N);
    dq_ = amalloc(3 * N); dp_ = amalloc(3 * N);
    dcoef = amalloc(3 * (N + 2)); dout = amalloc(3 * (N + 2));
    int n1[1] = {N}, n2[1] = {N + 2}; fftw_r2r_kind k1[1] = {FFTW_RODFT00}, k2[1] = {FFTW_REDFT00};
    pdst = fftw_plan_many_r2r(1, n1, 3, q, NULL, 1, N, p, NULL, 1, N, k1, FFTW_MEASURE);
    pdct = fftw_plan_many_r2r(1, n2, 3, dcoef, NULL, 1, N + 2, dout, NULL, 1, N + 2, k2, FFTW_ESTIMATE);
    dstnorm = 1.0 / sqrt(2.0 * (N + 1));
    FILE *fi = fopen(init, "rb"); if (!fi) { perror("init"); return 1; }
    if (fread(q, sizeof(double), 3 * N, fi) != (size_t)(3 * N) || fread(p, sizeof(double), 3 * N, fi) != (size_t)(3 * N)) { fprintf(stderr, "short init\n"); return 1; }
    fclose(fi);
    double tau0 = 0;
    if (argc > 15) { /* resume: Erec Qrec tau0 */ Erec = atof(argv[15]); Qrec = atof(argv[16]); tau0 = atof(argv[17]); }
    dst3(q, Qh); dst3(p, Ph);
    double W1 = 1.0 / (2.0 - cbrt(2.0)), W0 = 1.0 - 2.0 * W1;
    double drifts[3], kicks[3], first; int ns;
    if (scheme == 2) { ns = 1; drifts[0] = 1; kicks[0] = 1; first = 0.5; }
    else { ns = 3; drifts[0] = W1; drifts[1] = W0; drifts[2] = W1; kicks[0] = (W1 + W0) / 2; kicks[1] = (W0 + W1) / 2; kicks[2] = W1; first = W1 / 2; }
    double wmax = 0; for (int i = 0; i < 3 * N; i++) if (om[i] > wmax) wmax = om[i];
    double th = wmax * dt * (scheme == 2 ? 1.0 : fabs(W0));
    if (th > 0.9 * M_PI) { fprintf(stderr, "resonance risk theta=%g\n", th); return 2; }
    char fn[1024]; snprintf(fn, sizeof fn, "%s_history.bin", pre);
    FILE *fh = fopen(fn, tau0 > 0 ? "ab" : "wb");
    double row[NCOL];
    double tau = tau0;
    if (tau0 == 0) { diagnostics(0, row); fwrite(row, sizeof(double), NCOL, fh); }
    long every = lround(diag / dt), total = lround((tmax - tau0) / dt), done = 0;
    double E0 = row[1], Q0 = row[2];
    fprintf(stderr, "N=%d theta=%g E0=%.10f Q0=%.10f\n", N, th, E0, Q0);
    double nextck = ckpt > 0 ? (floor(tau0 / ckpt) + 1) * ckpt : 1e300;
    while (done < total) {
        long nb = every < total - done ? every : total - done;
        kick(first * dt);
        for (long i = 0; i < nb; i++) for (int j = 0; j < ns; j++) { drift(drifts[j] * dt); int last = (i == nb - 1) && (j == ns - 1); kick((last ? first : kicks[j]) * dt); }
        done += nb; tau = tau0 + done * dt;
        diagnostics(tau, row); fwrite(row, sizeof(double), NCOL, fh); fflush(fh);
        if ((done / every) % 200 == 0) fprintf(stderr, "tau=%.1f E=%.12f Q=%.12f Erec=%.6f Qrec=%.6f Q60=%.8f\n", tau, row[1], row[2], row[3], row[4], row[13]);
        int save = 0;
        if (strcmp(snaps, "-")) { char buf[4096]; strncpy(buf, snaps, sizeof buf - 1); buf[sizeof buf - 1] = 0;
            for (char *t = strtok(buf, ","); t; t = strtok(NULL, ",")) if (fabs(atof(t) - tau) < 0.5 * dt) save = 1; }
        int ck = tau >= nextck - 0.5 * dt; if (ck) nextck += ckpt;
        if (save || ck) {
            dst3(Qh, q); dst3(Ph, p);
            snprintf(fn, sizeof fn, save ? "%s_snap%ld.bin" : "%s_ckpt.bin", pre, lround(tau));
            FILE *fs = fopen(fn, "wb"); fwrite(q, sizeof(double), 3 * N, fs); fwrite(p, sizeof(double), 3 * N, fs);
            double meta[3] = {tau, Erec, Qrec}; fwrite(meta, sizeof(double), 3, fs); fclose(fs);
        }
    }
    fclose(fh);
    return 0;
}
