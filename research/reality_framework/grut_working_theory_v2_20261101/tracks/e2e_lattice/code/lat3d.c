/* lat3d.c -- track e2e_lattice (Claude, end-to-end pass, GitHub issue #13 context).
 *
 * 3D classical lattice of the HF hybrid after the waterfall: the clock sigma and the real biased waterfall chi
 * (mf_hybrid HFa/HFb potential, rescaling rule of mf_chain), in a self-consistently expanding FRW background, with
 *   - two spectators on the lattice: s_a with m^2 = (g_s chi)^2 (HFa mechanism) and s_b with m^2 = m0^2 + xi_s R(x)
 *     (HFb mechanism), R(x) = -T(x)/M_P^2 the local trace (closure: no metric perturbations);
 *   - a local, Pauli-limited transfer of sigma/chi energy into four axial singlet species (friction Gamma(x) phi_dot^2),
 *     with the half-period production rates of mf_hybrid (axial_tables.json) and the Pauli saturation F_sat(A) of the
 *     mf_hybrid skeptic (sk_axial), applied site by site; blocking energy B_j(x) local (or box-global, option);
 *     N energies/numbers, N -> L H decays and the SM bath are box-averaged fluids;
 *   - tensor perturbations u_ij sourced by grad phi grad phi (Dufaux et al. 0707.0875 method: TT projection at read-out).
 * Every term, closure and parameter is a model choice; nothing follows from closure of GRUT.
 *
 * Units: time and comoving length in H_d^-1 (the design Hubble rate; H_i/H_d = 1 + 2e-5), a = 1 at sigma = sigma_c of the
 * homogeneous attractor.  sigma, chi in units of v; s_a, s_b in units of s_*.  Energy densities rho~ = rho / v^2.
 * Integrator: leapfrog (kick at integer steps with exponential friction, drift at half steps), pi = a^3 phi_dot.
 * a(t) from the Friedmann constraint with a second-order predictor (H_n, Hdot_n); energy budget checked through
 * rho a^3 + int p d(a^3) = const.
 *
 * Build: gcc -O3 -march=native -fopenmp -o lat3d lat3d.c -lm      Run: OMP_NUM_THREADS=2 ./lat3d run.par
 * Outputs (in the run directory): diag.txt (time series), snap_XXX.bin (float32 fields), restart.bin, done.flag.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define NCH 4
#define NACC 96
#define MAXSNAP 400

static int N;
static long NT;
static double L, dx, dt, idx2;
static double *sg, *psg, *ch, *pch, *sa, *psa, *sb, *psb, *chav, *Bj[NCH], *uu[6], *pu[6];
static float *Rx;

/* model (v units) */
static double Lv2, Gv2, sct, cbv, c1v, m2c, mu2, Vc, q, srel, s0t, vH, m0b2, xib, gsa2;
static int loopon;
/* closure */
static int fric, pert_free, Bglobal, own_s, gwon;
static double Gth;
/* tables in u = ln(rho~) */
static int NTAB;
static double tu0, tdu, *tG[NCH], *tF[NCH], *tGb[NCH], *tEi[NCH];
static double Gp[NCH], Grest[NCH], Mch[NCH];
static int wch[NCH];
/* run control */
static double t0, t_end, wall_limit, tau_av, t_fric_on;
static int ndiag, nsnap, snap_gw[MAXSNAP];
static double tsnap[MAXSNAP];
static char rundir[512], initfile[512], tabfile[512];
static int restart;

static double *dalloc(void) {
    double *p = (double *)aligned_alloc(64, ((NT * 8 + 63) / 64) * 64);
    if (!p) { fprintf(stderr, "alloc failed\n"); exit(1); }
    return p;
}

static inline double Vpot(double s, double x) {
    double d = s - sct, x2 = x * x;
    double V = 0.25 * Lv2 * (x2 - 1.0) * (x2 - 1.0) + 0.5 * Gv2 * s * s * x2 - cbv * x + c1v * d - 0.5 * m2c * d * d + Vc;
    if (loopon && s > sct) {
        double M2 = Gv2 * s * s - Lv2;
        if (M2 != 0.0) V += M2 * M2 * (log(fabs(M2) / mu2) - 1.5) / (64.0 * M_PI * M_PI * vH * vH);
    }
    return V;
}
static inline void dVpot(double s, double x, double *Vs, double *Vx) {
    double d = s - sct, x2 = x * x;
    *Vs = Gv2 * s * x2 + c1v - m2c * d;
    if (loopon && s > sct) {
        double M2 = Gv2 * s * s - Lv2;
        if (M2 != 0.0) *Vs += M2 * 2.0 * Gv2 * s * (log(fabs(M2) / mu2) - 1.0) / (32.0 * M_PI * M_PI * vH * vH);
    }
    *Vx = Lv2 * (x2 - 1.0) * x + Gv2 * s * s * x - cbv;
}
static inline double fexpm(double y) { /* exp(-y), y >= 0, Taylor for small y */
    if (y < 0.02) return 1.0 - y * (1.0 - y * (0.5 - y * (1.0 / 6.0 - y * (1.0 / 24.0))));
    return exp(-y);
}
static inline double tab(const double *T, double u) {
    double f = (u - tu0) / tdu;
    if (f <= 0) return T[0];
    if (f >= NTAB - 1) return T[NTAB - 1];
    int k = (int)f;
    double w = f - k;
    return T[k] * (1 - w) + T[k + 1] * w;
}

/* ---------------------------------------------------------------- parameter file */
static void read_par(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "no par %s\n", fn); exit(1); }
    char key[128], line[4096];
    nsnap = 0;
    own_s = 0; gwon = 1; Bglobal = 0; pert_free = 1; restart = 0; tau_av = 0.05; t_fric_on = 0.0; loopon = 1;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || strlen(line) < 2) continue;
        if (sscanf(line, "%127s", key) != 1) continue;
        char *val = line + strlen(key);
#define RD(name, var) if (!strcmp(key, name)) { var = atof(val); continue; }
#define RI(name, var) if (!strcmp(key, name)) { var = atoi(val); continue; }
        RI("N", N) RD("L", L) RD("dt", dt) RD("t0", t0) RD("t_end", t_end) RD("wall_limit", wall_limit)
        RI("ndiag", ndiag) RD("Lv2", Lv2) RD("Gv2", Gv2) RD("sct", sct) RD("cbv", cbv) RD("c1v", c1v) RD("m2c", m2c)
        RD("mu2", mu2) RD("Vc", Vc) RD("q", q) RD("srel", srel) RD("s0t", s0t) RD("vH", vH) RD("m0b2", m0b2) RD("xib", xib)
        RD("gsa2", gsa2) RI("loop", loopon) RI("fric", fric) RI("pert_free", pert_free) RI("Bglobal", Bglobal) RD("Gth", Gth)
        RI("own_s", own_s) RI("gw", gwon) RD("tau_av", tau_av) RD("t_fric_on", t_fric_on) RI("restart", restart)
        if (!strcmp(key, "rundir")) { sscanf(val, "%511s", rundir); continue; }
        if (!strcmp(key, "init")) { sscanf(val, "%511s", initfile); continue; }
        if (!strcmp(key, "tables")) { sscanf(val, "%511s", tabfile); continue; }
        if (!strcmp(key, "snap")) {
            double ts; int g = 0;
            if (sscanf(val, "%lf %d", &ts, &g) >= 1 && nsnap < MAXSNAP) { tsnap[nsnap] = ts; snap_gw[nsnap] = g; nsnap++; }
            continue;
        }
        fprintf(stderr, "unknown key %s\n", key);
    }
    fclose(f);
}

static void read_tables(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "no tables %s\n", fn); exit(1); }
    if (fscanf(f, "%d %lf %lf", &NTAB, &tu0, &tdu) != 3) exit(2);
    for (int j = 0; j < NCH; j++) {
        if (fscanf(f, "%d %lf %lf %lf", &wch[j], &Gp[j], &Grest[j], &Mch[j]) != 4) exit(2);
        tG[j] = malloc(NTAB * 8); tF[j] = malloc(NTAB * 8); tGb[j] = malloc(NTAB * 8); tEi[j] = malloc(NTAB * 8);
        for (int k = 0; k < NTAB; k++)
            if (fscanf(f, "%lf %lf %lf %lf", &tG[j][k], &tF[j][k], &tGb[j][k], &tEi[j][k]) != 4) exit(2);
    }
    fclose(f);
}

/* ---------------------------------------------------------------- state scalars */
typedef struct {
    double t, a, anm, Hprev, Wp, rho_prev, a_prev;
    double rN[NCH], nN[NCH], rSM, Bmean[NCH], capmean[NCH], Emean[2];
    long step;
    int isnap;
    double rho0a3;
} State;
static State S;

static void write_restart(const char *fn) {
    FILE *f = fopen(fn, "wb");
    fwrite(&S, sizeof S, 1, f);
    double *arr[] = {sg, psg, ch, pch, sa, psa, sb, psb, chav};
    for (int k = 0; k < 9; k++) fwrite(arr[k], 8, NT, f);
    for (int j = 0; j < NCH; j++) fwrite(Bj[j], 8, NT, f);
    for (int c = 0; c < 6; c++) { fwrite(uu[c], 8, NT, f); fwrite(pu[c], 8, NT, f); }
    fclose(f);
}
static int read_restart(const char *fn) {
    FILE *f = fopen(fn, "rb");
    if (!f) return 0;
    size_t r = fread(&S, sizeof S, 1, f);
    double *arr[] = {sg, psg, ch, pch, sa, psa, sb, psb, chav};
    for (int k = 0; k < 9; k++) r += fread(arr[k], 8, NT, f);
    for (int j = 0; j < NCH; j++) r += fread(Bj[j], 8, NT, f);
    for (int c = 0; c < 6; c++) { r += fread(uu[c], 8, NT, f); r += fread(pu[c], 8, NT, f); }
    fclose(f);
    return 1;
}

static void write_snap(int k, double a) {
    char fn[640];
    snprintf(fn, sizeof fn, "%s/snap_%03d.bin", rundir, k);
    FILE *f = fopen(fn, "wb");
    int hdr[2] = {N, snap_gw[k] && gwon};
    double sc[3] = {S.t, a, L};
    fwrite(hdr, 4, 2, f);
    fwrite(sc, 8, 3, f);
    float *buf = malloc(NT * 4);
    double *arr[] = {ch, chav, sg, sa, sb};
    for (int m = 0; m < 5; m++) { for (long i = 0; i < NT; i++) buf[i] = (float)arr[m][i]; fwrite(buf, 4, NT, f); }
    fwrite(Rx, 4, NT, f);
    double ia3 = 1.0 / (a * a * a);
    for (long i = 0; i < NT; i++) buf[i] = (float)(pch[i] * ia3);
    fwrite(buf, 4, NT, f);
    if (snap_gw[k] && gwon) {
        for (int c = 0; c < 6; c++) { for (long i = 0; i < NT; i++) buf[i] = (float)(pu[c][i] * ia3); fwrite(buf, 4, NT, f); }
    }
    free(buf);
    fclose(f);
}

/* accumulator slots */
enum { A_KS, A_KX, A_GS, A_GX, A_V, A_ES, A_EX, A_SKIN, A_SGRAD, A_SPOT, A_INJ0, A_INJN0 = A_INJ0 + NCH,
       A_CAP0 = A_INJN0 + NCH, A_SDOT = A_CAP0 + NCH, A_XDOT, A_SG, A_XX,
       /* diag-only */
       A_R, A_R2, A_X2, A_XAV, A_S2, A_SA, A_SA2, A_SB, A_SB2, A_MA, A_MA2, A_MB, A_MB2, A_FM, A_FMAV, A_B0,
       A_GAMS = A_B0 + NCH, A_GAMX, A_HH, A_XD2, A_XABS, A_GAMXW, A_KXW, A_GB0, A_NONE = A_GB0 + NCH };

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: lat3d run.par\n"); return 1; }
    read_par(argv[1]);
    read_tables(tabfile);
    NT = (long)N * N * N;
    dx = L / N;
    idx2 = 1.0 / (dx * dx);
    sg = dalloc(); psg = dalloc(); ch = dalloc(); pch = dalloc(); sa = dalloc(); psa = dalloc(); sb = dalloc(); psb = dalloc();
    chav = dalloc();
    for (int j = 0; j < NCH; j++) Bj[j] = dalloc();
    for (int c = 0; c < 6; c++) { uu[c] = dalloc(); pu[c] = dalloc(); }
    Rx = malloc(NT * 4);
    char rfn[640];
    snprintf(rfn, sizeof rfn, "%s/restart.bin", rundir);
    int resumed = 0;
    if (restart) resumed = read_restart(rfn);
    if (!resumed) {
        FILE *f = fopen(initfile, "rb");
        if (!f) { fprintf(stderr, "no init %s\n", initfile); return 1; }
        double hdr[2];
        size_t r = fread(hdr, 8, 2, f);
        double *arr[] = {sg, psg, ch, pch, sa, psa, sb, psb};
        for (int k = 0; k < 8; k++) r += fread(arr[k], 8, NT, f);
        fclose(f);
        if (r != 2 + 8 * (size_t)NT) { fprintf(stderr, "init short read\n"); return 1; }
        memset(&S, 0, sizeof S);
        S.t = t0; S.a = hdr[0]; S.anm = hdr[0]; S.a_prev = hdr[0]; S.Hprev = hdr[1];
        for (long i = 0; i < NT; i++) { chav[i] = ch[i]; for (int j = 0; j < NCH; j++) Bj[j][i] = 0.0; }
        for (int c = 0; c < 6; c++) for (long i = 0; i < NT; i++) { uu[c][i] = 0.0; pu[c][i] = 0.0; }
        for (int j = 0; j < NCH; j++) { S.rN[j] = 0; S.nN[j] = 1e-300; S.rSM = 0; }
        while (S.isnap < nsnap && tsnap[S.isnap] < t0 - 0.5 * dt) S.isnap++;
    }
    char dfn[640];
    snprintf(dfn, sizeof dfn, "%s/diag.txt", rundir);
    FILE *fd = fopen(dfn, resumed ? "a" : "w");
    if (!resumed)
        fprintf(fd, "# t a H Hdot eps w rho fsig fchi fgrad fN1 fN2 fN3 fN4 fSM fs R_mean R_rms chi_mean chi_std chiav_mean "
                    "sig_mean sig_std fminus fminus_av area_av area_inst sa_mean sa_std sb_mean sb_std ms2a_mean ms2a_std "
                    "ms2b_mean ms2b_std B1 B2 B3 B4 cap1 cap2 cap3 cap4 Gam_sig Gam_chi budget hdot_uu fcoh xdot_std "
                    "chiabs_mean Gam_chi_w rN1 rN2 rN3 rN4 rSM nN1 nN2 nN3 nN4 Ks Kx Gs Gx Vm\n");
    double tstart_wall = omp_get_wtime();
    int exitcode = 0;
    double acc[NACC];

    while (S.t < t_end - 0.5 * dt) {
        int dg = (S.step % ndiag == 0);
        int dosnap = (S.isnap < nsnap && S.t >= tsnap[S.isnap] - 0.5 * dt);
        int diag = dg || dosnap;
        double a = S.a, a3 = a * a * a, anm3 = S.anm * S.anm * S.anm;
        double ia2 = 1.0 / (a * a), ianm3 = 1.0 / anm3;
        double H = S.Hprev;
        int fon = (fric > 0) && (S.t >= t_fric_on - 0.5 * dt);
        double Bm[NCH], Cm[NCH];
        for (int j = 0; j < NCH; j++) {
            Bm[j] = S.Bmean[j]; Cm[j] = S.capmean[j];
            if (Bglobal == 2) { /* skeptic's saturation on the box-mean condensate energy (homogeneous F_sat(A) model) */
                double Em = S.Emean[wch[j]];
                Cm[j] = Em > 1e-300 ? tab(tF[j], log(Em)) * Em : 0.0;
            }
        }
        for (int k = 0; k < NACC; k++) acc[k] = 0.0;

#pragma omp parallel for collapse(2) schedule(static) reduction(+ : acc[:NACC])
        for (int x = 0; x < N; x++)
            for (int y = 0; y < N; y++) {
                int xp = (x + 1) % N, xm = (x + N - 1) % N, yp = (y + 1) % N, ym = (y + N - 1) % N;
                for (int z = 0; z < N; z++) {
                    int zp = (z + 1) % N, zm = (z + N - 1) % N;
                    long i = ((long)x * N + y) * N + z;
                    long ixp = ((long)xp * N + y) * N + z, ixm = ((long)xm * N + y) * N + z;
                    long iyp = ((long)x * N + yp) * N + z, iym = ((long)x * N + ym) * N + z;
                    long izp = ((long)x * N + y) * N + zp, izm = ((long)x * N + y) * N + zm;
                    double s = sg[i], c = ch[i];
                    double lap_s = (sg[ixp] + sg[ixm] + sg[iyp] + sg[iym] + sg[izp] + sg[izm] - 6 * s) * idx2;
                    double lap_c = (ch[ixp] + ch[ixm] + ch[iyp] + ch[iym] + ch[izp] + ch[izm] - 6 * c) * idx2;
                    /* forward-difference gradient energy (consistent with the 7-point Laplacian) */
                    double fs1 = sg[ixp] - s, fs2 = sg[iyp] - s, fs3 = sg[izp] - s;
                    double fc1 = ch[ixp] - c, fc2 = ch[iyp] - c, fc3 = ch[izp] - c;
                    double gs = 0.5 * (fs1 * fs1 + fs2 * fs2 + fs3 * fs3) * idx2 * ia2;
                    double gx = 0.5 * (fc1 * fc1 + fc2 * fc2 + fc3 * fc3) * idx2 * ia2;
                    double V = Vpot(s, c), Vs, Vx;
                    dVpot(s, c, &Vs, &Vx);
                    double ps0 = psg[i], pc0 = pch[i];
                    /* local rates (phi_dot lagged by dt/2) */
                    double Gs_ = 0.0, Gx_ = 0.0, Gf[NCH], Ei[NCH], Gbv[NCH];
                    double Es = 0, Ex = 0;
                    {
                        double sdl = ps0 * ianm3, xdl = pc0 * ianm3;
                        double Et = 0.5 * (sdl * sdl + xdl * xdl) + V;
                        Es = 0.5 * sdl * sdl + 0.5 * Gv2 * c * c * (s - s0t) * (s - s0t);
                        if (Es < 0) Es = 0;
                        Ex = Et - Es;
                        if (Ex < 0) Ex = 0;
                    }
                    for (int j = 0; j < NCH; j++) { Gf[j] = 0.0; Ei[j] = 1.0; Gbv[j] = 0.0; }
                    if (fon) {
                        double us = Es > 1e-300 ? log(Es) : -700, ux = Ex > 1e-300 ? log(Ex) : -700;
                        for (int j = 0; j < NCH; j++) {
                            double E = wch[j] ? Ex : Es, uE = wch[j] ? ux : us;
                            if (E <= 1e-300) continue;
                            double G = tab(tG[j], uE);
                            Ei[j] = tab(tEi[j], uE);
                            Gbv[j] = tab(tGb[j], uE);
                            acc[A_GB0 + j] += Gbv[j];
                            double bl = 1.0;
                            if (fric == 2) {
                                double cap = tab(tF[j], uE) * E;
                                acc[A_CAP0 + j] += cap;
                                if (Bglobal) bl = Cm[j] > 0 ? 1.0 - Bm[j] / Cm[j] : 0.0;
                                else bl = cap > 0 ? 1.0 - Bj[j][i] / cap : 0.0;
                                if (bl < 0) bl = 0;
                            }
                            double gnp = G * bl, gp = Gp[j] * ((pert_free || fric != 2) ? 1.0 : bl);
                            Gf[j] = gnp > gp ? gnp : gp;
                            if (wch[j]) Gx_ += Gf[j]; else Gs_ += Gf[j];
                        }
                    }
                    double Fs = a * lap_s - a3 * Vs, Fx = a * lap_c - a3 * Vx;
                    double ds = fexpm(Gs_ * dt), dsh = fexpm(0.5 * Gs_ * dt), dxx = fexpm(Gx_ * dt), dxh = fexpm(0.5 * Gx_ * dt);
                    double ps1 = ps0 * ds + dt * Fs * dsh, pc1 = pc0 * dxx + dt * Fx * dxh;
                    psg[i] = ps1; pch[i] = pc1;
                    double sd = 0.5 * (ps0 + ps1) / a3, xd = 0.5 * (pc0 + pc1) / a3;
                    acc[A_KS] += 0.5 * sd * sd; acc[A_KX] += 0.5 * xd * xd; acc[A_GS] += gs; acc[A_GX] += gx; acc[A_V] += V;
                    acc[A_SDOT] += sd; acc[A_XDOT] += xd; acc[A_SG] += s; acc[A_XX] += c;
                    {
                        double Esn = 0.5 * sd * sd + 0.5 * Gv2 * c * c * (s - s0t) * (s - s0t);
                        acc[A_ES] += Esn; acc[A_EX] += 0.5 * (sd * sd + xd * xd) + V - Esn;
                    }
                    if (fon) {
                        for (int j = 0; j < NCH; j++) {
                            double vel2 = wch[j] ? xd * xd : sd * sd;
                            double dE = Gf[j] * vel2 * dt;
                            acc[A_INJ0 + j] += dE;
                            acc[A_INJN0 + j] += dE / Ei[j];
                            if (fric == 2) {
                                double bb = Bj[j][i] * fexpm((6.0 * H + Gbv[j] + Gth) * dt) + dE;
                                Bj[j][i] = bb;
                                if (diag) acc[A_B0 + j] += bb;
                            }
                        }
                    }
                    /* local Ricci scalar R = -T/M_P^2 (H^2 units) */
                    double Rl = q * (-(sd * sd + xd * xd) + 2.0 * (gs + gx) + 4.0 * V);
                    /* spectators */
                    double sav = sa[i], sbv = sb[i];
                    double lap_a = (sa[ixp] + sa[ixm] + sa[iyp] + sa[iym] + sa[izp] + sa[izm] - 6 * sav) * idx2;
                    double lap_b = (sb[ixp] + sb[ixm] + sb[iyp] + sb[iym] + sb[izp] + sb[izm] - 6 * sbv) * idx2;
                    double msa = gsa2 * c * c, msb = m0b2 + xib * Rl;
                    double pa0 = psa[i], pb0 = psb[i];
                    double pa1 = pa0 + dt * (a * lap_a - a3 * msa * sav), pb1 = pb0 + dt * (a * lap_b - a3 * msb * sbv);
                    psa[i] = pa1; psb[i] = pb1;
                    {
                        double sdd = own_s == 1 ? 0.5 * (pa0 + pa1) / a3 : 0.5 * (pb0 + pb1) / a3;
                        double sv = own_s == 1 ? sav : sbv, ms = own_s == 1 ? msa : m0b2;
                        double *arr = own_s == 1 ? sa : sb;
                        double f1 = arr[ixp] - sv, f2 = arr[iyp] - sv, f3 = arr[izp] - sv;
                        if (own_s) {
                            acc[A_SKIN] += srel * 0.5 * sdd * sdd;
                            acc[A_SGRAD] += srel * 0.5 * (f1 * f1 + f2 * f2 + f3 * f3) * idx2 * ia2;
                            acc[A_SPOT] += srel * 0.5 * ms * sv * sv;
                        }
                    }
                    /* tensor perturbations: pi_u += dt (a lap u + 2 q a d_i phi d_j phi), central gradients */
                    if (gwon) {
                        double cs1 = 0.5 * (sg[ixp] - sg[ixm]) / dx, cs2 = 0.5 * (sg[iyp] - sg[iym]) / dx, cs3 = 0.5 * (sg[izp] - sg[izm]) / dx;
                        double cc1 = 0.5 * (ch[ixp] - ch[ixm]) / dx, cc2 = 0.5 * (ch[iyp] - ch[iym]) / dx, cc3 = 0.5 * (ch[izp] - ch[izm]) / dx;
                        double src[6] = {cs1 * cs1 + cc1 * cc1, cs2 * cs2 + cc2 * cc2, cs3 * cs3 + cc3 * cc3,
                                         cs1 * cs2 + cc1 * cc2, cs1 * cs3 + cc1 * cc3, cs2 * cs3 + cc2 * cc3};
                        double hh = 0.0;
                        for (int m = 0; m < 6; m++) {
                            double *U = uu[m];
                            double lu = (U[ixp] + U[ixm] + U[iyp] + U[iym] + U[izp] + U[izm] - 6 * U[i]) * idx2;
                            double p0 = pu[m][i], p1 = p0 + dt * (a * lu + 2.0 * q * a * src[m]);
                            pu[m][i] = p1;
                            if (diag) { double hd = 0.5 * (p0 + p1) / a3; hh += (m < 3 ? 1.0 : 2.0) * hd * hd; }
                        }
                        if (diag) acc[A_HH] += hh;
                    }
                    if (diag) {
                        Rx[i] = (float)Rl;
                        acc[A_R] += Rl; acc[A_R2] += Rl * Rl; acc[A_X2] += c * c; acc[A_XAV] += chav[i]; acc[A_S2] += s * s;
                        acc[A_SA] += sav; acc[A_SA2] += sav * sav; acc[A_SB] += sbv; acc[A_SB2] += sbv * sbv;
                        acc[A_MA] += msa; acc[A_MA2] += msa * msa; acc[A_MB] += msb; acc[A_MB2] += msb * msb;
                        acc[A_FM] += (c < 0); acc[A_FMAV] += (chav[i] < 0);
                        acc[A_GAMS] += Gs_; acc[A_GAMX] += Gx_; acc[A_XD2] += xd * xd; acc[A_XABS] += fabs(c);
                        acc[A_GAMXW] += Gx_ * xd * xd; acc[A_KXW] += xd * xd;
                    }
                }
            }
        /* box means */
        double inv = 1.0 / (double)NT;
        for (int k = 0; k < NACC; k++) acc[k] *= inv;
        double rF[NCH], rN = 0;
        for (int j = 0; j < NCH; j++) { rF[j] = S.rN[j]; rN += S.rN[j]; }
        double rrad = rN + S.rSM;
        double Kin = acc[A_KS] + acc[A_KX], Grd = acc[A_GS] + acc[A_GX];
        double rho = Kin + Grd + acc[A_V] + rrad + acc[A_SKIN] + acc[A_SGRAD] + acc[A_SPOT];
        double p = Kin - Grd / 3.0 - acc[A_V] + rrad / 3.0 + acc[A_SKIN] - acc[A_SGRAD] / 3.0 - acc[A_SPOT];
        H = sqrt(q * rho / 3.0);
        double Hd = -0.5 * q * (rho + p);
        /* fluids */
        for (int j = 0; j < NCH; j++) {
            double Em = S.nN[j] > 1e-299 ? S.rN[j] / S.nN[j] : Mch[j];
            double Gl = Grest[j] * Mch[j] / (Mch[j] > Em ? Mch[j] : Em);
            double dec = Gl * S.rN[j] * dt;
            S.rN[j] = S.rN[j] * exp(-(4.0 * H + Gl) * dt) + acc[A_INJ0 + j];
            S.nN[j] = S.nN[j] * exp(-(3.0 * H + Gl) * dt) + acc[A_INJN0 + j];
            S.rSM += dec;
            if (fric == 2) S.Bmean[j] = Bglobal ? S.Bmean[j] * exp(-(6.0 * H + acc[A_GB0 + j] + Gth) * dt) + acc[A_INJ0 + j] : 0.0;
            S.capmean[j] = Bglobal == 2 ? Cm[j] : acc[A_CAP0 + j];
        }
        S.Emean[0] = acc[A_ES]; S.Emean[1] = acc[A_EX] > 0 ? acc[A_EX] : 0.0;
        S.rSM *= exp(-4.0 * H * dt);
        /* scale factor */
        double an = S.a;
        double ah = an * exp(0.5 * H * dt + 0.125 * Hd * dt * dt);
        double a1 = an * exp(H * dt + 0.5 * Hd * dt * dt);
        if (S.step == 0 && !resumed) S.rho0a3 = rho * an * an * an;
        S.Wp += 0.5 * p * (a1 * a1 * a1 - S.a_prev * S.a_prev * S.a_prev);
        double budget = (rho * an * an * an + S.Wp - 0.5 * p * (a1 * a1 * a1 - an * an * an)) / S.rho0a3 - 1.0;
        if (diag) {
            /* wall area from sign changes of chi_av and chi on links (estimator (2/3) n_cross dx^2 / V, comoving) */
            double nav = 0, nin = 0;
#pragma omp parallel for collapse(2) schedule(static) reduction(+ : nav, nin)
            for (int x = 0; x < N; x++)
                for (int y = 0; y < N; y++) {
                    int xp = (x + 1) % N, yp = (y + 1) % N;
                    for (int z = 0; z < N; z++) {
                        int zp = (z + 1) % N;
                        long i = ((long)x * N + y) * N + z;
                        long j1 = ((long)xp * N + y) * N + z, j2 = ((long)x * N + yp) * N + z, j3 = ((long)x * N + y) * N + zp;
                        int s0 = chav[i] < 0;
                        nav += (s0 != (chav[j1] < 0)) + (s0 != (chav[j2] < 0)) + (s0 != (chav[j3] < 0));
                        int s1 = ch[i] < 0;
                        nin += (s1 != (ch[j1] < 0)) + (s1 != (ch[j2] < 0)) + (s1 != (ch[j3] < 0));
                    }
                }
            double area_av = (2.0 / 3.0) * nav / ((double)NT * dx);
            double area_in = (2.0 / 3.0) * nin / ((double)NT * dx);
            double xm = acc[A_XX], sm_ = acc[A_SG];
            double Rm = acc[A_R], Rs = sqrt(fmax(acc[A_R2] - Rm * Rm, 0));
            double E0 = 0.5 * (acc[A_SDOT] * acc[A_SDOT] + acc[A_XDOT] * acc[A_XDOT]) + Vpot(sm_, xm);
            double fcoh = E0 / (Kin + Grd + acc[A_V]);
            double Gxw = acc[A_KXW] > 0 ? acc[A_GAMXW] / acc[A_KXW] : 0;
            fprintf(fd, "%.8g %.10g %.10g %.6g %.6g %.6g %.10g", S.t, an, H, Hd, -Hd / (H * H), p / rho, rho);
            fprintf(fd, " %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g", (acc[A_ES]) / rho, (acc[A_EX]) / rho, Grd / rho,
                    rF[0] / rho, rF[1] / rho, rF[2] / rho, rF[3] / rho, S.rSM / rho, (acc[A_SKIN] + acc[A_SGRAD] + acc[A_SPOT]) / rho);
            fprintf(fd, " %.8g %.6g %.8g %.6g %.8g %.8g %.6g %.6g %.6g %.6g %.6g", Rm, Rs, xm, sqrt(fmax(acc[A_X2] - xm * xm, 0)),
                    acc[A_XAV], sm_, sqrt(fmax(acc[A_S2] - sm_ * sm_, 0)), acc[A_FM], acc[A_FMAV], area_av, area_in);
            fprintf(fd, " %.10g %.6g %.10g %.6g %.8g %.6g %.8g %.6g", acc[A_SA], sqrt(fmax(acc[A_SA2] - acc[A_SA] * acc[A_SA], 0)),
                    acc[A_SB], sqrt(fmax(acc[A_SB2] - acc[A_SB] * acc[A_SB], 0)), acc[A_MA], sqrt(fmax(acc[A_MA2] - acc[A_MA] * acc[A_MA], 0)),
                    acc[A_MB], sqrt(fmax(acc[A_MB2] - acc[A_MB] * acc[A_MB], 0)));
            for (int j = 0; j < NCH; j++) fprintf(fd, " %.6g", (fric == 2 ? (Bglobal ? S.Bmean[j] : acc[A_B0 + j]) : 0.0) / rho);
            for (int j = 0; j < NCH; j++) fprintf(fd, " %.6g", acc[A_CAP0 + j] / rho);
            fprintf(fd, " %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g", acc[A_GAMS], acc[A_GAMX], budget, acc[A_HH], fcoh,
                    sqrt(fmax(acc[A_XD2] - acc[A_XDOT] * acc[A_XDOT], 0)), acc[A_XABS], Gxw);
            for (int j = 0; j < NCH; j++) fprintf(fd, " %.6g", rF[j]);
            fprintf(fd, " %.6g", S.rSM);
            for (int j = 0; j < NCH; j++) fprintf(fd, " %.6g", S.nN[j]);
            fprintf(fd, " %.6g %.6g %.6g %.6g %.6g\n", acc[A_KS], acc[A_KX], acc[A_GS], acc[A_GX], acc[A_V]);
            fflush(fd);
        }
        if (dosnap) { write_snap(S.isnap, an); S.isnap++; }
        /* drift */
        double iah3 = 1.0 / (ah * ah * ah), wav = dt / tau_av;
#pragma omp parallel for schedule(static)
        for (long i = 0; i < NT; i++) {
            sg[i] += dt * psg[i] * iah3;
            ch[i] += dt * pch[i] * iah3;
            sa[i] += dt * psa[i] * iah3;
            sb[i] += dt * psb[i] * iah3;
            if (gwon) for (int m = 0; m < 6; m++) uu[m][i] += dt * pu[m][i] * iah3;
            chav[i] += wav * (ch[i] - chav[i]);
        }
        S.a_prev = an;
        S.anm = ah;
        S.a = a1;
        S.Hprev = H;
        S.rho_prev = rho;
        S.t += dt;
        S.step++;
        if (S.step % 50 == 0 && omp_get_wtime() - tstart_wall > wall_limit) { exitcode = 3; break; }
    }
    fclose(fd);
    write_restart(rfn);
    if (exitcode == 0) {
        char ffn[640];
        snprintf(ffn, sizeof ffn, "%s/done.flag", rundir);
        FILE *ff = fopen(ffn, "w");
        fprintf(ff, "%.10g %ld %.3f\n", S.t, S.step, omp_get_wtime() - tstart_wall);
        fclose(ff);
    }
    printf("t=%.6g step=%ld a=%.6g wall=%.1fs exit=%d\n", S.t, S.step, S.a, omp_get_wtime() - tstart_wall, exitcode);
    return exitcode;
}
