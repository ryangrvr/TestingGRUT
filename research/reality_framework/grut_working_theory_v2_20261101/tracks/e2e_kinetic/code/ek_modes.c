/* ek_modes.c -- e2e_kinetic (Claude, end-to-end pass, GitHub issue #13 context).  Own code.
 *
 * Energy transfer from the oscillating sigma / chi condensates of HFa/HFb into the Standard Model (SM) bath through
 * the axial singlet channels -(M_j/2) Nbar_j exp(i gamma5 phi/f_j) N_j, with
 *   (1) Pauli-blocked fermion production computed mode by mode (two-level system per comoving momentum q and
 *       helicity h, rotated frame: H = b.sigma, b = (M, 0, h q/a - mu5(t)), mu5 = phidot/2f = -A m sin(m(t - t_TP)));
 *   (2) N removal by the seesaw decay N -> L H (vacuum width with thermal-mass kinematics) and by in-medium
 *       N + SM -> SM + SM scatterings (Besak-Bodeker 1202.1288 leading-order total, damping-rate form
 *       Gamma = c h^2 Tcal^2 / omega), with detailed balance against the thermal part of the SM bath;
 *   (3) the SM cascade: a hard (non-thermal) component (rho_h, n_h) injected by the N removal, which relaxes into a
 *       thermal component rho_th at the Kurkela-Lu (1405.6318) under-occupied equilibration rate
 *       1/t_eq, t_eq = (34 + 21 ln(Q/T)) (Q/T)^{1/2} / [(1 + 0.037 ln(1/lam)) lam^2 T], T the temperature of the whole
 *       SM energy, Q = rho_h/n_h.  (Kinetic-theory result with LPM-suppressed splitting; Harigaya-Mukaida 1312.3097 give
 *       the same (Q/T)^{1/2}/(alpha^2 T) scaling parametrically.)
 * Background: flat FRW, H^2 = rho/(3 mP^2); condensates are harmonic, cycle-averaged energy rho_P (w = 0), amplitude
 * Sigma = sqrt(2 rho_P)/m_P, A_j = Sigma/(2 f_j) (mu5,max in units of m_P).
 *
 * Two modes of N evolution (every choice is a model choice, see NOTES.md):
 *   COHERENT: exact Bloch-vector rotation r' = R(2|b|dt) r for every active mode (unitary, exact for piecewise constant
 *             b at the step midpoint), collisions applied at every turning point of the parent (mu5 = 0, where the
 *             free basis is the instantaneous eigenbasis): n -> n e^{-G dt} + f_eff (1 - e^{-G dt}), the transverse
 *             (coherence) part damped by e^{-G dt/2}.  Back-reaction: the condensate loses the free-basis N energy
 *             produced between turning points.  Modes above k_cut = (1.3 A + 4 M/m + 2) m are adiabatic: frozen
 *             (occupation kept, coherence dropped), still subject to collisions.
 *   DEPHASED: rate equation per mode, dn/dt = R (1 - 2n) - G n + G_eq f_eq, R = (m/2 pi) pp(k; A), pp = p_- + p_+ the
 *             single half-period transition probabilities of both helicities from the free vacuum (table mode below);
 *             for the resonant channels (2M < m) at A < A_pert the golden-rule shell at k_res = (m/2) beta.
 *   TABLE:    half-period transition probabilities pp(k; A, M/m) (no expansion, constant amplitude), units m = 1.
 *
 * Usage:  ek_modes table OUTFILE Mh A1 A2 ... (prints k, p_minus, p_plus)
 *         ek_modes run INFILE OUTPREFIX     (input format: see write_input() in ek_common.py)
 * Units: H_i = 1 (H at sigma = sigma_c), a = 1 there.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NSP 4
#define TWOPI2 (2.0 * M_PI * M_PI)

typedef struct {
    int par;                 /* parent: 0 sigma, 1 chi */
    double M, f, Grest, h2;  /* mass, decay constant, vacuum rest width, sum_i |h_ji|^2 */
    double kres_m, Gpert;    /* resonance k/m for 2M < m (else 0), perturbative width (M/f)^2 m beta/16pi */
    int nq;
    double *q, *w;           /* comoving momenta, weights q^2 dq / (2 pi^2) (per helicity) */
    double *r;               /* Bloch vectors [ (h*nq + i)*3 + c ] (coherent) */
    double *n;               /* occupations [h*nq + i] (free basis, after the last collision step) */
    int *frozen;             /* [h*nq + i] */
    int nA, nx;
    double *lA, *pp;         /* table: log A grid, pp[iA*nx + ix] on x in [0,1], k = x (1.25 A + 3 M/m + 2) m */
    double E, Nn;            /* physical energy and number density of the species (last evaluation) */
    double Pprod_acc;        /* energy produced since the last record (per comoving volume a^3) */
    double Gband, Gall;      /* occupation-weighted removal rates (band, all) at the last collision step */
    double nband;            /* mean occupation in the band k < A m (both helicities) */
    double Ecoll_acc;        /* energy moved N -> SM since the last record (comoving) */
} Spec;

/* ------------------------------------------------------------------ global state and parameters */
static int coherent = 1, noexp = 0, fixamp = 0, nocoll = 0;
static double mP, mpar[2], rhoP[2];
static double tph_g[2] = {0, 0};      /* phase reference of the parents: mu5 = -A m sin(m (t - tph)) */
static double tnext_in[2] = {-1, -1}, tlastc_in[2] = {-1, -1};   /* coherent-stage restart (negative: fresh start) */
static double t, a, H;
static double rho_h, n_h, rho_th;      /* SM: hard energy and number, thermal energy (physical) */
static double lam_eff, gstar, kH, kL, cBB, wthermal_inv = 1.0, KLc0 = 34.0, KLc1 = 21.0, KLl = 0.037, KLscale = 1.0;
static double gsum_f = 73.0;           /* g_b + g_f/2 of the SM (Tcal^2 normalisation) */
static double Apert = 1.0, cdt = 0.04, ckcut = 1.3;
#ifndef LCUT
#define LCUT 40.0               /* Lorentzian shell cut |x| < LCUT (model choice; -DLCUT=... variants) */
#endif
static double dshell = 0.0;            /* log width of the perturbative shell; 0: one bin; >0 Gaussian sigma; <0 Lorentzian HWHM */
static double inj_E = 0, inj_N = 0, inv_E = 0;   /* SM injections since the last SM update (physical) */
static double t_end;
static Spec S[NSP];
static FILE *frec;
static double Ecoll_tot = 0, Eprod_tot = 0;    /* comoving bookkeeping */

static double Tof(double rho) { return rho > 0 ? pow(30.0 * rho / (M_PI * M_PI * gstar), 0.25) : 0.0; }
static double Tcal2(void) {
    double T = Tof(rho_th), hard = (rho_h > 1e-300 && n_h > 0) ? (12.0 / gsum_f) * n_h * n_h / rho_h : 0.0;
    return T * T + hard;
}
static double Phi_kin(double M, double Tc2) {
    double mL2 = kL * Tc2, mH2 = kH * Tc2;
    if (M <= sqrt(mL2) + sqrt(mH2)) return 0.0;
    double xL = mL2 / (M * M), xH = mH2 / (M * M);
    double lam = (1 - xL - xH) * (1 - xL - xH) - 4 * xL * xH;
    return (1 + xL - xH) * sqrt(lam > 0 ? lam : 0);
}
static double tau_KL(void) {
    double rs = rho_h + rho_th;
    if (rho_h <= 0 || n_h <= 0 || rs <= 0) return 1e300;
    double T = Tof(rs), Q = rho_h / n_h, x = Q / T;
    if (x < 1) x = 1;
    double lam = lam_eff;
    return KLscale * (KLc0 + KLc1 * log(x)) / (1.0 + KLl * log(1.0 / lam)) * sqrt(x) / (lam * lam * T);
}

/* ------------------------------------------------------------------ input */
static void die(const char *m) { fprintf(stderr, "ek_modes: %s\n", m); exit(1); }
static double rd(FILE *fp) { double x; if (fscanf(fp, "%lf", &x) != 1) die("read"); return x; }
static int ri(FILE *fp) { int x; if (fscanf(fp, "%d", &x) != 1) die("readi"); return x; }

static void read_input(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("open input");
    coherent = ri(fp); noexp = ri(fp); fixamp = ri(fp); nocoll = ri(fp);
    mP = rd(fp); t = rd(fp); t_end = rd(fp); a = rd(fp);
    tph_g[0] = rd(fp); tph_g[1] = rd(fp); tnext_in[0] = rd(fp); tnext_in[1] = rd(fp); tlastc_in[0] = rd(fp); tlastc_in[1] = rd(fp);
    mpar[0] = rd(fp); mpar[1] = rd(fp); rhoP[0] = rd(fp); rhoP[1] = rd(fp);
    rho_h = rd(fp); n_h = rd(fp); rho_th = rd(fp);
    lam_eff = rd(fp); gstar = rd(fp); kH = rd(fp); kL = rd(fp); cBB = rd(fp); gsum_f = rd(fp);
    KLc0 = rd(fp); KLc1 = rd(fp); KLl = rd(fp); KLscale = rd(fp); wthermal_inv = rd(fp);
    Apert = rd(fp); cdt = rd(fp); ckcut = rd(fp); dshell = rd(fp);
    for (int j = 0; j < NSP; j++) {
        Spec *s = &S[j];
        s->par = ri(fp); s->M = rd(fp); s->f = rd(fp); s->Grest = rd(fp); s->h2 = rd(fp);
        s->kres_m = rd(fp); s->Gpert = rd(fp);
        s->nq = ri(fp);
        int nq = s->nq;
        s->q = malloc(nq * sizeof(double)); s->w = malloc(nq * sizeof(double));
        s->n = malloc(2 * nq * sizeof(double)); s->frozen = calloc(2 * nq, sizeof(int));
        s->r = malloc(6 * nq * sizeof(double));
        for (int i = 0; i < nq; i++) { s->q[i] = rd(fp); s->w[i] = rd(fp); }
        for (int i = 0; i < 2 * nq; i++) s->n[i] = rd(fp);
        int hasr = ri(fp);
        if (hasr) for (int i = 0; i < 6 * nq; i++) s->r[i] = rd(fp);
        else s->r[0] = NAN;
        for (int i = 0; i < 2 * nq; i++) s->frozen[i] = ri(fp);
        s->nA = ri(fp); s->nx = ri(fp);
        s->lA = malloc((s->nA > 0 ? s->nA : 1) * sizeof(double));
        s->pp = malloc((s->nA * s->nx > 0 ? s->nA * s->nx : 1) * sizeof(double));
        for (int i = 0; i < s->nA; i++) s->lA[i] = rd(fp);
        for (int i = 0; i < s->nA * s->nx; i++) s->pp[i] = rd(fp);
        s->Pprod_acc = s->Ecoll_acc = 0; s->Gband = s->Gall = 0; s->nband = 0;
    }
    fclose(fp);
    double rho = rhoP[0] + rhoP[1] + rho_h + rho_th;
    H = sqrt(rho / (3 * mP * mP));
}

static double Aof(Spec *s) {           /* mu5,max / m of the parent oscillation */
    double m = mpar[s->par], Sig = sqrt(2.0 * (rhoP[s->par] > 0 ? rhoP[s->par] : 0)) / m;
    return Sig / (2.0 * s->f);
}

static double pdot_of(Spec *s) {       /* dp/dt = -dmu5/dt = A m^2 cos(m (t - tph)) at the current time */
    double m = mpar[s->par];
    return Aof(s) * m * m * cos(m * (t - tph_g[s->par]));
}

/* free-basis axis for helicity sign hs at physical momentum k: bhat = (M, 0, hs k)/E */
static inline void free_axis(double M, double hk, double *bx, double *bz) {
    double E = sqrt(M * M + hk * hk); *bx = M / E; *bz = hk / E;
}

/* first-order superadiabatic vacuum direction at a turning point (mu5 = 0): r_vac = -bhat + delta, delta =
 * (bhat x dbhat/dt)/(2|b|) = (0, -pdot M/(2 E^3), 0), normalised.  Projecting on it removes the reversible boundary
 * "dressing" that a plain free-basis projection at a point with mudot5 != 0 would count as particles. */
static inline void vac_axis(double M, double hk, double pdot, double *vx, double *vy, double *vz) {
    double E = sqrt(M * M + hk * hk);
    double d = -pdot * M / (2.0 * E * E * E);
    /* modes with |delta| > 0.3 are mid-crossing at the turning point (small k): the correction is tapered off there
       (continuous at 0.3, ~1/delta beyond); this concerns k < ~sqrt(A) m, < 1% of the band phase space at A > 10 */
    if (fabs(d) > 0.3) d = (d > 0 ? 0.09 : -0.09) / fabs(d);
    double x = -M / E, y = d, z = -hk / E, nrm = sqrt(x * x + y * y + z * z);
    *vx = x / nrm; *vy = y / nrm; *vz = z / nrm;
}
static double pdot_of(Spec *s);

/* ------------------------------------------------------------------ energies and rates */
static void species_energy(Spec *s) {
    double E = 0, N = 0, ia3 = 1.0 / (a * a * a);
    for (int h = 0; h < 2; h++)
        for (int i = 0; i < s->nq; i++) {
            double k = s->q[i] / a, om = sqrt(k * k + s->M * s->M), nn = s->n[h * s->nq + i];
            E += s->w[i] * om * nn; N += s->w[i] * nn;
        }
    s->E = E * ia3; s->Nn = N * ia3;
}

static void rates(Spec *s, double k, double Tc2, double Tth, double *Gtot, double *Geq, double *feq) {
    double om = sqrt(k * k + s->M * s->M);
    double Gvac = s->Grest * (s->M / om) * Phi_kin(s->M, Tc2);
    double hard = (rho_h > 1e-300 && n_h > 0) ? (12.0 / gsum_f) * n_h * n_h / rho_h : 0.0;
    double Gth = cBB * s->h2 * Tth * Tth / om, Gh = cBB * s->h2 * hard / om;
    double wth = (rho_th + rho_h > 0) ? rho_th / (rho_th + rho_h) : 0.0;
    *Gtot = Gvac + Gth + Gh;
    *Geq = Gvac * (wthermal_inv > 0 ? wth : 0.0) + Gth;
    *feq = (Tth > 0) ? 1.0 / (exp(fmin(om / Tth, 700.0)) + 1.0) : 0.0;
}

/* collisions over dtc for species s (both modes: acts on n; in coherent mode also on the Bloch vectors) */
static void collide(Spec *s, double dtc) {
    double Tc2 = Tcal2(), Tth = Tof(rho_th), ia3 = 1.0 / (a * a * a);
    double A = Aof(s), kband = A * mpar[s->par];
    double gb = 0, nb = 0, ga = 0, na = 0, nbw = 0, wb = 0;
    for (int h = 0; h < 2; h++) {
        double hs = h ? 1.0 : -1.0;
        for (int i = 0; i < s->nq; i++) {
            int id = h * s->nq + i;
            double k = s->q[i] / a, om = sqrt(k * k + s->M * s->M), G, Geq, feq;
            rates(s, k, Tc2, Tth, &G, &Geq, &feq);
            double nn = s->n[id];
            if (k < kband) { gb += s->w[i] * nn * G; nb += s->w[i] * nn; nbw += s->w[i] * nn; wb += s->w[i]; }
            ga += s->w[i] * nn * G; na += s->w[i] * nn;
            if (nocoll || G <= 0) continue;
            double e = exp(-G * dtc), feff = Geq * feq / G;
            double loss = nn * (1 - e), gain = feff * (1 - e);
            double n2 = nn - loss + gain;
            inj_E += s->w[i] * om * loss * ia3; inj_N += 2.0 * s->w[i] * loss * ia3;
            inv_E += s->w[i] * om * gain * ia3;
            s->Ecoll_acc += s->w[i] * om * (loss - gain);
            Ecoll_tot += s->w[i] * om * (loss - gain);
            s->n[id] = n2;
            if (coherent && !s->frozen[id]) {
                double vx, vy, vz; vac_axis(s->M, hs * k, pdot_of(s), &vx, &vy, &vz);
                double *r = &s->r[3 * id];
                double par = -(r[0] * vx + r[1] * vy + r[2] * vz);          /* = 2n - 1 */
                double px = r[0] + par * vx, py = r[1] + par * vy, pz = r[2] + par * vz, d = exp(-0.5 * G * dtc);
                double npar = 2 * n2 - 1;
                r[0] = -npar * vx + d * px; r[1] = -npar * vy + d * py; r[2] = -npar * vz + d * pz;
            }
        }
    }
    s->Gband = nb > 0 ? gb / nb : 0; s->Gall = na > 0 ? ga / na : 0;
    s->nband = wb > 0 ? nbw / wb / 2.0 : 0;
}

/* ------------------------------------------------------------------ SM update over dt (injections applied first) */
static void sm_update(double dt, double aold, double anew) {
    double r = aold / anew;                      /* redshift first, then inject (injections are at anew) */
    rho_h *= r * r * r * r; rho_th *= r * r * r * r; n_h *= r * r * r;
    rho_h += inj_E; n_h += inj_N;
    if (inv_E > 0) {
        if (inv_E <= rho_th) rho_th -= inv_E;
        else { double rest = inv_E - rho_th; rho_th = 0; rho_h -= rest; if (rho_h < 0) rho_h = 0; }
    }
    inj_E = inj_N = inv_E = 0;
    double tau = tau_KL();
    if (dt > 0 && tau < 1e299) {
        double e = exp(-dt / tau);
        double dE = rho_h * (1 - e);
        rho_h -= dE; n_h *= e; rho_th += dE;
    }
}

static double rho_tot(void) {
    double r = rhoP[0] + rhoP[1] + rho_h + rho_th;
    for (int j = 0; j < NSP; j++) r += S[j].E;
    return r;
}

/* ------------------------------------------------------------------ record */
static void record(void) {
    double rho = rho_tot(), EN = 0;
    for (int j = 0; j < NSP; j++) EN += S[j].E;
    fprintf(frec, "%.10e %.10e %.10e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e", t, a, H, rhoP[0] / rho, rhoP[1] / rho,
            EN / rho, rho_h / rho, rho_th / rho, Tof(rho_th), Tcal2(), Tof(rho_h + rho_th), tau_KL());
    for (int j = 0; j < NSP; j++)
        fprintf(frec, " %.6e %.6e %.6e %.6e %.6e %.6e", S[j].E / rho, S[j].Gband, S[j].Gall, S[j].nband, Aof(&S[j]),
                S[j].Nn);
    fprintf(frec, " %.6e %.6e %.10e\n", n_h, rho, rho * a * a * a * a);
    fflush(frec);
}

/* ------------------------------------------------------------------ coherent mode */
static void init_bloch(Spec *s) {
    for (int h = 0; h < 2; h++) {
        double hs = h ? 1.0 : -1.0;
        for (int i = 0; i < s->nq; i++) {
            int id = h * s->nq + i;
            double vx, vy, vz; vac_axis(s->M, hs * s->q[i] / a, pdot_of(s), &vx, &vy, &vz);
            double npar = 2 * s->n[id] - 1;
            s->r[3 * id] = -npar * vx; s->r[3 * id + 1] = -npar * vy; s->r[3 * id + 2] = -npar * vz;
        }
    }
}

/* advance active modes of species s from t1 to t2: a(t) = a1 exp(H1 (t - t1)), mu5 = -A m sin(m (t - tTP)) */
static void advance(Spec *s, double t1, double t2, double a1, double H1, double tph) {
    double m = mpar[s->par], A = Aof(s), M = s->M;
    double kmax = 0;
    for (int i = 0; i < s->nq; i++) if (!s->frozen[i] || !s->frozen[s->nq + i]) kmax = fmax(kmax, s->q[i] / a1);
    double wmax = kmax + A * m + M;
    int nst = (int)ceil((t2 - t1) * wmax / cdt);
    if (nst < 1) nst = 1;
    double dt = (t2 - t1) / nst;
    for (int st = 0; st < nst; st++) {
        double tm = t1 + (st + 0.5) * dt;
        double am = noexp ? a1 : a1 * exp(H1 * (tm - t1));
        double mu5 = -A * m * sin(m * (tm - tph));
        for (int h = 0; h < 2; h++) {
            double hs = h ? 1.0 : -1.0;
            double *r = &s->r[3 * h * s->nq];
            int *fz = &s->frozen[h * s->nq];
            for (int i = 0; i < s->nq; i++) {
                if (fz[i]) continue;
                double p = hs * s->q[i] / am - mu5;
                double b = sqrt(M * M + p * p), ux = M / b, uz = p / b;
                double th = 2.0 * b * dt, c = cos(th), sn = sin(th);
                double *v = &r[3 * i];
                double dot = ux * v[0] + uz * v[2];
                /* u x v with u = (ux, 0, uz): (0*v2 - uz*v1, uz*v0 - ux*v2, ux*v1 - 0*v0) */
                double cx = -uz * v[1], cy = uz * v[0] - ux * v[2], cz = ux * v[1];
                double nx = v[0] * c + cx * sn + ux * dot * (1 - c);
                double ny = v[1] * c + cy * sn;
                double nz = v[2] * c + cz * sn + uz * dot * (1 - c);
                v[0] = nx; v[1] = ny; v[2] = nz;
            }
        }
    }
}

/* at a turning point of the parent: free-basis occupations, production energy, freezing */
static double tp_update(Spec *s) {
    double m = mpar[s->par], A = Aof(s), kcut = (ckcut * A + 4.0 * s->M / m + 2.0) * m, dE = 0;
    for (int h = 0; h < 2; h++) {
        double hs = h ? 1.0 : -1.0;
        for (int i = 0; i < s->nq; i++) {
            int id = h * s->nq + i;
            if (s->frozen[id]) continue;
            double k = s->q[i] / a, om = sqrt(k * k + s->M * s->M), vx, vy, vz;
            vac_axis(s->M, hs * k, pdot_of(s), &vx, &vy, &vz);
            double *r = &s->r[3 * id];
            double nn = 0.5 * (1.0 - (r[0] * vx + r[1] * vy + r[2] * vz));
            dE += s->w[i] * om * (nn - s->n[id]);
            s->n[id] = nn;
            if (k > kcut) { s->frozen[id] = 1; }
        }
    }
    return dE / (a * a * a);          /* physical energy density produced */
}

static void run_coherent(void) {
    double tnext[2], tlastc[2] = {t, t};
    for (int P = 0; P < 2; P++) {
        if (tnext_in[P] < 0) { tph_g[P] = t; tnext[P] = t + M_PI / mpar[P]; tlastc[P] = t; }   /* fresh: t is a TP */
        else { tnext[P] = tnext_in[P]; tlastc[P] = tlastc_in[P]; }                            /* restart */
    }
    for (int j = 0; j < NSP; j++) { if (isnan(S[j].r[0])) init_bloch(&S[j]); species_energy(&S[j]); }
    double rho = rho_tot(); H = sqrt(rho / (3 * mP * mP));
    double trec = t; int nrec = 0;
    record();
    double tlast_sm = t;
    while (t < t_end - 1e-15) {
        int P = tnext[0] < tnext[1] ? 0 : 1;
        double t2 = fmin(tnext[P], t_end);
        double a1 = a, H1 = H;
        for (int j = 0; j < NSP; j++) {
            int pj = S[j].par;
            /* each species advanced with its own parent phase reference */
            advance(&S[j], t, t2, a1, H1, tph_g[pj]);
        }
        double a2 = noexp ? a1 : a1 * exp(H1 * (t2 - t));
        /* condensates redshift between events */
        if (!noexp) { double r = a1 / a2; rhoP[0] *= r * r * r; rhoP[1] *= r * r * r; }
        a = a2; t = t2;
        if (t2 >= tnext[P] - 1e-15) {
            /* turning point of parent P */
            double dEprod = 0;
            for (int j = 0; j < NSP; j++) if (S[j].par == P) dEprod += tp_update(&S[j]);
            if (!fixamp) { rhoP[P] -= dEprod; if (rhoP[P] < 0) rhoP[P] = 0; }
            Eprod_tot += dEprod * a * a * a;
            for (int j = 0; j < NSP; j++) if (S[j].par == P) {
                collide(&S[j], t - tlastc[P]);
                species_energy(&S[j]);
            }
            tlastc[P] = t;
            tnext[P] += M_PI / mpar[P];
            sm_update(t - tlast_sm, a1, a);   /* redshift, injections, KL relaxation over the interval */
            tlast_sm = t;
            for (int j = 0; j < NSP; j++) if (S[j].par != P) {
                /* other species' energies redshift only through omega(k): re-evaluate (occupations unchanged) */
                species_energy(&S[j]);
            }
        } else {
            sm_update(t - tlast_sm, a1, a); tlast_sm = t;
        }
        rho = rho_tot(); H = noexp ? H : sqrt(rho / (3 * mP * mP));
        tnext_in[0] = tnext[0]; tnext_in[1] = tnext[1]; tlastc_in[0] = tlastc[0]; tlastc_in[1] = tlastc[1];
        nrec++;
        if (t >= trec * 1.01 || t >= t_end - 1e-15 || nrec % 400 == 0) { record(); trec = t; }
    }
    record();
}

/* ------------------------------------------------------------------ dephased mode */
static double pp_lookup(Spec *s, double kphys, double A) {
    if (s->nA <= 0) return 0;
    double m = mpar[s->par], lA = log(A);
    if (lA < s->lA[0]) return 0;
    int iA = 0;
    if (lA >= s->lA[s->nA - 1]) { iA = s->nA - 2; lA = s->lA[s->nA - 1]; }
    else { while (iA < s->nA - 2 && s->lA[iA + 1] <= lA) iA++; }
    double u = (lA - s->lA[iA]) / (s->lA[iA + 1] - s->lA[iA]);
    double val = 0;
    /* both bracketing tables are read at the same relative position x = k / ks(A), ks at the actual A, so that the
       band edge (~1.1 A m) moves continuously with A (interpolating at fixed k would smear the wider band of the
       upper table into the lower amplitude) */
    double ks = (1.25 * A + 3.0 * s->M / m + 2.0) * m;
    double x = kphys / ks;
    for (int d = 0; d < 2; d++) {
        if (x >= 1.0 || x < 0) continue;
        double fx = x * (s->nx - 1);
        int ix = (int)fx; double ux = fx - ix;
        if (ix >= s->nx - 1) { ix = s->nx - 2; ux = 1.0; }
        double v = s->pp[(iA + d) * s->nx + ix] * (1 - ux) + s->pp[(iA + d) * s->nx + ix + 1] * ux;
        val += (d ? u : 1 - u) * v;
    }
    return val > 0 ? val : 0;
}

static void run_dephased(void) {
    for (int j = 0; j < NSP; j++) species_energy(&S[j]);
    double rho = rho_tot(); H = sqrt(rho / (3 * mP * mP));
    record();
    double dt = 1e-4, trec = t;
    long it = 0;
    double *R[NSP];
    for (int j = 0; j < NSP; j++) R[j] = malloc(S[j].nq * sizeof(double));
    while (t < t_end) {
        /* production rates per mode at the current state */
        double Pc[2] = {0, 0};
        for (int j = 0; j < NSP; j++) {
            Spec *s = &S[j]; double m = mpar[s->par], A = Aof(s);
            int res = (s->kres_m > 0 && A < Apert);
            for (int i = 0; i < s->nq; i++) R[j][i] = 0;
            if (rhoP[s->par] <= 0) continue;
            if (res) {
                double kr = s->kres_m * m, qr = kr * a;
                double Ntot = s->Gpert * rhoP[s->par] / m;           /* decays per volume per time; each fills the pair (k,h),(-k,h): R_i = g_i Ntot / (w_i/a^3) per entry (resumed run: a factor 1/2 removed) */
                if (dshell == 0) {
                    int ib = -1;
                    for (int i = 0; i < s->nq - 1; i++) if (s->q[i] <= qr && qr < s->q[i + 1]) { ib = (qr - s->q[i] < s->q[i + 1] - qr) ? i : i + 1; break; }
                    if (ib >= 0) R[j][ib] = Ntot / (s->w[ib] / (a * a * a));
                } else {
                    /* products spread as a Gaussian in ln k of width dshell (fragmented / moving parents); each bin
                       receives the share g_i of the decays: R_i = g_i Ntot / (2 w_i / a^3) */
                    double gs = 0, hw = fabs(dshell);
                    if (dshell <= -1.0e5) {
                        /* omega-Lorentzian variant (resumed run): the pair-energy line 2 omega - m of HWHM fac G(k_res)
                           written in omega (HWHM G/2), weight per bin L(x) d omega, d omega = (k^2/omega) d ln k, kept
                           over the kinematic range of chi -> N + N*(-> L H): M < omega < (m^2 + M^2 - (m_L + m_H)^2)/(2m)
                           (k = 0 to s* = (m_L + m_H)^2), no cut in units of the width; fac = -dshell/1e5 */
                        double Gk, Geqk, feqk, Tc2 = Tcal2(); rates(s, kr, Tc2, Tof(rho_th), &Gk, &Geqk, &feqk);
                        double hwo = 0.5 * (-dshell / 1.0e5) * Gk; if (hwo < 1e-9 * m) hwo = 1e-9 * m;
                        double smin = pow(sqrt(kL * Tc2) + sqrt(kH * Tc2), 2.0);
                        double ommax = (m * m + s->M * s->M - smin) / (2.0 * m);
                        for (int i = 0; i < s->nq; i++) {
                            double k = s->q[i] / a, om = sqrt(k * k + s->M * s->M);
                            if (om > ommax) continue;
                            double dl = (i + 1 < s->nq) ? log(s->q[i + 1] / s->q[i]) : log(s->q[i] / s->q[i - 1]);
                            double x = (om - 0.5 * m) / hwo;
                            R[j][i] = (k * k / om) * dl / (1.0 + x * x); gs += R[j][i];
                        }
                        if (gs > 0) for (int i = 0; i < s->nq; i++) if (R[j][i] > 0) R[j][i] = (R[j][i] / gs) * Ntot / (s->w[i] / (a * a * a));
                        goto shell_done;
                    }
                    if (dshell <= -100.0) {
                        /* self-consistent natural line (added in the resumed run): both final-state N have removal
                           rate G(k_res); the pair-energy Lorentzian has HWHM G, i.e. G/2 in omega, and in ln k
                           HWHM = (G/2) omega / k^2 = G / (m beta^2), beta^2 = 4 kres_m^2; times fac = -dshell/1000 */
                        double Gk, Geqk, feqk; rates(s, kr, Tcal2(), Tof(rho_th), &Gk, &Geqk, &feqk);
                        hw = (-dshell / 1000.0) * Gk / (m * 4.0 * s->kres_m * s->kres_m);
                        if (hw < 1e-6) hw = 1e-6;
                    }
                    for (int i = 0; i < s->nq; i++) {
                        /* dshell > 0: Gaussian of width dshell; dshell < 0: Lorentzian of HWHM |dshell| (natural line of
                           damped final states), cut at |x| < LCUT */
                        double x = log(s->q[i] / qr) / hw;
                        if (fabs(x) < (dshell > 0 ? 5 : LCUT)) { double dl = (i + 1 < s->nq) ? log(s->q[i + 1] / s->q[i]) : log(s->q[i] / s->q[i - 1]);
                            R[j][i] = (dshell > 0 ? exp(-0.5 * x * x) : 1.0 / (1.0 + x * x)) * dl; gs += R[j][i]; }
                    }
                    for (int i = 0; i < s->nq; i++) if (R[j][i] > 0) R[j][i] = (R[j][i] / gs) * Ntot / (s->w[i] / (a * a * a));
                }
                shell_done: ;
            } else {
                for (int i = 0; i < s->nq; i++) R[j][i] = (m / (2 * M_PI)) * pp_lookup(s, s->q[i] / a, A);
            }
            for (int h = 0; h < 2; h++)
                for (int i = 0; i < s->nq; i++) {
                    double k = s->q[i] / a, om = sqrt(k * k + s->M * s->M);
                    Pc[s->par] += s->w[i] * om * R[j][i] * (1 - 2 * s->n[h * s->nq + i]) / (a * a * a);
                }
        }
        /* time step */
        double lim = 0.02 / H;
        for (int P = 0; P < 2; P++) if (Pc[P] > 0 && rhoP[P] > 1e-30 * rho) lim = fmin(lim, 0.01 * rhoP[P] / Pc[P]);
        double tk = tau_KL(); if (tk < 1e299 && rho_h > 1e-6 * rho) lim = fmin(lim, 0.2 * tk);
        double dtn = fmin(lim, dt * 1.2);
        if (dtn < 1e-7) dtn = 1e-7;
        if (t + dtn > t_end) dtn = t_end - t;
        dt = dtn;
        /* production (exact for dn/dt = R (1 - 2n)) */
        double dE[2] = {0, 0};
        for (int j = 0; j < NSP; j++) {
            Spec *s = &S[j];
            for (int h = 0; h < 2; h++)
                for (int i = 0; i < s->nq; i++) {
                    if (R[j][i] <= 0) continue;
                    int id = h * s->nq + i;
                    double k = s->q[i] / a, om = sqrt(k * k + s->M * s->M);
                    double n1 = 0.5 + (s->n[id] - 0.5) * exp(-2 * R[j][i] * dt);
                    dE[s->par] += s->w[i] * om * (n1 - s->n[id]) / (a * a * a);
                    s->n[id] = n1;
                }
        }
        for (int P = 0; P < 2; P++) { if (!fixamp) { rhoP[P] -= dE[P]; if (rhoP[P] < 0) rhoP[P] = 0; } Eprod_tot += dE[P] * a * a * a; }
        /* collisions */
        for (int j = 0; j < NSP; j++) collide(&S[j], dt);
        /* expansion */
        double aold = a, anew = noexp ? a : a * exp(H * dt);
        if (!noexp) { double r = aold / anew; rhoP[0] *= r * r * r; rhoP[1] *= r * r * r; }
        a = anew; t += dt;
        sm_update(dt, aold, anew);
        for (int j = 0; j < NSP; j++) species_energy(&S[j]);
        rho = rho_tot(); H = noexp ? H : sqrt(rho / (3 * mP * mP));
        it++;
        if (t >= trec * 1.02 || t >= t_end || it % 2000 == 0) { record(); trec = t; }
        /* stop when the condensates and N are gone */
        double EN = 0; for (int j = 0; j < NSP; j++) EN += S[j].E;
        if ((rhoP[0] + rhoP[1]) < 1e-9 * rho && EN < 1e-9 * rho && t > 10) break;
    }
    record();
}

/* ------------------------------------------------------------------ dump final state */
static void dump(const char *pre) {
    char fn[1024]; snprintf(fn, sizeof fn, "%s_state.txt", pre);
    FILE *fp = fopen(fn, "w");
    fprintf(fp, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", t, a, H, rhoP[0], rhoP[1], rho_h, n_h,
            rho_th, Ecoll_tot, Eprod_tot, rho_tot());
    for (int j = 0; j < NSP; j++) {
        Spec *s = &S[j];
        fprintf(fp, "%d\n", s->nq);
        for (int i = 0; i < 2 * s->nq; i++) fprintf(fp, "%.10e %d\n", s->n[i], s->frozen[i]);
    }
    fclose(fp);
    if (coherent) {                    /* restart file of the coherent stage: event times and Bloch vectors */
        snprintf(fn, sizeof fn, "%s_bloch.txt", pre);
        fp = fopen(fn, "w");
        fprintf(fp, "%.17g %.17g %.17g %.17g %.17g %.17g\n", tph_g[0], tph_g[1], tnext_in[0], tnext_in[1], tlastc_in[0], tlastc_in[1]);
        for (int j = 0; j < NSP; j++) for (int i = 0; i < 6 * S[j].nq; i++) fprintf(fp, "%.15e\n", S[j].r[i]);
        fclose(fp);
    }
}

/* ------------------------------------------------------------------ table mode */
static void table(const char *out, double Mh, int nA, double *As) {
    FILE *fp = fopen(out, "w");
    for (int ia = 0; ia < nA; ia++) {
        double A = As[ia];
        double kmax = 1.25 * A + 3 * Mh + 2;
        int nk = (int)fmax(800, 40 * kmax);
        if (nk > 6000) nk = 6000;
        fprintf(fp, "A %.10g %d\n", A, nk);
        double *pm = malloc(nk * sizeof(double)), *ppl = malloc(nk * sizeof(double));
        for (int h = 0; h < 2; h++) {
            double hs = h ? 1.0 : -1.0;
            for (int i = 0; i < nk; i++) {
                double k = 1e-3 + (kmax - 1e-3) * i / (nk - 1.0);
                double vx, vy, vz; vac_axis(Mh, hs * k, +A, &vx, &vy, &vz);       /* t = 0: pdot = +A */
                double v[3] = {vx, vy, vz};
                double wmax = k + A + Mh;
                int nst = (int)ceil(M_PI * wmax / 0.02); if (nst < 4000) nst = 4000;
                double dt = M_PI / nst;
                for (int st = 0; st < nst; st++) {
                    double tm = (st + 0.5) * dt;
                    double p = hs * k + A * sin(tm);
                    double b = sqrt(Mh * Mh + p * p), ux = Mh / b, uz = p / b;
                    double th = 2 * b * dt, c = cos(th), sn = sin(th);
                    double dot = ux * v[0] + uz * v[2];
                    double cx = -uz * v[1], cy = uz * v[0] - ux * v[2], cz = ux * v[1];
                    double nx = v[0] * c + cx * sn + ux * dot * (1 - c), ny = v[1] * c + cy * sn,
                           nz = v[2] * c + cz * sn + uz * dot * (1 - c);
                    v[0] = nx; v[1] = ny; v[2] = nz;
                }
                vac_axis(Mh, hs * k, -A, &vx, &vy, &vz);                           /* t = pi: pdot = -A */
                double nn = 0.5 * (1 - (v[0] * vx + v[1] * vy + v[2] * vz));
                if (h == 0) pm[i] = nn; else ppl[i] = nn;
            }
        }
        for (int i = 0; i < nk; i++) fprintf(fp, "%.8e %.8e %.8e\n", 1e-3 + (kmax - 1e-3) * i / (nk - 1.0), pm[i], ppl[i]);
        free(pm); free(ppl);
    }
    fclose(fp);
}

int main(int argc, char **argv) {
    if (argc < 3) die("usage");
    if (!strcmp(argv[1], "table")) {
        double Mh = atof(argv[3]); int nA = argc - 4; double *As = malloc(nA * sizeof(double));
        for (int i = 0; i < nA; i++) As[i] = atof(argv[4 + i]);
        table(argv[2], Mh, nA, As);
        return 0;
    }
    read_input(argv[2]);
    char fn[1024]; snprintf(fn, sizeof fn, "%s_rec.txt", argv[3]);
    frec = fopen(fn, "w");
    if (coherent) run_coherent(); else run_dephased();
    fclose(frec);
    dump(argv[3]);
    return 0;
}
