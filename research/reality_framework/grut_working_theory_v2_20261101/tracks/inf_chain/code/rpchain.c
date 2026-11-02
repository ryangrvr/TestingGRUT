/* inf_chain track (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): rpchain.c is a copy of
 * tracks/inf_potential/code/rpinf.c, which is itself a copy of tracks/dq_entropy/code/rpdq.c (the M1r release engine) plus
 * MODE 4 (an explicit inflaton potential).  ONE addition here (marked CH1), everything else untouched:
 *   CH1. optional relaxation of the step cap.  The engine caps the step at 0.05/max(1, b/2) so that the inflaton
 *        oscillation (m_I = b H_i) is resolved.  With the environment variable RPCH_CAPT = t_c (> 0) the cap becomes 0.05
 *        for tau > t_c (use only after the inflaton has decayed, rho_I/rho < 1e-12; long X <-> Sigma freeze-out runs).
 *        Unset (default): byte-identical to rpinf.c.
 * --- header of rpinf.c follows ---
 */
/* inf_potential track (Claude, own code; issue #13 inflaton-sector pass, 3 November 2026): rpinf.c is a copy of
 * tracks/dq_entropy/code/rpdq.c (the M1r release engine) with ONE addition, MODE 4: an explicit inflaton potential in
 * place of the mode-0 quadratic oscillator.  MODE 0..3 are untouched (bit-for-bit regression in inf_potential NOTES).
 *   MODE 4: I is the HOLDER field phi_I = phi - phi_min in M_r units (phihat = I/sqrt(eta2) as in mode 0), and
 *     rho_I = I'^2/2 + v(I),  p_I = I'^2/2 - v(I),  I'' + (3H + Gamma) I' + v'(I) + sqrt(eta2) dV_h/dphihat = 0,
 *     R-trace trI = -I'^2 + 4 v,  Rdot = 2(3H + Gamma) I'^2 + 6 v'(I) I'  (same approximation as mode 0: the hidden
 *     back-reaction is dropped from Rdot only), radiation source Gamma I'^2.
 *   v(I) = AMP * shape(phi), phi = I + PHIMIN (M_r units), AMP = V0/eta2 (V0 in M_r^4), from the environment:
 *     RPINF_KIND = 1 (T: tanh^2(phi/P1)), 2 (E: (1 - exp(-phi/P1))^2), 3 (HT: (1 - P3 x - x^P2)^2/FTOP^2, x = phi/P1,
 *                  FTOP from the environment as RPINF_FTOP);
 *     RPINF_AMP, RPINF_P1, RPINF_P2, RPINF_P3, RPINF_FTOP, RPINF_PHIMIN;
 *     F spec: RPINF_I0 (holder field at the slice, M_r units) and RPINF_SIGN (sign of I'0); |I'0| is fixed by the
 *     constraint 3 f H^2 + 3 H fd = rho at H = 1 (the same iteration as mode 0).  If the constraint leaves no kinetic
 *     energy (need < v(I0)) the run stops with an error.
 *   The step cap uses b = m_I/H_i (pass m_I/H_i as the "b" argument in MODE 4; it is used for nothing else).
 * --- original header of rpdq.c follows ---
/* dq_entropy track (Claude, own code; issue #13 decisive pass, 2 November 2026): rpdq.c is a copy of
 * tracks/fx_holder/code/rptri.c with three changes (marked DQ):
 *   DQ1. the forced-Y diagnostic uses c = 2 conj(u) Sigma (the Y equation is M^2 Y + kap X^3 + 2 conj(u) S conj(Y) = 0 at
 *        rest), its time derivative includes udot (u' = G phihat' [+ LOCK R']) and the portal mass; the light/heavy split
 *        rotates by arg(c)/2 at the rate Im(cdot/c)/2.  (fx_holder SKEPTIC N6: rptri used B = 2 u Sigma.)  Diagnostic only.
 *        DQ1': the forced point includes the Y quartic self-consistently, and the light/heavy split uses the second
 *        variation A' = A + 4 ly|yp|^2, c' = c + 2 ly yp^2 (exact zero free-Y number at the prepared minimum).
 *   DQ2. optional phase-blind induced-portal counterterm  V += cport (|u|^2 - mu^2)(|Sigma|^2 + 2|Y|^2)  (fx_holder SKEPTIC
 *        N1; cport = c(M_Y), 1/(4 pi^2) per unit ln Q), with its exact inflaton back-reaction for the trilinear:
 *        d/dphihat (|u|^2 - mu^2) = 2 Re(conj(G) u).  Not allowed together with LOCK (Option 1A).
 *   DQ3. conformal-time bins (weights d eta = dt/a): lines "C k eta0 eta1 t0 t1 <NQ means>", width cbin in eta.
 *        Conformal time eta is accumulated by the trapezoid rule on accepted steps (eta(0) = 0).
 * Usage: rpdq mode b gamma xi eta2 cH lx ly ls kap M h mu lockre lockim gre gim cport end rtol atol nsamp binw cbin SPEC ...
 * --- original header of rptri.c follows ---
 * fx_holder track (Claude, own code; issue #13 working-theory repair pass, 1 November 2026).
 *
 * Copy of our m1_holder/code/rpholder.c (itself a copy of rp_reheating/code/rpreheat.c), with the phase holder made a
 * renormalizable inflaton trilinear.  Changes relative to rpholder.c:
 *   1. the cubic coefficient is u = mu + LOCK R + G phihat, with
 *        LOCK  complex, Option 1A (xi e^{i alpha}/M_Y) R Sigma* Y^2 + h.c.   (kept for the comparison runs; 0 otherwise)
 *        G     complex dimensionless, the trilinear  L ⊃ -(g_I e^{i alpha_I} phi_I Sigma* Y^2 + h.c.),  G = g_I e^{i alpha_I},
 *        phihat = phi_I / H_i = I / sqrt(eta2)   (I is the inflaton in M_r units; mode 0 only).
 *      V_h = ... + 2 Re[u Sigma* Y^2];  dV_h/dphihat = 2 Re[G Sigma* Y^2].
 *   2. the trilinear acts back on the inflaton exactly:  I'' + (3H + Gamma) I' + b^2 I + sqrt(eta2) dV_h/dphihat = 0
 *      (phi_I = M_r I, V_h = H_i^4 v, so d(eta2 v)/dI = sqrt(eta2) dv/dphihat).  The interaction energy is part of V_h, so
 *      it enters rho, the trace equation (4 V_h) and the constraint exactly; no explicit time dependence is left for G.
 *      (The 1A term keeps rpholder.c's treatment: explicit R-dependence, its F contribution dropped.)
 *   3. the hidden-work ledger adds (dV_h/dphihat) phihat' with phihat' = I'/sqrt(eta2) (and keeps dV/dR Rdot for 1A).
 *   4. the F-spec start iterates the constraint for I0 because V_h now depends on I (three passes; converged to 1e-16).
 *   5. bins: time-averaged diagnostics over consecutive bins of width binw (trapezoid on accepted steps; binw must be a
 *      multiple of end/nsamp so that no step straddles a bin edge): lines "B k t0 t1 <NQ means>".
 *   6. extra diagnostic columns: Re u, Im u, phihat, mL2 = M_Y^2 - |2 u Sigma| (light-Y mass^2 at the forced point).
 * Units: H_i = 1, tau = H_i t, fields in H_i, I in M_r, eta2 = (H_i/M_r)^2.
 * Usage: rptri mode b gamma xi eta2 cH lx ly ls kap M h mu lockre lockim gre gim end rtol atol nsamp binw SPEC ...
 *   SPEC = F Xre Xim Yre Yim Sre Sim Vxre Vxim Vyre Vyim Vsre Vsim      (a = 1, H = 1, inflaton from the constraint)
 *   SPEC = S y0 ... y18                                               (full state: lna I J fields vel rho_r Haux W rhoI)
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double complex cplx;
#define ND 19
#define NQ 30
static int MODE;
static int PKIND = 0; static double PAMP = 0, PP1 = 1, PP2 = 6, PP3 = 0, PFTOP = 1, PPHIMIN = 0;   /* MODE 4 (inf_potential) */
static double vI(double q) {
  double f = q + PPHIMIN;
  if (PKIND == 1) { double t = tanh(f / PP1); return PAMP * t * t; }
  if (PKIND == 2) { double e = -expm1(-f / PP1); return PAMP * e * e; }
  if (PKIND == 3) { double x = f / PP1; double fx = 1 - PP3 * x - pow(x, PP2); return PAMP * fx * fx / (PFTOP * PFTOP); }
  return 0.0;
}
static double dvI(double q) {
  double f = q + PPHIMIN;
  if (PKIND == 1) { double y = f / PP1; double c = cosh(y); return PAMP * 2 * tanh(y) / (c * c) / PP1; }
  if (PKIND == 2) { double e = exp(-f / PP1); return PAMP * 2 * (1 - e) * e / PP1; }
  if (PKIND == 3) { double x = f / PP1; double fx = 1 - PP3 * x - pow(x, PP2); double fp = -PP3 - PP2 * pow(x, PP2 - 1);
                    return PAMP * 2 * fx * fp / PP1 / (PFTOP * PFTOP); }
  return 0.0;
}
static double b, gam, xi, eta2, sqe, cH, lx, ly, ls, kap, MM, hh, mu, CPORT;
static cplx LOCK, GTRI;

static void potential(const cplx *z, double R, double ph, double *V, cplx *F, double *dVdph) {
  cplx x = z[0], y = z[1], s = z[2], u = mu + LOCK * R + GTRI * ph;
  double x2 = creal(x * conj(x)), y2 = creal(y * conj(y)), s2 = creal(s * conj(s));
  double du = CPORT * (creal(u * conj(u)) - mu * mu);                /* DQ2 portal counterterm */
  *V = lx * x2 * x2 + MM * MM * y2 + ly * y2 * y2 + ls * (s2 - hh * hh / 2) * (s2 - hh * hh / 2)
       + 2 * creal(kap * conj(y) * x * x * x + u * conj(s) * y * y) + du * (s2 + 2 * y2);
  F[0] = 2 * lx * x2 * x + 3 * kap * y * conj(x) * conj(x);
  F[1] = MM * MM * y + 2 * ly * y2 * y + kap * x * x * x + 2 * conj(u) * s * conj(y) + 2 * du * y;
  F[2] = 2 * ls * (s2 - hh * hh / 2) * s + u * y * y + du * s;
  *dVdph = 2 * creal(GTRI * conj(s) * y * y) + CPORT * 2 * creal(conj(GTRI) * u) * (s2 + 2 * y2);
}

typedef struct { double H, R, rho, f, fd, m2c, rhoI, pI, K, V, w, VRlock, ph, dVdph; cplx z[3], v[3], F[3]; } geo;

static void geometry(const double *y, geo *g) {
  for (int k = 0; k < 3; k++) { g->z[k] = y[3 + 2 * k] + I * y[4 + 2 * k]; g->v[k] = y[9 + 2 * k] + I * y[10 + 2 * k]; }
  g->K = 0; for (int k = 0; k < 3; k++) g->K += creal(g->v[k] * conj(g->v[k]));
  double x2 = creal(g->z[0] * conj(g->z[0]));
  g->f = 1 - 2 * xi * eta2 * x2;
  g->fd = -4 * xi * eta2 * creal(conj(g->z[0]) * g->v[0]);
  g->ph = (MODE == 0 || MODE == 4) ? y[1] / sqe : 0.0;
  double trI;
  if (MODE == 0) { g->rhoI = 0.5 * y[2] * y[2] + 0.5 * b * b * y[1] * y[1]; g->pI = 0.5 * y[2] * y[2] - 0.5 * b * b * y[1] * y[1]; trI = g->rhoI - 3 * g->pI; }
  else if (MODE == 4) { double vv = vI(y[1]); g->rhoI = 0.5 * y[2] * y[2] + vv; g->pI = 0.5 * y[2] * y[2] - vv; trI = g->rhoI - 3 * g->pI; }
  else if (MODE == 3) { g->rhoI = 0; g->pI = 0; trI = 0; }
  else { g->rhoI = y[18]; g->pI = 0; trI = g->rhoI; }
  double vx2 = creal(g->v[0] * conj(g->v[0]));
  double den = g->f + 12 * xi * xi * eta2 * x2;
  double R = trI / den;                                   /* zeroth pass: eta2 -> 0 */
  for (int pass = 0; pass < 2; pass++) {
    potential(g->z, R, g->ph, &g->V, g->F, &g->dVdph);
    R = (trI + eta2 * (4 * g->V - 2 * g->K + 12 * xi * vx2 - 12 * xi * creal(conj(g->z[0]) * g->F[0]))) / den;
  }
  potential(g->z, R, g->ph, &g->V, g->F, &g->dVdph);
  g->R = R;
  g->rho = g->rhoI + eta2 * (g->K + g->V) + y[15];
  g->H = (-3 * g->fd + sqrt(9 * g->fd * g->fd + 12 * g->f * g->rho)) / (6 * g->f);
  g->m2c = (MODE == 2) ? -cH * g->rhoI / 3.0 : xi * g->R;
  g->w = (g->rhoI > 1e-300) ? g->pI / g->rhoI : 0;
  g->VRlock = 2 * creal(LOCK * conj(g->z[2]) * g->z[1] * g->z[1]);
}

static void rhs(const double *y, double *d) {
  geo g; geometry(y, &g);
  double H = g.H;
  d[0] = H;
  if (MODE == 0) { d[1] = y[2]; d[2] = -(3 * H + gam) * y[2] - b * b * y[1] - sqe * g.dVdph; d[18] = 0; d[15] = gam * y[2] * y[2] - 4 * H * y[15]; }
  else if (MODE == 4) { d[1] = y[2]; d[2] = -(3 * H + gam) * y[2] - dvI(y[1]) - sqe * g.dVdph; d[18] = 0; d[15] = gam * y[2] * y[2] - 4 * H * y[15]; }
  else if (MODE == 3) { d[1] = d[2] = d[18] = 0; d[15] = -4 * H * y[15]; }
  else { d[1] = d[2] = 0; d[18] = -(3 * H + gam) * y[18]; d[15] = gam * y[18] - 4 * H * y[15]; }
  cplx acc[3];
  for (int k = 0; k < 3; k++) acc[k] = -3 * H * g.v[k] - g.F[k];
  acc[0] -= g.m2c * g.z[0];
  for (int k = 0; k < 3; k++) { d[3 + 2 * k] = creal(g.v[k]); d[4 + 2 * k] = cimag(g.v[k]); d[9 + 2 * k] = creal(acc[k]); d[10 + 2 * k] = cimag(acc[k]); }
  double vx2 = creal(g.v[0] * conj(g.v[0]));
  double fdd = -4 * xi * eta2 * (vx2 + creal(conj(g.z[0]) * acc[0]));
  double rhoppI = (MODE == 0 || MODE == 4) ? y[2] * y[2] : g.rhoI;
  if (MODE == 3) rhoppI = 0;
  d[16] = -(rhoppI + 2 * eta2 * g.K + 4 * y[15] / 3 + fdd - H * g.fd) / (2 * g.f);
  double Rdot = (MODE == 0) ? 6 * b * b * y[1] * y[2] + 2 * (3 * H + gam) * y[2] * y[2]
               : (MODE == 4) ? 6 * dvI(y[1]) * y[2] + 2 * (3 * H + gam) * y[2] * y[2] : 0.0;
  double phdot = (MODE == 0 || MODE == 4) ? y[2] / sqe : 0.0;
  d[17] = -6 * H * g.K - 2 * g.m2c * creal(conj(g.z[0]) * g.v[0]) + g.VRlock * Rdot + g.dVdph * phdot;
}

static const char *QN[NQ] = {"tau", "a", "H", "R", "NX", "NY", "NS", "Gauss", "EX", "eps", "NYH", "NYL", "Nsig", "absX", "absY",
                             "S_over_vev", "rho_r", "rho_I", "radfrac", "Hcheck", "Wcheck", "w_I", "m2curv", "KV", "Haux", "rho",
                             "ure", "uim", "phihat", "mL2"};
static double KV0 = 0;

static void diag(double t, const double *y, double *q) {
  geo g; geometry(y, &g);
  double a = exp(y[0]), a3 = a * a * a, a4 = a3 * a;
  double n[3]; for (int k = 0; k < 3; k++) n[k] = -2 * cimag(conj(g.z[k]) * g.v[k]) * a3;
  double EX = a4 * (creal(g.v[0] * conj(g.v[0])) + lx * pow(cabs(g.z[0]), 4));
  double Nc = 2 * sqrt(2 * lx) * pow(EX / (3 * lx), 0.75);
  cplx u = mu + LOCK * g.R + GTRI * g.ph;
  /* DQ1: forced point of  A Y + c conj(Y) = -J,  A = M^2 + 2 du,  c = 2 conj(u) Sigma,  J = kap X^3 */
  double phdot = (MODE == 0 || MODE == 4) ? y[2] / sqe : 0.0;
  double Rdot = (MODE == 0) ? 6 * b * b * y[1] * y[2] + 2 * (3 * g.H + gam) * y[2] * y[2]
               : (MODE == 4) ? 6 * dvI(y[1]) * y[2] + 2 * (3 * g.H + gam) * y[2] * y[2] : 0.0;
  cplx ud = GTRI * phdot + LOCK * Rdot;
  double du = CPORT * (creal(u * conj(u)) - mu * mu), dud = CPORT * 2 * creal(conj(u) * ud);
  double A = MM * MM + 2 * du, Ad = 2 * dud;
  cplx B = 2 * conj(u) * g.z[2], Bd = 2 * conj(ud) * g.z[2] + 2 * conj(u) * g.v[2], J = kap * g.z[0] * g.z[0] * g.z[0];
  double bb = cabs(B), D = A * A - bb * bb, NYH = NAN, NYL = NAN;
  if (D > 0 && bb > 1e-20) {
    /* DQ1': self-consistent forced point with the Y quartic, (A + 2 ly |yp|^2) yp + c conj(yp) = -J (fixed-point
       iteration), second variation A' |d|^2 + Re(conj(c') d^2) with A' = A + 4 ly |yp|^2, c' = c + 2 ly yp^2 */
    cplx yp = (-A * J + B * conj(J)) / D;
    for (int it = 0; it < 60; it++) {
      double Ae = A + 2 * ly * creal(yp * conj(yp)), De = Ae * Ae - bb * bb;
      cplx yn = (-Ae * J + B * conj(J)) / De;
      if (cabs(yn - yp) <= 1e-15 * (1 + cabs(yn))) { yp = yn; break; }
      yp = yn;
    }
    double A2 = A + 4 * ly * creal(yp * conj(yp)); cplx C2 = B + 2 * ly * yp * yp;
    double b2 = cabs(C2), D2 = A2 * A2 - b2 * b2;
    cplx Jd = 3 * kap * g.z[0] * g.z[0] * g.v[0] + Bd * conj(yp) + Ad * yp;
    cplx ypd = (-A2 * Jd + C2 * conj(Jd)) / D2;
    cplx C2d = Bd + 4 * ly * yp * ypd;
    double th = cimag(C2d / C2);                          /* d arg(c')/dt */
    cplx rot = cexp(-0.5 * I * carg(C2));
    cplx dy = (g.z[1] - yp) * rot, dv = (g.v[1] - ypd - 0.5 * I * th * (g.z[1] - yp)) * rot;
    if (D2 <= 0) { A2 = A; b2 = bb; }
    double wH = sqrt(A2 + b2), wL = sqrt(fmax(A2 - b2, 1e-300));
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
  q[24] = y[16]; q[25] = g.rho; q[26] = creal(u); q[27] = cimag(u); q[28] = g.ph; q[29] = A - bb;
}

int main(int argc, char **argv) {
  if (argc < 27) { fprintf(stderr, "usage: see header (%d args given)\n", argc); return 1; }
  int k = 1;
  MODE = atoi(argv[k++]); b = atof(argv[k++]);
  if (MODE == 4) {                                              /* inf_potential: explicit potential from the environment */
    const char *e;
    PKIND = (e = getenv("RPINF_KIND")) ? atoi(e) : 0; PAMP = (e = getenv("RPINF_AMP")) ? atof(e) : 0;
    PP1 = (e = getenv("RPINF_P1")) ? atof(e) : 1; PP2 = (e = getenv("RPINF_P2")) ? atof(e) : 6; PP3 = (e = getenv("RPINF_P3")) ? atof(e) : 0;
    PFTOP = (e = getenv("RPINF_FTOP")) ? atof(e) : 1; PPHIMIN = (e = getenv("RPINF_PHIMIN")) ? atof(e) : 0;
    if (PKIND < 1 || PKIND > 3 || PAMP <= 0) { fprintf(stderr, "MODE 4 needs RPINF_KIND in 1..3 and RPINF_AMP > 0\n"); return 1; }
  } gam = atof(argv[k++]); xi = atof(argv[k++]); eta2 = atof(argv[k++]); cH = atof(argv[k++]);
  sqe = sqrt(eta2);
  lx = atof(argv[k++]); ly = atof(argv[k++]); ls = atof(argv[k++]); kap = atof(argv[k++]); MM = atof(argv[k++]); hh = atof(argv[k++]); mu = atof(argv[k++]);
  double lre = atof(argv[k++]), lim = atof(argv[k++]); LOCK = lre + I * lim;
  double gre = atof(argv[k++]), gim = atof(argv[k++]); GTRI = gre + I * gim;
  CPORT = atof(argv[k++]);
  if (CPORT != 0 && cabs(LOCK) != 0) { fprintf(stderr, "cport with LOCK not implemented\n"); return 1; }
  double end = atof(argv[k++]), rtol = atof(argv[k++]), atol = atof(argv[k++]); int ns = atoi(argv[k++]);
  double binw = atof(argv[k++]);
  double cbin = atof(argv[k++]);
  double y[ND] = {0};
  const char *spec = argv[k++];
  if (spec[0] == 'F') {
    if (argc != k + 12) { fprintf(stderr, "F spec needs 12 numbers\n"); return 1; }
    for (int i = 0; i < 6; i++) y[3 + i] = atof(argv[k++]);
    for (int i = 0; i < 6; i++) y[9 + i] = atof(argv[k++]);
    /* initial state: a = 1, H = 1 fixed by the inflaton amplitude: 3 f H^2 + 3 H fd = rho -> rho_I = 3f + 3fd - eta2 (K+V(I)) */
    cplx z[3], v[3]; double V, K = 0, dVd; cplx F[3];
    for (int i = 0; i < 3; i++) { z[i] = y[3 + 2 * i] + I * y[4 + 2 * i]; v[i] = y[9 + 2 * i] + I * y[10 + 2 * i]; K += creal(v[i] * conj(v[i])); }
    double f = 1 - 2 * xi * eta2 * creal(z[0] * conj(z[0]));
    double fd = -4 * xi * eta2 * creal(conj(z[0]) * v[0]);
    double R0 = 12.0;                                           /* de Sitter start; eta2 corrections below 1e-10 */
    double Ival = 0.0;
    if (MODE == 4) Ival = atof(getenv("RPINF_I0") ? getenv("RPINF_I0") : "0");
    double Jval = 0.0;
    for (int it = 0; it < 4; it++) {
      potential(z, R0, Ival / sqe, &V, F, &dVd);
      double need = 3 * f + 3 * fd - eta2 * (K + V);
      if (MODE == 0) Ival = sqrt(2 * need) / b;
      else if (MODE == 4) { double kin = need - vI(Ival); if (kin < 0) { fprintf(stderr, "MODE 4: no kinetic energy left (%g)\n", kin); return 5; }
                            double sg = atof(getenv("RPINF_SIGN") ? getenv("RPINF_SIGN") : "1"); Jval = (sg < 0 ? -1 : 1) * sqrt(2 * kin); }
      else if (MODE == 3) y[15] = need;
      else y[18] = need;
    }
    if (MODE == 0) y[1] = Ival;
    if (MODE == 4) { y[1] = Ival; y[2] = Jval; }
    y[16] = 1.0; y[17] = 0.0;
  } else if (spec[0] == 'S') {
    if (argc != k + ND) { fprintf(stderr, "S spec needs %d numbers\n", ND); return 1; }
    for (int i = 0; i < ND; i++) y[i] = atof(argv[k++]);
    y[17] = 0.0;
  } else { fprintf(stderr, "unknown spec\n"); return 1; }
  { KV0 = 0; geo g; geometry(y, &g); KV0 = g.K + g.V;
    cplx u = mu + LOCK * g.R + GTRI * g.ph;
    printf("I0 %.17g %.17g %.17g %.17g %.17g\n", y[1], g.ph, creal(u), cimag(u), g.R); }
  static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
    a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
    a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
    b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84,
    e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
  double t = 0, h = 1e-4, q[NQ], qo[NQ];
  diag(0, y, q);
  printf("S"); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n");
  long steps = 0, rejects = 0; double maxG = 0, maxH = 0, maxW = 0, minmL2 = 1e300;
  double s80[NQ] = {0}, s90[NQ] = {0}, sa[NQ] = {0}, sb[NQ] = {0}, w80 = 0, w90 = 0, wa = 0, wb = 0;
  int si = 1, bin = 0; double hmax_cap = 0.05 / fmax(1.0, b / 2.0);
  double capT = getenv("RPCH_CAPT") ? atof(getenv("RPCH_CAPT")) : 0.0;     /* CH1 */
  double eta = 0, sc_[NQ] = {0}, wc = 0, c_eta0 = 0, c_t0 = 0; int cb = 0;     /* DQ3 */
  while (t < end * (1 - 1e-14)) {
    double target = end * si / ns; if (target > end) target = end;
    if (capT > 0 && t > capT) hmax_cap = 0.05;                              /* CH1 */
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
      if (q[29] < minmL2) minmL2 = q[29];
      double dt = t - told;
      if (told >= 0.8 * end * (1 - 1e-12)) { w80 += dt; for (int j = 0; j < NQ; j++) s80[j] += 0.5 * dt * (qo[j] + q[j]); }
      if (told >= 0.9 * end * (1 - 1e-12)) { w90 += dt; for (int j = 0; j < NQ; j++) s90[j] += 0.5 * dt * (qo[j] + q[j]); }
      if (qo[1] >= 30 && q[1] <= 40 + 1e-9) { double da = q[1] - qo[1]; wa += da; for (int j = 0; j < NQ; j++) sa[j] += 0.5 * da * (qo[j] + q[j]); }
      {                                                    /* DQ3 conformal-time bins */
        double de = 0.5 * dt * (1.0 / qo[1] + 1.0 / q[1]);
        eta += de;
        if (cbin > 0) {
          wc += de; for (int j = 0; j < NQ; j++) sc_[j] += 0.5 * de * (qo[j] + q[j]);
          if (eta >= (cb + 1) * cbin) {
            printf("C %d %.17g %.17g %.17g %.17g", cb, c_eta0, eta, c_t0, t); for (int j = 0; j < NQ; j++) printf(" %.17g", sc_[j] / wc); printf("\n");
            cb++; wc = 0; c_eta0 = eta; c_t0 = t; for (int j = 0; j < NQ; j++) sc_[j] = 0;
          }
        }
      }
      if (binw > 0) {
        wb += dt; for (int j = 0; j < NQ; j++) sb[j] += 0.5 * dt * (qo[j] + q[j]);
        if (t >= (bin + 1) * binw * (1 - 1e-12)) {
          printf("B %d %.17g %.17g", bin, bin * binw, t); for (int j = 0; j < NQ; j++) printf(" %.17g", sb[j] / wb); printf("\n");
          bin++; wb = 0; for (int j = 0; j < NQ; j++) sb[j] = 0;
        }
      }
      if (fabs(t - target) < 1e-12 * end) {
        printf("S"); for (int j = 0; j < NQ; j++) printf(" %.17g", q[j]); printf("\n"); si++;
      }
      if (hstep >= h * (1 - 1e-12)) h = hstep * fac; else h = fmax(h, hstep * fac);
    } else { rejects++; h = hstep * fac; }
    if (steps > 400000000L) { fprintf(stderr, "step limit\n"); return 4; }
  }
  if (w80 > 0) { printf("W80 %.17g", w80); for (int j = 0; j < NQ; j++) printf(" %.17g", s80[j] / w80); printf("\n"); }
  if (w90 > 0) { printf("W90 %.17g", w90); for (int j = 0; j < NQ; j++) printf(" %.17g", s90[j] / w90); printf("\n"); }
  if (wa > 0) { printf("A3040 %.17g", wa); for (int j = 0; j < NQ; j++) printf(" %.17g", sa[j] / wa); printf("\n"); }
  printf("ETA %.17g\n", eta);
  printf("Y"); for (int i = 0; i < ND; i++) printf(" %.17g", y[i]); printf("\n");
  printf("END %ld %ld %.6g %.6g %.6g %.6g\n", steps, rejects, maxG, maxH, maxW, minmL2);
  (void)QN;
  return 0;
}
