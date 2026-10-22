/* smdr_driver.c -- verifier v_vacua (q-theory construction stage 2, issue #13).

   Own driver around the SMDR library of S.P. Martin and D.G. Robertson (arXiv:1907.02500, SMDR v1.3,
   http://www.niu.edu/spmartin/SMDR, GPL).  The library is NOT in the repository; v3_build_smdr.sh downloads and
   builds it in the scratchpad.  SMDR conventions: V = m2 |H|^2 + k |H|^4 (k = lambda in the same normalization as
   Buttazzo et al. and track `vacua`), betaX = Q dX/dQ, kappa = 1/(16 pi^2).

   Modes (all output is "key value" lines on stdout):
     betas  L g3 g gp yt yb ytau k              -> beta_k, beta_yt, beta_g3, beta_g, beta_gp (Q d/dQ) at loop order L
     run    L Q0 g3 g gp yt yb ytau k  N  Qmin Qmax
                                                -> couplings at N log-spaced Q in [Qmin, Qmax] (MSbar inputs at Q0)
     fit    Mt Mh alphaS                         -> SMDR full-precision MSbar couplings at Q = Mt (pole-mass inputs;
                                                    all other inputs = SMDR's PDG-2025 reference values)
     fitrun L Mt Mh alphaS N Qmin Qmax          -> fit, then as 'run' (also prints beta_k at each Q)
*/
#include "smdr.h"
#include <string.h>
#include <stdlib.h>

/* internal SMDR globals/functions (declared in src/smdr_internal.h, common symbols in libsmdr.a) */
extern SMDR_REAL SMDR_beta_k, SMDR_beta_yt, SMDR_beta_g3, SMDR_beta_g, SMDR_beta_gp, SMDR_beta_v;
int SMDR_Betas (int loopOrder);

static void zero_light (void)
{
  SMDR_yc_in = SMDR_ys_in = SMDR_yu_in = SMDR_yd_in = 0;
  SMDR_ymu_in = SMDR_ye_in = 0;
  SMDR_v_in = 246.0; SMDR_m2_in = -8600.0; SMDR_Lambda_in = 0;
  SMDR_Delta_alpha_had_5_MZ_in = 0.02783;
}

static void print_couplings (const char *tag)
{
  printf ("%s Q %.12Le g3 %.12Le g %.12Le gp %.12Le yt %.12Le yb %.12Le ytau %.12Le k %.12Le v %.12Le\n", tag,
          SMDR_Q, SMDR_g3, SMDR_g, SMDR_gp, SMDR_yt, SMDR_yb, SMDR_ytau, SMDR_lambda, SMDR_v);
}

static void scan (int L, int N, SMDR_REAL Qmin, SMDR_REAL Qmax)
{
  int i;
  SMDR_Load_Inputs ();
  for (i = 0; i < N; i++) {
    SMDR_REAL Qi = Qmin * powl (Qmax / Qmin, (SMDR_REAL) i / (SMDR_REAL) (N - 1));
    SMDR_RGrun_SM (Qi, L);
    SMDR_Betas (L);
    printf ("pt Q %.12Le k %.14Le betak %.14Le yt %.12Le g3 %.12Le g %.12Le gp %.12Le v %.12Le\n",
            SMDR_Q, SMDR_lambda, SMDR_beta_k, SMDR_yt, SMDR_g3, SMDR_g, SMDR_gp, SMDR_v);
  }
}

static void do_fit (SMDR_REAL Mt, SMDR_REAL Mh, SMDR_REAL as)
{
  /* initial guess = SMDR reference model (MSbar values at Q = 172.4) */
  SMDR_Read_MSbar_Inputs ("ReferenceModel.dat");
  SMDR_Fit_Inputs (Mt, as, SMDR_alpha_EXPT, SMDR_GFermi_EXPT, SMDR_MZ_EXPT, Mh, Mt,
                   SMDR_mbmb_EXPT, SMDR_mcmc_EXPT, SMDR_ms_2GeV_EXPT, SMDR_md_2GeV_EXPT, SMDR_mu_2GeV_EXPT,
                   SMDR_Mtau_EXPT, SMDR_Mmuon_EXPT, SMDR_Melectron_EXPT, SMDR_Delta_alpha_had_5_MZ_EXPT,
                   1.0e-11);
  SMDR_Load_Inputs ();
  print_couplings ("fit");
  printf ("fitMW MWpole %.10Le MWPDG %.10Le\n", SMDR_MW_pole, SMDR_MW_PDG);
}

int main (int argc, char *argv[])
{
  if (argc < 2) { fprintf (stderr, "usage: see header\n"); return 1; }
  if (!strcmp (argv[1], "betas")) {
    int L = atoi (argv[2]);
    SMDR_g3 = strtold (argv[3], NULL); SMDR_g = strtold (argv[4], NULL); SMDR_gp = strtold (argv[5], NULL);
    SMDR_yt = strtold (argv[6], NULL); SMDR_yb = strtold (argv[7], NULL); SMDR_ytau = strtold (argv[8], NULL);
    SMDR_lambda = strtold (argv[9], NULL);
    SMDR_yc = SMDR_ys = SMDR_yu = SMDR_yd = SMDR_ymu = SMDR_ye = 0;
    SMDR_v = 246.0; SMDR_m2 = -8600.0; SMDR_Lambda = 0; SMDR_Q = 1000.0;
    SMDR_Betas (L);
    printf ("beta_k %.15Le\nbeta_yt %.15Le\nbeta_g3 %.15Le\nbeta_g %.15Le\nbeta_gp %.15Le\nbeta_v %.15Le\n",
            SMDR_beta_k, SMDR_beta_yt, SMDR_beta_g3, SMDR_beta_g, SMDR_beta_gp, SMDR_beta_v);
    return 0;
  }
  if (!strcmp (argv[1], "run")) {
    int L = atoi (argv[2]);
    zero_light ();
    SMDR_Q_in = strtold (argv[3], NULL);
    SMDR_g3_in = strtold (argv[4], NULL); SMDR_g_in = strtold (argv[5], NULL); SMDR_gp_in = strtold (argv[6], NULL);
    SMDR_yt_in = strtold (argv[7], NULL); SMDR_yb_in = strtold (argv[8], NULL);
    SMDR_ytau_in = strtold (argv[9], NULL); SMDR_lambda_in = strtold (argv[10], NULL);
    scan (L, atoi (argv[11]), strtold (argv[12], NULL), strtold (argv[13], NULL));
    return 0;
  }
  if (!strcmp (argv[1], "fit")) {
    do_fit (strtold (argv[2], NULL), strtold (argv[3], NULL), strtold (argv[4], NULL));
    return 0;
  }
  if (!strcmp (argv[1], "fitrun")) {
    int L = atoi (argv[2]);
    do_fit (strtold (argv[3], NULL), strtold (argv[4], NULL), strtold (argv[5], NULL));
    scan (L, atoi (argv[6]), strtold (argv[7], NULL), strtold (argv[8], NULL));
    return 0;
  }
  fprintf (stderr, "unknown mode\n");
  return 1;
}
