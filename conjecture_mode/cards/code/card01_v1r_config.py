"""CARD #1 v1R — POST-DATA CONVERGENCE VERIFICATION: frozen numerical protocol (owner ruling v1R).

Everything scientific is inherited unchanged from v1 (card01_run_config.py at c9899abc + CR-1..CR-4, b50f5224):
physics, eps < 0 primary branch, eps grid, PRIMARY dataset, parameter ranges, CAMB accuracy, likelihoods/data,
S1-S8 thresholds, the S5 decide() function, the CPL comparator definition, eps = 0 null, chi2_eff definition.
Only the numerical optimization procedure below is new.  Frozen before any v1R optimization is executed.
"""
import os

V1_OUT = os.environ.get("CARD01_V1_OUT", "")        # directory of the frozen v1 raw outputs (read-only use)

# CR-5 (owner-authorized for v1R only)
RHOEND = 0.005                    # BOBYQA convergence tolerance (v1: 0.05)
CYCLE_IMPROVEMENT_STOP = 0.1      # converged when one full refinement cycle improves chi2 by < 0.1
MAX_CYCLES = 6                    # if not converged after this many cycles -> target NOT CONVERGED (no state forced)
COVMAT_LCDM = "v1r_covmat_lcdm.txt"   # from official DESI base/<primary> chains (checksum-verified)
COVMAT_W0WA = "v1r_covmat_w0wa.txt"   # from official DESI base_w_wa/<primary> chains (checksum-verified)

# Independent convergence check (owner ruling v1R s.4): retains the frozen v1 0.2 scale as acceptance scale
V1R_SEED_B = 90017                # arm-B seed base, distinct from v1 seeds 17 / 1017 / 2017
FREE_B_EPS_START = -0.20         # free-eps arm-B start: distinct admissible eps on the far side of the v1 minimum
ACCEPT_AB = 0.2                   # |chi2_best,A - chi2_best,B| <= 0.2 for every load-bearing target
MAX_LB_ITERATIONS = 3             # re-identify load-bearing set after arm-B results (best-of-A,B), at most 3 times

# Arm A order (all 22 targets); arm B only for load-bearing targets (identified from data, mechanically)
TARGETS = (["card:0.0", "cpl", "free", "card:-0.08", "card:-0.1", "card:-0.06", "card:-0.12", "card:-0.05",
            "card:-0.15", "card:-0.04", "card:-0.2", "card:-0.03", "card:-0.02", "card:-0.015", "card:-0.01",
            "card:-0.005", "card:-0.25", "card:-0.3", "card:-0.4", "card:-0.5", "card:-0.7", "card:-1.0"])
