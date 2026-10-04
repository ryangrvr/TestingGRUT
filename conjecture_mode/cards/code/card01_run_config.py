"""CARD-01 v1 RUN CONFIGURATION — frozen before any likelihood value is inspected.

Committed after the owner threshold ruling (86ccbb2e) and BEFORE the network/provenance preflight.
Implements the frozen spec (8c634145) and the owner ruling S1-S8 / O3. Nothing here may be changed after
data access (C-F6, owner ruling).

Software (pinned): python 3.11.15; camb 2.0.4; cobaya 3.6.2; act_dr6_lenslike 1.2.1 (data v1.2);
Py-BOBYQA 1.5.0; clipy-like 0.15 (benabed/clipy tag clipy_0.15 @ ad1aff3b); numpy 2.4.6; scipy 1.17.1; getdist 1.7.7.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- theory (accuracy frozen)
# ACT DR6 lensing requires lmax 4000 (trim 2998); accuracy follows that likelihood's documented needs.
CAMB_EXTRA_ARGS = {
    "lens_potential_accuracy": 4,
    "lens_output_margin": 1250,
    "AccuracyBoost": 1.0,
    "lSampleBoost": 1.0,
    "lAccuracyBoost": 1.0,
    "halofit_version": "mead2016",
    "num_massive_neutrinos": 1,
    "nnu": 3.044,
    "dark_energy_model": "ppf",
    # CR-2 (pre-evaluation conformance, see CARD_01_RUN_CORRECTIONS.md): DESI DR2 baseline BBN table
    "bbn_predictor": "PArthENoPE_880.2_standard.dat",
}


def theory_card01():
    return {"card01_cobaya_theory.Card01CAMB": {"python_path": HERE, "extra_args": dict(CAMB_EXTRA_ARGS)}}


def theory_camb():
    """Stock CAMB for LambdaCDM-by-cosmological-constant and the CPL comparator (PPF)."""
    return {"camb": {"extra_args": dict(CAMB_EXTRA_ARGS)}}


# ---------------------------------------------------------------- cosmological parameters (shared by all models)
BASE_PARAMS = {
    "logA": {"prior": {"min": 1.61, "max": 3.91}, "ref": 3.05, "proposal": 0.001, "drop": True, "latex": r"\log(10^{10} A_\mathrm{s})"},
    "As": {"value": "lambda logA: 1e-10*np.exp(logA)"},
    "ns": {"prior": {"min": 0.8, "max": 1.2}, "ref": 0.965, "proposal": 0.002},
    "H0": {"prior": {"min": 20, "max": 100}, "ref": 67.5, "proposal": 0.5},
    "ombh2": {"prior": {"min": 0.005, "max": 0.1}, "ref": 0.0224, "proposal": 0.0001},
    "omch2": {"prior": {"min": 0.001, "max": 0.99}, "ref": 0.120, "proposal": 0.0005},
    "tau": {"prior": {"min": 0.01, "max": 0.8}, "ref": 0.055, "proposal": 0.003},
    "mnu": 0.06,
    "rdrag": None,
}

# S7 (owner-approved): primary eps in [-1, 0]; control eps in [0, 1] (CONTROL - NOT A GRUT CLAIM)
EPS_PRIMARY = (-1.0, 0.0)
EPS_CONTROL = (0.0, 1.0)
# Profile grid frozen in the pre-data spec (8c634145, sec. C0.5(b)); control is the mirror image.
EPS_GRID_PRIMARY = [0.0, -0.005, -0.01, -0.015, -0.02, -0.03, -0.04, -0.05, -0.06, -0.08, -0.10, -0.12,
                    -0.15, -0.20, -0.25, -0.30, -0.40, -0.50, -0.70, -1.00]
EPS_GRID_CONTROL = [-e for e in EPS_GRID_PRIMARY]

# CPL comparator (C-F7), identical likelihood; standard DESI-style ranges with w0 + wa < 0.
CPL_PARAMS = {
    "w": {"prior": {"min": -3.0, "max": 1.0}, "ref": -0.9, "proposal": 0.02},
    "wa": {"prior": {"min": -3.0, "max": 2.0}, "ref": -0.3, "proposal": 0.05},
}
CPL_PRIOR_CONSTRAINT = {"w0wa_lt0": "lambda w, wa: 0 if w + wa < 0 else -np.inf"}

# ---------------------------------------------------------------- likelihood combinations (S6, owner-approved)
BAO = {"bao.desi_dr2.desi_bao_all": None}
CMB = {  # DESI DR2 paper baseline CMB
    # CR-1 (pre-evaluation conformance): the DESI DR2 baseline uses the clik versions (official PLA data via clipy)
    "planck_2018_lowl.TT_clik": None,                  # Planck PR3 low-ell TT Commander (commander_dx12_v3_2_29.clik)
    "planck_2018_lowl.EE_clik": None,                  # Planck PR3 low-ell EE SimAll (simall_100x143_offlike5_EE_Aplanck_B.clik)
    "planck_NPIPE_highl_CamSpec.TTTEEE": None,         # Planck PR4 NPIPE CamSpec high-ell TTTEEE
    "act_dr6_lenslike.ACTDR6LensLike": {               # Planck PR4 + ACT DR6 lensing, likelihood v1.2
        "variant": "actplanck_baseline", "lens_only": False, "lmax": 4000, "version": "v1.2"},
}
COMBINATIONS = {
    "PRIMARY_DESI_CMB": {**BAO, **CMB},                       # verdict-bearing (S6)
    "EXT_PANTHEONPLUS": {**BAO, **CMB, "sn.pantheonplus": None},
    "EXT_UNION3": {**BAO, **CMB, "sn.union3": None},
    "EXT_DESY5": {**BAO, **CMB, "sn.desy5": None},
}
PRIMARY = "PRIMARY_DESI_CMB"

# ---------------------------------------------------------------- minimizer protocol (frozen in spec)
# CR-3 (pre-evaluation): ignore_prior False so the Gaussian nuisance priors shipped with the likelihoods are kept;
# flat cosmological priors are removed analytically (chi2_eff in card01_d3_run.py). 3 starts = 3 separate seeded
# runs (CR-4) so start-to-start disagreement can be flagged.
MINIMIZER = {"minimize": {"method": "bobyqa", "best_of": 1, "ignore_prior": False}}
N_STARTS = 3
CHI2_START_DISAGREEMENT_FLAG = 0.2      # flag a grid point if the 3 starts disagree by more than this
CONVERGENCE_NOTE = "profile chi2 = 2*(-log L) minimized over all non-eps parameters incl. likelihood nuisances"

# ---------------------------------------------------------------- owner-locked thresholds (S1-S5; transcribed)
EPS_NULL = 0.03
I_MIN = 4.0            # S3/S5: improvement threshold for I_card and I_CPL
F_MIN = 0.5            # S3: recovery fraction
PROFILE_SET_DCHI2 = 3.84   # descriptive 95% profile set relative to branch best fit (NOT the null calibration)
EPS_RANGE_BOUNDARY = -1.0  # S7: optimum at this boundary => range-sensitivity repair, no final verdict


def decide(I_CPL, I_card, eps_hat, profile_set_min, profile_set_max, eps_hat_at_boundary):
    """S5 decision tree, applied mechanically (owner ruling). MODEL-KILLED is not available in v1."""
    if eps_hat_at_boundary:
        return "NO FINAL VERDICT - RANGE-SENSITIVITY REPAIR (S7)"
    if I_CPL < I_MIN:
        return "INCONCLUSIVE"
    if profile_set_min > -EPS_NULL and profile_set_max <= 0.0:
        return "EXPLANATORY-KILLED"
    if eps_hat < -EPS_NULL and I_card >= I_MIN and I_card / I_CPL >= F_MIN:
        return "EXPLANATORY-SURVIVES"
    if I_card >= I_MIN and I_card / I_CPL < F_MIN:
        return "EXPLANATORY-KILLED"
    return "INCONCLUSIVE"
