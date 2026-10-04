# CARD #1 — PRE-EVALUATION RUN CORRECTIONS (executor; for owner review)

**Timing.** These corrections were committed **after** product acquisition and **before any likelihood was evaluated**. No Card #1, ΛCDM or CPL likelihood or χ² value existed when they were made.

**What did not change.** None of them alters:
- the postulate, projection or primary branch;
- thresholds S1–S8 or the S5 decision logic;
- the owner-approved S6 dataset composition.

**Their purpose.** They bring the executor's implementation into line with the owner's S6 / O3 wording ("the versions used in the DESI DR2 analysis"). This conformance was checked against the **official** DESI DR2 release, which became reachable once the owner enabled full network access:
- `data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.md`;
- the official `cobaya/base_w_wa/desi-bao-all_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing/chain.input.yaml`.

| ID | Frozen config (`c9899abc`) | Official DESI DR2 baseline | Correction | Why it is not a result-driven change |
|---|---|---|---|---|
| **CR-1** | Planck low-ℓ TT/EE: Cobaya **native** Python versions (`planck_2018_lowl.TT`, `.EE`) | `planck2018-lowl-TT-clik`, `planck2018-lowl-EE-clik` (README); `planck_2018_lowl.TT_clik`, `.EE_clik` (input yaml) | Switch to `planck_2018_lowl.TT_clik` / `.EE_clik`, using the official PLA data `COM_Likelihood_Data-baseline_R3.00` (Commander `commander_dx12_v3_2_29.clik`; SimAll `simall_100x143_offlike5_EE_Aplanck_B.clik`) read by clipy-like 0.15 | Owner S6 names "Commander TT + simall EE", and O3 requires the versions used in the DESI DR2 analysis; DESI used clik. Cobaya documents the native EE as equivalent to clik SimAll, but the native TT is a separate Python translation. Made before any evaluation |
| **CR-2** | CAMB default BBN table (CAMB 2.0.4 default `PRIMAT_Yp_DH_ErrorMC_2024.dat`; gives Y_He = 0.24567 at a reference point) | `bbn_predictor: PArthENoPE_880.2_standard.dat` (Y_He = 0.24541 at the same point) | Set `bbn_predictor: PArthENoPE_880.2_standard.dat` | Conformance to the official baseline theory settings. Every other frozen accuracy setting already matched DESI exactly: `lens_potential_accuracy` 4, lens margin 1250 (named `lens_output_margin` in CAMB 2.x), `mead2016`, PPF, the Accuracy/lSample/lAccuracy boosts = 1, `lmax` 4000, one massive neutrino, N_eff 3.044, Σm_ν 0.06 |
| **CR-3** | minimizer `ignore_prior: True` | DESI posterior maximization | `ignore_prior: False`. The flat cosmological priors are removed analytically: `chi2_eff = 2·minuslogpost − 2·Σ ln(width)` over the uniform-prior parameters | `ignore_prior: True` would drop the Gaussian nuisance priors shipped with the Planck likelihoods (calibration / foregrounds), which are part of the standard likelihood treatment. With this correction, χ² is the standard profile quantity, identical in treatment for ε-grid, ΛCDM and CPL runs within each combination |
| **CR-4** | `best_of: 3` with fixed (point) `ref` values, which would give 3 identical starts | DESI input uses normal `ref` distributions | Run 3 separately seeded single-start minimizations per point, with `ref` drawn from normal distributions using DESI's ref widths | Implements the frozen "3 independent starts" protocol and its 0.2 χ² disagreement flag |

**Noted, not changed (pre-existing, frozen):**
- **Sampling coordinate.** The frozen config samples **H0**, while DESI samples `theta_MC_100` (with `theta_H0_range [20, 100]`). For a profile **minimization** this is a reparameterization with equivalent bounds, so it does not change the minimum.
- **CPL prior.** The frozen comparator imposes `w0 + wa < 0`. The official input yaml shows no such constraint, though DESI's DR2 paper text describes it. It is kept frozen and flagged: it can matter only if the CPL best fit had `w0 + wa ≥ 0`.

Under ER-2, the owner may reject any of CR-1 to CR-4. Rejecting a correction means rerunning with the original frozen configuration.
