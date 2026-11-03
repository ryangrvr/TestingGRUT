# Brief and pre-registration: the reheating closure from fermions on the lattice

This file is committed **before** any of the calculations it describes. Claude, for the owner (issue #13). Read
`tracks/pr_principle/NOTES.md`, `tracks/e2e_card/END_TO_END.md` §5–§7, `tracks/e2e_lattice/NOTES.md`,
`tracks/e2e_kinetic/NOTES.md` and `tracks/e2e_chain/NOTES.md` and its SKEPTIC.md first. The rules of
[BRIEF_E2E.md](BRIEF_E2E.md) apply.

## Why this calculation comes first

The M2 windows depend on how Pauli blocking acts during reheating:
- **P0 (local blocking):** K = 8.8–9.2, and M2-Q is robust.
- **GM (box-mean blocking):** K = 15.3–16.4, and M2-Q survives only in a narrow pocket.

This is not a free parameter; Fermi statistics decides it. Estimated scales [E]: the singlet wavelength is about
0.05 H_i⁻¹ and the χ fragments are about 0.03 H_i⁻¹, so the truth should lie between P0 and GM. No principle should be
judged against windows that move with an uncomputed closure.

## Pre-registered questions and hypotheses (answered and reported whatever the outcome)

- **Q1.** Where does the computed closure lie? Give an interpolation parameter λ_PB, defined so that 0 = P0 and 1 = GM.
  Is it unique, or does it depend on the history (time, species, H_i)?
- **Q2.** The new K, H_i and windows for M2-X and M2-Q in (c0, δ), with the branch-aware fixpoint of the e2e_chain
  skeptic. Do they move, and by how much?
- **Q3.** Re-evaluate candidates B (maximal CP torque), D (maximal yield), E (edge of stability) and F (extremal S/Q,
  branch-defining) from `tracks/pr_principle` against the new windows, with the same code and definitions.
- **H1 (softest held mode, stated now as a hypothesis).** The softest held hidden mode during inflation sits at
  c0_soft = 1.08 / 1.13 / 1.17 for δ = −4 / 0 / +4°. That point does not depend on reheating.
  - H1 **succeeds** if c0_soft at the window's own δ lies inside a self-consistent window, M2-X or M2-Q, on the
    computed closure.
  - H1 **fails** otherwise.
  - Chance baseline: about 3% for a log-uniform c0 over two decades.
- **Q4.** Does any principle become testable only once the closure is known? List any such principle **before** the
  windows are computed, in `tracks/lf_fermion/NOTES.md` §0, or state that there is none.

## Tracks

### lf_fermion: the closure from fermions on the lattice

**The fermions.** Evolve the axially coupled singlet fermions −(M/2)N̄ e^{iγ5φ/f} N (mf_hybrid / e2e_sm couplings) on
the inhomogeneous σ(x, t) and χ(x, t) of the e2e_lattice runs (HFa, design H_i). Use a method that keeps the full
Fermi statistics:
- stochastic "low-cost" fermion pairs, or an exact mode-function method on a smaller box;
- Wilson or an equivalent doubler treatment;
- verified literature for the method.

**Regression first.** With a homogeneous background, the method must reproduce mf_hybrid's exact two-level
half-period production and the skeptic's saturation F_sat.

**Non-backreacting test.** On the same recorded background, compare the actual fermion production and occupation per
species with what P0 and GM predict. Extract λ_PB(t), or an effective coarse-graining length if that is the better
description.

**Then the closure.** Implement the calibrated closure in the lattice energy drain and rerun to the end of the
coherent stage.

**Exports.** `results/hist_HFa_LF.json`, and `results/fluid` state moments for e2e_kinetic's handover. Give the
convergence in box size, spacing, realisations and time step.

### lf_fermion skeptic

Attack the method, the regression, λ_PB, the doubler and lattice artefacts, the noise of the stochastic estimator, and
the backreaction.

### lf_chain: the new windows

- Splice the new history with e2e_chain's handover (energy-conserving fixes included), then run the kinetic network.
- Run the (c0, δ) scans and the branch-aware fixpoint for M2-X and M2-Q.
- Report the new points, windows, H_i, r and e_crit.
- Answer Q2, Q3 and H1 with `tracks/pr_principle/code` (unchanged definitions).
