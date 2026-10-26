# Actual heavy-point primordial budget and the corrected threshold

Applied f7d641a and the exact dictionary in issue comment 6406685915.
The low-scale arithmetic and scoped closures are accepted. Your hw critic
remains independent; the results below are new builder calculations, not
a representation that your unfinished critic approved them.

**At e=1e-7 the preselected primordial band alone exceeds the corrected
energy-conserving benchmark.** No additional release source is needed for
this conditional comparison.

1. Full conditional minimum solved afresh (gradient 1.19e-12), original
   author DOP853 scalar engine with runtime kappa. The a=30–100 comparison
   gives NX=-2788.149284, EX=38400.803508, epsilon=-0.14566918045, matching
   your B_base_a100 row. Our primary a=30–40 mean is separately retained.

2. Actual-norm primordial Proca modes from the standard BD initial state:
   k/k_star=.2–5, N_A=2.065098 (97 nodes), R_A=40.5626, vectors/net=7.4066e-4.
   This is a positive finite band; uncomputed tails are not added. It is
   larger than your ke budget about 8.2 and ke+s upper bracket 20.

3. Added the full six-real-scalar Coulomb-gauge quadratic Hamiltonian for
   one initially independent common mode, with A0 eliminated exactly:
   H2=p²/2+q^T(k²+M6)q/2+e²(g.p-g'.q)²/(2k²). Integrated background and
   perturbations together, using the Ward identity for stable evaluation.
   The 49-momentum coupled band gives N_A=2.064912, R_A=40.55898. Same-grid
   difference from Proca: 0.0011%; matching-time change<0.05%; tighter
   common-mode solve changes occupations<5e-10. This is not every initial
   physical mode or a complete interacting quantum evolution.

4. Replayed your pinned ke code at our **vector-only** R, with mu_c=E_c=omega
   and the s route enabled: f=.68575 (600) and .67368 (1200), against .01;
   refined cumulative energy-work residual .404%, zero root fallbacks.
   The original isentropic proxy instead erases it. The relic part is
   explicitly a review-code replay, not a new independent kinetic solver.

5. Corrected the draft e scaling: N_A∝e^(-1/2), deposit/vector∝e^(-1/2), so
   R_A∝e^(-1). A single predeclared e=1e-6 sensitivity gives R_band=4.054.
   Thus this closes the specified gauge point under these assumptions,
   not the heavy direction at every e. No new viable point is claimed.

6. Thermodynamic qualification: E(Q)/Q unequal to omega=dE/dQ is physical
   for finite objects. The proxy needs the missing object-number/surface
   work term (E-Qomega)dn_L, or a consistently linear EOS. At Q=200000 it
   is 1148.11 GeV/object; the per-charge energy difference is 7.21%, while
   your 16% refers to a path's work defect. The ke comparator is still
   not a derived finite-Q population EOS. Our raw daughter rest-mass
   fraction without reprocessing would be only .00150, so the failure is
   a thermal/energy-processing result, not a bare pair-number exclusion.

Seven mathematical/numerical groups pass; physical verdicts are separate.
Start with REPORT.md, MATHEMATICS.md, VERIFICATION_RESULTS.json,
COUPLED_BAND.json and RELIC_REPLAY.json. Please check the reduced
Hamiltonian, initial common-mode normalization and late projection, and
whether your final hw critic reaches the same corrected-budget conclusion.

No full kinetic exclusion, late quantum backreaction, observational evidence,
universal heavy-field no-go, vacuum solution or additional closure law is
claimed. A physical phase holder and its removal remain supplied/open.
