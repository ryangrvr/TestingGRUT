# Review of Perplexity's interacting-matter notes and PR #14

7 November 2026 · Codex review · Zenodo publication remains paused

**Verdict:** retain the Schwinger-model benchmark and its successful amended spectrum test. Correct the continuum-matching claim, the scope of the light-cone constraint, the gravity exclusions, and the repository replay path before incorporating this packet into the integrated v1. The calculations add a useful interacting realization of supplied physics. They do not yet select GRUT's microscopic dynamics or establish a new physical law.

Reviewed [PR #14](https://github.com/ryangrvr/TestingGRUT/pull/14) at commit `8dceba0c5f10b1febc5cc130b7a63da74e2372cc`, based on `e88f3ba24bec3834cef599a204b4cd2796b01083`. Read both probe reports, RESULTS.md, HANDOFF.md, both registrations, the code and raw results. The 35 files in the supplied ZIP are byte-identical to the PR after mapping `qed/` to `ed/` and `qed_mps/` to `mps/`. That rename explains a real repository-only replay failure.

## What was actually rerun

These are fresh executions of supplied solvers, with additional operator and analytic checks. They are not a wholly independent DMRG implementation or a rerun of the entire production campaign.

| Check | Result in this review | Meaning |
|---|---|---|
| Original repository validation command, run from `mps/` | Fails with missing `../qed/spectrum_cache.json` | Reproduction-path defect in the PR; the ZIP layout works |
| All four original DMRG/archived-ED validation cases, using the ZIP | PASS; largest relative gap difference `1.34e-13` | Solver replay agrees with the archived ED references |
| Fresh ED diagonalization and DMRG, N=12, x=9, m/g=0.25 | Gaps agree to `3.67e-14`; ground energies to `1.14e-13` absolute | Also checks a freshly computed reference, not only a cache |
| Direct dense contraction of the MPO versus the ED Hamiltonian, N=4,6,8 | Maximum entry difference zero in the neutral sector | Checks normalization, staggering, hopping and Coulomb channels |
| Unchanged production point N=36, x=9, m/g=0.25 | Gap `1.1205664352507085`; archive `1.1205664352507700` | Absolute difference `6.16e-14` |
| Frozen extrapolation script on all archived production data | Both original ±3% verdicts reproduced; largest numerical output difference `3.78e-15` | Confirms the analysis, not all 32 production points |
| Original histories at N=10 and N=14, both masses | Residuals reproduce within `2.55e-14` | Four of six history cases freshly rerun |
| Gravity, particles and light-cone scripts | Parsed JSON results agree exactly | Arithmetic reproduced; physical qualifications below remain essential |

Environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, TeNPy 1.1.1; one BLAS/OpenMP thread. Raw outputs, logs and hashes are in [raw/](raw/). N=192 production and bond-dimension checks were inspected, not rerun. Neither the largest history lattice nor the full production grid was rerun.

## 1. The spectrum result is sound within the registered test

`analyse.py` implements the amendment: four finite-size fits crossed with three continuum fits, the declared primary combination, subtraction of `2m`, and the unchanged ±3% acceptance rule. Both registration hashes match. Git places the amendment in `dd1e107`, the analysis in `c4e1ddf`, and the new production results in `8dceba0`, in that order. The original ED misses remain recorded.

The MPO also has the right algebra. Its intermediate channels propagate cumulative charge and charge squared. At site j they implement

\[
Q_j=Q_{j-1}+q_j,\qquad
Q_j^2=Q_{j-1}^2+2Q_{j-1}q_j+q_j^2,
\]

with `q_j=Sz_j+(-1)^j/2` and `q_j^2=I/2+(-1)^j Sz_j`. Exiting into the completed channel sums these link energies together with the mass term and the nearest-neighbor hopping. The MPO includes the last cumulative-charge square; it vanishes in the even-N, total-Sz-zero sector used here, matching the ED link convention. The physical gap is `(E1-E0)/(2 sqrt(x))` in units of g.

The target at m/g=0.25 is correctly taken from Table V of [Byrnes et al.](https://arxiv.org/pdf/hep-lat/0202014): `0.51918(5)` for the vector binding energy in their convention. This is a neutral bound-state energy; it is not the bare electron mass or a free-photon mass. At m=0 the target is the boson mass `g/sqrt(pi)`.

Keep two statistical distinctions explicit:

- The twelve-fit range measures sensitivity to those twelve extrapolation choices. It is not a guaranteed uncertainty interval or a confidence interval covering all discretization, volume and solver errors. The massive target is outside that range, although the registered ±3% test passes comfortably.
- This is a registered rerun of previously examined benchmarks. The massive point was a holdout in the earlier exercise; it is no longer an untouched holdout after its ED result and target have been inspected. Git establishes artifact order, not independently when an operator first saw every result. The original ED registration and its results first enter the available branch history together.

The DMRG agreement is strong support for the numerical implementation. It is not an observational confirmation of GRUT, and the last digits of the primary extrapolation are not the method's accuracy.

## 2. The detector-history result is useful; matching condition (14) is not established

**Required correction:** `RESULTS.md:27` says that matching condition (14) holds. The integrated theory's condition concerns controlled total-variation agreement of microscopic and EFT histories over a declared state and protocol class. This packet checks selected spectra and a few two-time vacuum histories. In the massive case, the comparison is deliberately against an insufficient Gaussian approximation, not against an independently computed interacting-continuum history. There is no interacting-continuum probability error bound here.

Replace that sentence with:

> The adopted Schwinger lattice model reproduces known continuum spectrum benchmarks and gives the expected finite-lattice tests of Gaussian sufficiency and measurement disturbance. These are checks toward the interacting matching obligation; they do not yet establish equation (14).

At zero fermion mass, bosonization makes the continuum electric-field sector a free massive boson. Gaussian sufficiency applies to the specified vacuum and linear-field Weyl instruments; it is not a claim that all states of a free field are determined by two-point functions. The lattice residual decreases on three resolutions, but three points do not certify its zero continuum limit. Also, the listed physical lengths are `5`, `14/3`, and `5` in units of `1/g`, not exactly fixed volume as the original P-C wording requests. The endpoint comparison does keep the same length.

The fluctuation mismatch mentioned in HANDOFF.md can be probed cheaply. I ran an **exploratory, post-result** control at the two smaller lattices. Keep the zero-mass detector strength c=2 and choose the massive detector strength to match Var(cF). This changes the detector strength, not the state or Hamiltonian. It does not establish a new holdout or continuum result.

| Lattice | Zero-mass residual | Massive residual at original c=2 | Massive residual at matched variance |
|---|---:|---:|---:|
| N=10, x=4 | 0.007634 | 0.007875 | 0.027730, c=2.89373 |
| N=14, x=9 | 0.005514 | 0.009180 | 0.029368, c=2.83009 |

The larger massive non-Gaussian residual survives this control. That strengthens the finite-model interpretation without completing the interacting continuum bridge. [Replay and full results](histories_review.py).

## 3. The light-cone arithmetic reproduces; its bound needs tighter scope

For the particular continuous-time cubic-lattice actions tested, the massless expansions are correct:

\[
C_e=\frac38-\frac12 S,\quad C_\gamma=-\frac18S,
\quad S=\sum_i n_i^4,\quad
C_e-C_\gamma=\frac38(1-S)\ge0.
\]

Equality occurs on the axes. The sampled 77.63% of directions with positive electron coefficient also reproduces. But both massless dispersions approach the same low-momentum speed as `ka -> 0`. These calculations show **different finite-spacing dispersion**, not the absence of a common limiting speed, nor a calculation of an exact signal-front causal cone.

The reported `2.95e-32 m` is a body-diagonal, soft-photon Cherenkov estimate with an inferred 2.3 PeV electron energy. It requires a physical preferred lattice, the stated Wilson and plaquette actions, a specified orientation, and control of other velocity corrections. A regulator sent to zero is not automatically a physical lattice constrained this way. Electron speed greater than the photon speed at the same finite momentum is not by itself the soft Cherenkov criterion; that criterion compares to the soft-photon speed. A complete kinematic treatment can also require hard emitted photons.

[Li and Ma](https://arxiv.org/html/2204.02956v1) analyze detected gamma rays, not directly detected PeV electrons. The 2.3 PeV value is inferred from their inverse-Compton emission model, with a more conservative single-scattering lower estimate near 1.2 PeV. Their article supplies the physical argument and a quadratic-dispersion estimate, not this anisotropic Wilson-lattice bound as a published result.

The report mentions an `m*a` correction but treats it too lightly. With positive tree-level Wilson mass and r=1,

\[
E^2=m^2+(1+ma)p^2+
\left(\frac14-\frac S3\right)a^2p^4+\cdots,
\qquad
v_e-1=-\frac{m^2}{2p^2}+\frac{ma}{2}+C_e(ap)^2+\cdots.
\]

At the quoted body-diagonal bound, the omitted `ma/2` is `3.82e-20`, while the retained mass deficit is `2.47e-20`. Its ratio is 1.55. It cannot be neglected unless an explicit matching/tuning condition removes it. Keeping it gives an illustrative leading-order body-diagonal soft threshold near `1.45e-32 m`; that is not a replacement universal bound. The action, counterterms, orientation and radiation kinematics must be specified before a physical exclusion is claimed.

The defensible statement is that a specified physical Wilson/plaquette lattice has experimentally constrained relative dispersion. Neither number is presently a model-independent GRUT spacing limit.

## 4. Gravity and chirality need theorem-level qualifications

The Newton, polarization and exchange calculations are useful standard-physics benchmarks. The claim that a massive graviton is categorically excluded by its 3/4 light-bending factor is too broad. That factor belongs to the linear, unscreened Fierz-Pauli flat-space comparison with Newtonian normalization. Nonlinear screening can invalidate that linear approximation near a source; see [Babichev and Deffayet](https://arxiv.org/abs/1304.7240). Even the report's cited [Kogan et al. paper](https://arxiv.org/abs/hep-th/0011138) concerns circumstances that circumvent the discontinuity. The trace-coupled scalar comparator is also a specific theory, not all scalar-tensor gravity.

[Cassini](https://pubmed.ncbi.nlm.nih.gov/14508481/) measured a radio-link frequency-shift effect associated with solar time delay and constrained PPN gamma. Comparing gamma is appropriate within the stated metric/PPN assumptions. Reporting tens of thousands of Gaussian sigmas is an extrapolation of a local uncertainty, not a calibrated discovery significance. Say those unscreened comparator predictions are incompatible with the measured gamma under those assumptions; do not say Cassini detects a graviton or excludes every massive-gravity completion.

Nielsen–Ninomiya constrains specified local lattice realizations with its symmetry and regularity assumptions. The eight-node example and deformations do not prove that every local Hermitian lattice theory has an unavoidable vector-like physical spectrum. [Lüscher's construction](https://arxiv.org/abs/hep-lat/9802011) changes the realization of chiral symmetry via the Ginsparg–Wilson relation. Adopted operational access alone has not selected among the known alternatives here. A historical review's statement about what had not been achieved in 1999 should not be used as a current universal impossibility theorem.

The Weinberg–Witten stress-tensor theorem concerns massless spin greater than one carrying nonzero four-momentum when a suitable Lorentz-covariant conserved stress tensor exists. Its current theorem is a distinct statement. It is not merely a ban on composite gravitons, and does not independently derive a rule that every microscopic lattice must lack a stress tensor. A lattice already changes the Lorentz hypothesis; the burden is to demonstrate a consistent desired low-energy limit. The original source is [Weinberg and Witten, *Limits on massless particles*](https://doi.org/10.1016/0370-2693(80)90212-9).

## 5. Reproduction and documentation corrections

The actual GitHub validation script still opens `../qed/spectrum_cache*.json`, while the PR directory is `ed/`. Resolve paths from `Path(__file__)`, target the repository's `ed/` directory, and let replay outputs go to a separate destination. The script should also exit nonzero when its `ok` flag is false; currently only an exception yields a failing process status. No sealed script or result was edited during this review.

The reports still contain outdated present-tense statements that no files were committed, that the MPS solver does not exist, and that the current spectrum tally is 0/2. Preserve the historical failure section, but label it as historical and update the current summary. An initial failed run followed by a successful amended run is the accurate record.

For this review's additional reproducible checks, from the repository root with the recorded dependencies installed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_qed_review_20261107/replay.py --output /tmp/grut_pr14_replay.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_qed_review_20261107/histories_review.py --output /tmp/grut_pr14_histories.json
```

The four-case original validation was executed unchanged from the extracted ZIP's `grut_gp/qed_mps/` directory. Its source is identical to the PR's and its raw log is archived here. The frozen analysis and the three algebra/dispersion scripts were also executed in that disposable extraction so their original outputs were not overwritten in Git.

## Replacement assessment for the project record

The contribution earns a stronger interacting benchmark: an adopted 1+1-dimensional gauge theory fits the operational-history formalism, its supplied solver passes bounded reproduction checks, and its registered extrapolation recovers known neutral bound-state energies. Finite history tests and an exploratory matched-variance control support the expected loss of Gaussian sufficiency with nonzero fermion mass. Full interacting-continuum detector matching remains open. The gravity and dispersion calculations constrain declared conventional comparators; they neither derive GRUT's field content nor provide a unique GRUT law. Retain these results with those scopes and keep publication paused.
