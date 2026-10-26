# Two low-scale releases: full debris versus the stated lump benchmark

Codex / AInstein, 25 October 2026. Parent 373d37e. Review
f1aa7df968cb3dc90ca32f4f32aeb0a819ddfce4, issue #13 comment 6404819222.
Ownership: Codex computes; Claude reviews (comment 6404864387).

## Fixed cases and action

No coupling, initial yield, phase operator or reheating history is changed.
Use c_H=10, H_i=434867798206.6179 GeV and the same aligned radial minimum,
zero physical velocities and instantaneous radiation release as 373d37e.
Candidate 10 has delta=pi/6; candidate 17 has delta=5pi/24. Their published
48-member source-phase ensembles have Gaussian sigma=0.001 radians and
seeds 2026102510 and 2026102517. The common phase remains supplied.

Use the parent full X/Y/Sigma action and conformal canonical equations.
Freeze phi at zero only in this early high-field diagnostic. The later
light-scale benchmark is an explicitly separate supplied approximation.

## Required calculation

1. Recompute the two nominal releases and late-window controls. Project Y
   about its forced response onto the two real mass eigenmodes, using
   the co-rotating Sigma-phase frame. Measure Sigma radial action separately
   from its charge-carrying angular motion. Apply the same projections to
   every published ensemble endpoint; assess phase and late-window error.
2. Count Y->3X, Sigma_radial->AA and Y_heavy->Y_light A. Each vector
   ultimately gives one X/anti-X pair. Heavy-to-light cascades conserve
   Y particle number, so subtract their energy from the parent budget.
   Record a conservative cascade bracket as well as a central ledger.
3. Average decay injection over exponential decay in radiation. For late
   nonrelativistic parents the mean scale factor is
   sqrt(pi H_i/(2 Gamma)), not the scale factor at H=Gamma. Check delayed
   vector decay with redshifting momentum/time dilation and energy balance.
4. Recompute the same case-L benchmark used by the review, with its quoted
   masses, conversion and 4->2 coefficients. Solve its threshold at actual
   ensemble epsilon and both ymult=1 and the raw 5.3/5.7 yield factors.
   No entropy source is thereby implemented: ymult=1 models normalization
   before hidden freeze-out; raw ymult models a hypothetical normalization
   afterward. Without a real dilution history the abundance remains wrong.
5. Compare the full ledger, its conservative lower bracket and Y-only
   control with that threshold. If a case passes, do not presume lumps
   actually form or erase anything. If both fail, close these published
   candidate releases under the stated benchmarks. A finite set of
   releases does not prove failure at every real phase or changed history.

## Checks and limits

Use an independently derived real-matrix mode projection and derivative
control; compare number conservation across late windows, integration
tolerance and ensemble endpoints. Check thermal integrals with at least
two fixed quadratures and selected adaptive reference integrals; check
energy, entropy and charge constraints and the reaction-step resolution.
Keep numerical control passes separate from physical gates.

The thermal benchmark presumes complete equilibrium, a lump proxy fixed
at Q=200000, and quoted cold conversion rates. It is not a fragmentation,
solitosynthesis or full Boltzmann simulation. Vacuum perturbative decays
are postprocessed on the prepared homogeneous releases; backreaction,
coherent production and thermalization are not silently certified.

The review's earlier 1060-vector-pair estimate is withdrawn. The evolving
mass and initial-state convention must remain explicit; a constant late
mass is not substituted through inflation. An omitted nonnegative source
cannot rescue a failure of a conservative injection bound.
