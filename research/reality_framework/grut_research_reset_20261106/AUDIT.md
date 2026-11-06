# GRUT research reset

6 November 2026. GRUT names a research program. A conventional
mechanism or a scoped negative result is a valid outcome.
The audit uses [X11 at 5d40afd](https://github.com/ryangrvr/TestingGRUT/commit/5d40afdd87f124953b8a83348e17460f49350604),
its raw production, spectra, kinetic and release files, the X2
acceptance rules, and the current working theory. The small
evaluator cases were [registered](PROTOCOL.md) before calculation.
No new cosmological production history was evolved.
[Working theory](../grut_working_theory_v2_20261101/WORKING_THEORY.md) ·
[X11 report](../grut_working_theory_v2_20261101/tracks/x11_h55/REPORT.md)

## Three assumptions that determine the verdict

1. **The source represents spatial fermion production.** X11
   evolves scalar gradients but confines each fermion continuum
   to one site. Its live grids retain 0.0087% and 0.0617% of the
   archived initial waterfall field's gradient power. The
   homogeneous subtraction and specified vacuum preparation
   do not control a nonlocal stress in that geometry. The first
   bulk radial refinement changes the energy fraction by 3.86%.
   These are separate uncertainties.
2. **The handoff preserves the stress and remaining energy exchange.**
   Native occupations preserve number and energy at the splice.
   The conditional kinetic engine discards fermion coherence and
   replaces the evolving scalar stress with pressureless parents.
   It supplies collision, thermalization and bath evolution after
   that boundary. Both the pressure and work criteria fail on
   the raw records.
3. **The hidden-sector state and abundance are independently
   identifiable.** The fields, masses, couplings, holder state,
   CP texture, release and population map are choices. Earlier
   histories fit abundance through the inflation scale, and
   the spectator amplitude was screened using isocurvature.
   X11 fixes the scale for its comparison but inherits the chosen
   holder and population model. A matching density alone cannot
   identify the mechanism.

## Evidence ledger

| Category | What the record establishes | Earliest dependency |
|---|---|---|
| Adopted physics | Quantum dynamics, the Majorana particle dictionary, Pauli bounds, the FRW metric and standard visible-sector reaction laws enter as inputs. | The action, quantum statistics and choice of spacetime description. |
| Supplied model | Hybrid scalar background, four axial sterile masses, seed, cutoff and subtraction, dark fields and holder, rate tables and population map. | A specified EFT and preparation. No tested closure argument selects them. |
| Mathematical effects | Scalar work and reciprocal force follow from the same vertex. At matched number, energy and pressure, transverse coherence can reverse work. The positive-parity real-Yukawa mass bound is exact at tree matching. | A quadratic driven Hamiltonian, or the specified seesaw texture class and widths. |
| Implementation | X11's one-cell reproduction of X5 exercises the same microscopic kernel. X10's reduced bath passes its finite controls. The new ladder uses a direct complex ODE, an invariant integral and a separate outgoing-momentum integral. | The declared finite models and reference equations. |
| Numerical convergence | The X11 timestep effect on the conditional relic is 0.0053%. Radial production changes it by 10.86% and moving the handoff changes it by 6.57%. X8's converged synthetic spatial result is conditional on a fixed ultraviolet preparation. | The same source, state, regulator and interface for each individual refinement. |
| Uncontrolled physics | Spatial fermion propagation, physical UV matching, sterile–bath dynamics through production, coherent-to-kinetic transfer and the hidden population remain open. | The approximations above, before any abundance comparison. |
| Observational comparison | The fixed X11 reference gives a conditional $\Omega_Xh^2=0.06468$ against 0.120. Its handoff fails. X11's axial CP versions miss the tested mass or baryon benchmark under their respective assumptions. | The full conditional chain and the data model. No fresh observational holdout. |
| Novelty | The contribution is a quantitative diagnosis of when the tested reductions lose information. Axial production, Pauli blocking, spatial phase space and memory from initial correlations have conventional precedents. | The finite-model checks and comparison to prior physics. |

The **0/135** result is the older M2-X free-X acceptance rule,
including the $e_{\rm crit}$ window, on the bracket formed by the
computed nodes at $8.4\times10^{10}$ and $1.3\times10^{11}$ GeV.
H55 is below that bracket. The older audit also records two
M2-Q points under its different rule. The ledger does not
classify H55 or the wider candidate family.
[X2 preregistration](../grut_working_theory_v2_20261101/tracks/x2_strip/PREREGISTRATION.md) ·
[X2 results](../grut_working_theory_v2_20261101/tracks/x2_strip/RESULTS.md)

The wider record retains other mathematical negatives. In
frozen CP-12, a connected $Q=N$ graph has one independent
cycle. At most $R/N$ of its vertices can lie on a cycle of
length at most $R$, which rules out a homogeneous abundance
of microscopic plaquettes under that graph interpretation.
The closure synthesis finds no new selector in the tested
families. The SCOUT handoff places K1 inside a luminal Horndeski
degeneracy, so the relation has an existing comparator. These
were reviewed for scope without fresh replication here.
[CP-12 audit](../cp12_holonomy_locality_gate/REPORT.md) ·
[Closure synthesis](../grut_closure_forcing_20261027/SYNTHESIS.md) ·
[SCOUT pin](https://github.com/ryangrvr/TestingGRUT/blob/ab2da47407bb670d94e2a52c87599fa13fd8ab99/playground/SCOUT_0/SCOUT_0_HANDOFF.md)

## What the checks eliminated

The raw-file replay uses $q$ as comoving momentum and
$(q_+^3-q_-^3)/(6\pi^2)$ for each helicity. Physical momentum
is $q/a$, and densities carry $a^{-3}$. The source and spectra
digests match the pinned kinetic inputs. The free-particle
moments and conditional abundance arithmetic reproduce from
the raw excerpts. The population factor remains an inherited
input to this replay.

| $tH_i$ | $\lvert\Delta\rho\rvert/\rho$ | Fermion-only $\Delta P/\rho$ | Full $\Delta P/\rho$ | Conditional $\Omega_Xh^2$ |
|---:|---:|---:|---:|---:|
| 0.350 | $4.35\times10^{-14}$ | $-8.23\times10^{-5}$ | −28.85% | 0.064189 |
| 0.550 | $2.44\times10^{-15}$ | $+2.36\times10^{-5}$ | +14.34% | 0.064685 |
| 2.9981 | $1.38\times10^{-14}$ | $-2.58\times10^{-5}$ | −38.84% | 0.060435 |

The $10^{-6}$ pressure gate fails for the fermions alone and
more strongly for the complete splice. The raw work is still
active in all three windows. These comparisons eliminate an
energy-only acceptance rule and this implementation of the
handoff at these boundaries.
[Raw excerpts](results/raw_slices.npz) · [Recomputed values](results/checks.json)

For the weak spatial drive, linearizing the phase gives the
pseudoscalar vertex $-(M\theta/2)\bar N i\gamma_5 N$.
For frequency $\nu$ and spatial momentum $k$, the two-body
invariant is $s=\nu^2-k^2$. The Born energy response is
proportional to
$$
\int_{\sqrt{4M^2+k^2}}^\infty d\nu\,\nu |\widetilde f(\nu)|^2
\,s\sqrt{1-4M^2/s}.
$$
The separate lab-frame integral fixes the second momentum
by $\mathbf p_1+\mathbf p_2=\mathbf k$ and evaluates
the pair energy $E_1+E_2$ for every quadrature point. The site-local rule
uses $s=\nu^2$ at the same RMS drive. At $\omega=2.4M$
the central line is allowed for $k=M$ and below threshold
for $k=1.6M$.

| Wave at $\tau M=24$ | Propagating energy / site-local energy | Interpretation |
|---|---:|---|
| $k=M$ | 0.59633 | Ordinary phase-space suppression |
| $k=1.6M$ | $2.14\times10^{-9}$ | The finite pulse's high-frequency tail contributes |

The independent quadratures agree within $4\times10^{-12}$
relative. The long-duration homogeneous response reproduces
$\Gamma_{\phi\to NN}=g^2m_\phi\sqrt{1-4M^2/m_\phi^2}/(16\pi)$
within 0.0145% at $\tau M=48$. Omitting the identical-particle
factor doubles the width. The volume derivative checks
recover the massless radiation limit and finite vacuum and
interaction pressure. The broken pressure check keeps the
vacuum term but loses 75.8% of the pressure. These are evaluator
controls in their specified regimes.
[PDG phase space](https://pdg.lbl.gov/2024/download/db2024.pdf) ·
[Results](results/checks.json)

The exact rotating-mode example also gives two admissible
states with the same number, energy and pressure and opposite
work. With mass axis $\mathbf b=(M,0,p)$, the three matched
observables fix combinations of $r_x,r_z$. The axial torque
depends on $r_y$. Reversing that component leaves the matched
quantities invariant. This excludes a universal instantaneous
closure on those moments for the allowed state class.
A restricted vacuum preparation could still support a useful
late-time closure. X9's specific early reflection already
suggests weak sensitivity of the integrated bulk count, while
its channel diagnostics retain memory.
[X9 report](../grut_working_theory_v2_20261101/tracks/x9_memory/REPORT.md) ·
[Initial correlations in reduced dynamics](https://doi.org/10.1103/PhysRevLett.73.1060)

**With the GRUT label removed**, the portable results are a
checked counterexample to an occupation-only work rule, a
quantified failure of a site-local production approximation,
and the failed physical splice on these particular histories.
The spatial threshold is established relativistic kinematics.
What the audit adds is its use as a matched control for this
pipeline. The finite-state counterexample is a mathematical
effect. The weak vacuum response is a controlled physical
benchmark in its stated regime. The audited results supply
no new evidence about nature.

## Competing mechanisms

Use the same production target first. Compare net energy
delivered to fermions and their outgoing momentum and
helicity spectra when a prescribed scalar drive has ended.
Match the vacuum, background, source duration, input mass
and physically appropriate drive normalization. The axial
and real-Yukawa couplings have different vertices, so a
difference between them cannot be assigned to spatial
propagation alone.

| Mechanism | Evolution, preparation and free inputs | Prediction and cheapest falsifier |
|---|---|---|
| Perturbative pseudoscalar decay | $-(g/2)\phi\bar N i\gamma_5N$. A weak classical drive or dilute parent population, with vacuum fermions. Supply $M$, $g$ and the parent spectrum. | Recover the vacuum width and the invariant threshold at small amplitude. Reject rates that violate the threshold or the identical-particle normalization. |
| Coherent axial production used by H55 | $-(M/2)\bar N e^{i\gamma_5\theta(x)}N$ with a chosen scalar action and state. The cosmology also supplies the holder and the four masses and phase scales. | Nonperturbative energy transfer and possible helicity contrast. First keep the waveform and initial covariance fixed and compare propagation against the local reduction. Then test reciprocal feedback and the physical bath separately. |
| Real Yukawa preheating without a holder | $-(m_0+y\phi)\bar N N/2$. Use an oscillating real scalar, with vacuum fermions before the crossing. Supply $m_0$, $y$ and the trajectory. | A resolved linear mass crossing gives $f_p\simeq e^{-\pi p^2/|\dot m_*|}$ in its single-crossing regime. Recover that spectrum, Pauli bounds and helicity symmetry for a parity-symmetric preparation before testing repeated crossings. |

The weak test already shows that bath coupling, scalar
feedback and expansion are unnecessary for the **existence**
of pair production. The real-Yukawa comparator needs no
CP-odd input for a bulk count. Those simplifications say
less about a net asymmetry, cosmological redshifting or a
self-consistent relic. No necessity test for all H55 inputs
has passed.
[Axial comparator](https://arxiv.org/abs/1508.00891) ·
[Real Yukawa comparator](https://doi.org/10.1103/PhysRevD.62.123516)

For the separate **relic** target compare the final energy
per visible entropy and the perturbation spectrum on the
same visible cosmology. Use adiabatic cold dark matter as
the observational null. Thermal freeze-out starts in
equilibrium with
$\dot n+3Hn=-\langle\sigma v\rangle(n^2-n_{\rm eq}^2)$.
Its mass and annihilation rate must give the abundance
without conflicting with the same rate's other constraints.
Equilibrium-parent freeze-in starts with negligible X and
uses $\dot n_X+3Hn_X=N_X\Gamma_B n_B^{\rm eq}\langle m_B/E_B\rangle$
in the dilute decay limit. With parent and portal fixed,
check that X never equilibrates and compute whether the
resulting spectrum is too warm. Both can be adiabatic
under single-clock conditions. H55's separate release and
population calculation would have to earn its comparison
on that same target. Equal abundance alone cannot
distinguish them.
[PDG freeze-out](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-dark-matter.pdf) ·
[Freeze-in](https://arxiv.org/abs/0911.1120)

## A constraining relation and its comparator

The scalar force and fermion work are linked by
$\partial_\phi H_N$ once the Hamiltonian and state are fixed.
That removes an independent source coefficient in a conserving
implementation. It adds no coupling and is a standard
Hamiltonian identity. A fitted local drain with the same
initial energy is a comparator that can violate the relation
to the computed state. Finite matching and counterterms
still enter both sides.

X11's real-Yukawa postulate removes a continuous phase from
the sterile matching ansatz. It leaves the real texture,
parities, widths, scalar vacuum phases and matching scale
as inputs. Its invariant pair ratios are
$\eta_i\eta_j\sin(\theta_i-\theta_j)$. At fixed positive
masses, column rephasings of the Yukawa matrix can change
those invariants while keeping all four diagonal decay
widths. This comparator shows that the widths and
collisionless source do not select the postulate.
The positive-parity tree mass bound of 0.19888 eV and the
parity-preserving fixed-texture baryon miss are retained
at their stated scope. Running, thresholds and a flavoured
thermal treatment are distinct from the tree relation.
[X11 axial derivation](../grut_working_theory_v2_20261101/tracks/x11_h55/MODEL.md#axial-cp-matching-postulate)

## Observational access

The clearest potential test of the spectator preparation is
CDM–photon isocurvature, $S_{X\gamma}=\delta\ln(\rho_X/s_\gamma)$
on uniform photon density after X is cold. A fluctuation of the
holder can change the mass yield on that slice. A fully
adiabatic one-clock relic gives zero intrinsic isocurvature
in this limit. A detection would distinguish a nonadiabatic
component from that baseline. Other spectator mechanisms
can generate the same signal.

Planck's 95% bound $\beta_{\rm iso}<0.038$ at
$k=0.05\,\mathrm{Mpc}^{-1}$ is for uncorrelated,
scale-invariant CDM isocurvature with the paper's
assumptions. H55 would need the response including entropy,
the dark-matter fraction and the curvature correlation to
predict the spectrum. The earlier work used the Planck
bound to screen the spectator. The mean abundance and
scalar spectrum have also entered earlier fits. Use those
as calibration or consistency checks. Reserve a genuinely
unused polarization likelihood or blinded modes with their
covariance before choosing any new candidate or response.
No untouched observational holdout has been specified.
[Planck inflation analysis](https://arxiv.org/abs/1807.06211) ·
[Adiabatic modes](https://doi.org/10.1103/PhysRevD.67.123504)

## Decisions and the next calculation

- **Retain** the benchmark ladder, the conserving source
  construction, the observable-specific memory tests and
  the exact negative results. Use conventional baselines
  as live competitors.
- **Suspend** H55 and M2 abundance claims and large parameter
  scans. Reopen physical production only with a matched
  propagating source and controlled preparation and cutoff.
  A relic calculation also needs a physical handoff and
  independent observational access.
- **Abandon** this occupation-only handoff at the tested
  boundaries and the claimed uniqueness of a closure or K1
  signature. Close the two axial CP candidates on their
  tested assumptions. A materially different physical input
  would be a new hypothesis.

**Prioritize one small matched propagation experiment.**
First require the propagating quadratic solver to recover
the weak Gaussian benchmark. Reserve a different, compact
pulse for the frozen strong-drive test
$$
\theta(t,z)=1.2\cos^4\!\left(\frac{\pi t}{2T}\right)
\cos(2.4Mt)\cos(1.6Mz),\qquad |t|<T,\quad TM=12 ,
$$
and zero elsewhere. Use flat space, a common free vacuum
before the pulse and a common free particle dictionary
after it. Compare a propagating fermion with the same
mass-phase drive and weights in the cell-local model.
Keep the transverse phase space and match the homogeneous
limit. Hold feedback, bath and expansion fixed outside
both arms. Compare accumulated work and final number,
with the spectrum as a diagnostic. Independently vary
timestep, momentum spacing and cutoffs. Require a
difference to exceed 1% of the larger magnitude and three
times the numerical uncertainty. If matching or cutoff
ambiguity is that large, pause the comparison.

This tests whether the local approximation's weak-limit
failure persists when the phase drive is order one.
Nonlinear harmonics can open additional channels, so the
weak threshold supplies a calibration and does not fix the
strong outcome. A stable large difference closes the local
source route for the tested regime. A small difference
establishes a restricted numerical regime for later tests.
The proposed pulse has not been evaluated in this pass.

## Replay and scope

Run from the TestingGRUT root. The portable run reads the
saved raw slices. The optional full-source path rechecks
the pinned history and spectrum hashes and extracts those
slices again. They contain native spectra, the history tails,
rows around $t=0$, and the kinetic boundary and final records.

~~~bash
RESET=research/reality_framework/grut_research_reset_20261106
OPENBLAS_NUM_THREADS=1 python -B "$RESET/code/evaluate.py" --output /tmp/grut-reset-replay
# With the original X11 work directory present:
OPENBLAS_NUM_THREADS=1 python -B "$RESET/code/evaluate.py" \
    --raw-root /tmp/grut-x11-package --output /tmp/grut-reset-from-raw
~~~

The two phase-space integrations use the same weak vacuum
action. The raw replay tests reconstruction and arithmetic.
Neither constitutes an outside replication or tests nature.
No physical matching or observational holdout was changed
after seeing the results.
[Raw file hashes](results/raw_provenance.json) · [Numerical evidence](results/checks.json)
