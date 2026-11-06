# X11 · H55 production, a fixed-scale relic, and axial CP matching

The H55 card now has a live three-axis scalar history with reciprocal
production of all four sterile species in a cell-local approximation.
It reaches the archived interface time. Native occupations carry the
number and energy into a conditional kinetic and relic run. The
pressure and work criteria fail, and the dark abundance falls short
of the benchmark on the tested continuations. A separate axial CP
postulate makes definite decay predictions and fails the tests below.
[Protocol](PROTOCOL.md) · [Equations](MODEL.md)

## Production through the waterfall

The source starts at the archived scale $H_i=5.5\times10^{10}$ GeV.
It uses the H55 Gaussian seed, potential and masses, with the
fermion and scalar forces calculated from one interaction. A common
FRW metric uses the volume-averaged stress. The model is collisionless
in the sterile sector during this stage. It retains three-dimensional
scalar gradients and site-local sterile covariance. Fermion streaming
and the lepton–Higgs bath are separate physical inputs.

| Time $tH_i$ | $3^3$ sterile energy fraction | Calculation |
|---:|---:|---|
| 0.160 | about 0.01% | Before the waterfall |
| 0.200 | 47.34% | Main production |
| 0.220 | 59.66% | Continuing coherent transfer |
| 0.350 | 70.63% | First tested kinetic boundary |
| 0.550 | 69.58% | Reference kinetic boundary |
| 2.9981 | 53.33% | Long three-axis evolution |
| 2.9981 | 51.08% | Refined homogeneous comparison |

The quoted fractions use the free mass energy before reallocating
the derivative subtraction. They describe the coupled stress and
production history in this finite model. The direct X5 homogeneous
replay agrees with the one-cell version in all checked columns.
The reference normalized total and reciprocal quantum ledger errors
are $2.56\times10^{-10}$ and $7.20\times10^{-11}$. The long spatial
total ledger error is $4.07\times10^{-10}$. All terminal free spectra
are Pauli feasible and reproduce the native number and energy.

| Change at $tH_i=0.55$ | Sterile energy fraction | Relative accumulated transfer change | Signed work trace discrepancy |
|---|---:|---:|---:|
| Reference, $3^3$ and 360 radial bins | 69.58% | — | — |
| Half timestep | 69.57% | −0.013% | 6.51% |
| 720 radial bins | 66.89% | −3.18% | 76.40% |
| Wider cutoff at the smaller timestep | 69.63% | +0.077% | 12.84% |
| $5^3$ scalar cube | 69.60% | +0.021% | 5.65% |

The wider-cutoff row compares against the half-timestep baseline at
matched radial spacing. At the original timestep the wider run
fails the midpoint iteration before the waterfall. The signed
trace discrepancy is the maximum difference on common output
times divided by peak baseline work per density. The sampled
instantaneous work contains frequencies above the output cadence.
The discrete impulse and volume-work ledger provide a separate
check on the accumulated transfer.

The two live cubes retain only 0.00870% and 0.0617% of the initial
$\chi$ gradient power. The scalar-only $17^3$ comparison retains
7.76%. In the reference run a null ray could cross about 20 cells
of the original 64³ scalar grid between $tH_i=0.18$ and $0.35$.
The long run differs from the homogeneous result by 4.40% relative
in the late sterile energy fraction. The radial and spatial source
uncertainties prevent a controlled physical H55 production history.
[Production results](results/production.json)

## Conditional abundance at the fixed scale

The continuation starts with the native, volume-averaged occupations.
It keeps their shell weights and adds empty high-momentum bins for the
later kinetic production. The energy match follows from carrying the
full represented energy, including the spectator and dark radiation.
A harmonic partition supplies two pressureless scalar reservoirs.
At $tH_i=0.55$ the fermion-only pressure jump is
$2.36\times10^{-5}$ of total density. The complete pressure jump
is 14.34%, and the residual scalar work remains resolved. At
$tH_i=2.9981$ the fermion-only pressure jump is still
$-2.58\times10^{-5}$. These fail the registered $10^{-6}$ gate.

The first conditional run exposed a 0.022 error in $\ln a$ between
the release clock and its background. The corrected release uses
one interpolant for scale factor and Hubble rate. The kinetic
expansion also uses the pressure in its second-order update.
The corrected reference release reduces the clock error to
$1.03\times10^{-11}$. The native kinetic energy exchange error is
near roundoff, and the refined pressure-work budget error is
$6.71\times10^{-6}$. These are checks on the conditional model.

The hidden release holds the M2-X holder product $c_0=1.55$,
zero holder offset and $e=10^{-6}$ at H55. The holder product is
applied to the coupling at the archived spectator amplitude.
The Q-matter population map receives the calculated visible entropy
and charge yield. The reference result has $S/Q=2.9054$.

| Source and handoff time $tH_i$ | Complete pressure jump divided by density | $K_{\rm vis}$ | $Y_Q$ | Conditional $\Omega_Xh^2$ |
|---|---:|---:|---:|---:|
| $3^3$ at 0.350 | −28.85% | 10.7846 | $3.1241\times10^{-9}$ | 0.06419 |
| $3^3$ at 0.550 | +14.34% | 10.5336 | $3.1523\times10^{-9}$ | 0.06468 |
| Same source at 0.550, half kinetic step | +14.34% | 10.5332 | $3.1524\times10^{-9}$ | 0.06469 |
| Same source at 0.550, 720 production radial bins | +2.84% | 11.3431 | $2.8024\times10^{-9}$ | 0.05766 |
| $5^3$ at 0.550 | +14.32% | 10.5293 | $3.1542\times10^{-9}$ | 0.06472 |
| $3^3$ at 2.9981 | −38.84% | 10.9668 | $2.9471\times10^{-9}$ | 0.06044 |
| Homogeneous at 2.9981 | −28.00% | 11.6155 | $2.6926\times10^{-9}$ | 0.05531 |

Halving the kinetic step changes the reference abundance by 0.0053%.
Moving the boundary to the end of the long spatial run changes it by
−6.57%. Doubling the production radial bins changes it by −10.86%.
The reference integrated vector count changes by 0.31% between two
readouts. Its individual mode readout changes by up to 21%.
The Planck comparison for all cold dark matter is
$0.120\pm0.001$ under base $\Lambda$CDM. The production,
handoff and population assumptions remain part of this inference.
[Continuation results](results/continuation.json) · [Planck](https://arxiv.org/abs/1807.06209)

## Tested GRUT candidate law

The proposed matching law makes the sterile Yukawas real in the
specified inflaton interaction basis. The vacuum scalar-dependent
mass terms then supply the phases in the positive-mass basis.
The physical vacuum phases on H55 are about $-6.52\times10^{-5}$
and 0.24447 rad modulo $2\pi$ for the inflaton and waterfall
channels respectively. The normalized CP decay invariants between
N2, N3 and N4 vanish. The three connecting N1 to that triplet
have magnitude 0.242105 with signs set by the discrete parities.
These relations remove the independent fitted leptogenesis angle
within the postulate.

With all parities positive, the tree-level light-neutrino mass
sum is at least **0.19888 eV** for any real texture with the four
inherited widths. The fixed texture gives 0.19891 eV and the
repository's one-loop estimate gives 0.15625 eV. Both exceed the
Planck-plus-BAO comparison limit of 0.12 eV under its cosmological
assumptions. The bound is obtained from the trace of the seesaw
matrix, independent of the H55 production calculation.
[Derivation](MODEL.md#axial-cp-matching-postulate) · [Planck](https://arxiv.org/abs/1807.06209)

Keeping the existing N2 parity and using its later seesaw
normalization gives a tree mass sum of **0.05881 eV** on the
fixed texture. Carrying those widths and the axial CP invariants
into the corrected conditional reheating history gives
$\eta_B=1.87\times10^{-8}$, about **30.6 times** the measured
benchmark. At twice the inherited washout it is still 13.0
times the benchmark. The numerical baryon test uses the
unflavoured kinetic approximation. The mass sum and exact
invariant relation have separate checks.
[Axial tests](results/axial_law.json) · [Baryon runs](results/continuation.json)

The finite production history, the conditional relic numbers
and the tested matching postulate have distinct statuses.
The physical H55 handoff and dark abundance remain unresolved.
Neither axial matching variant supplies a viable law on the
tested input and kinetic model. The inherited strict candidate
ledger is **0/135**. The evidence verifier passes 108 checks,
including a fresh axial replay and the direct X5 homogeneous run.
[Replay](REPLAY.md) · [Evidence](results/verification.json)
