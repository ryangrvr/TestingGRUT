# Continue GRUT from X6

5 November 2026 · TestingGRUT branch codex/grut-h55-slab-transport-20261105.

Read [REPORT.md](REPORT.md), [MODEL.md](MODEL.md) and [REPLAY.md](REPLAY.md).
The [X5 handoff](../x5_coupled/HANDOFF.md) links the wider theory and
the H55 source state. Scientific acceptance is false and the strict
relic scan is 0/135.

The four-component slab recovers the X5 homogeneous helicity
equations on the same momentum quadrature. It evolves the
scalar gradient and all occupied fermion orbitals with a
reciprocal mass interaction and a common pressure-work budget.
The metric is an isotropic spatial average. This is a finite
intermediate reduction.

**Source gate.** Matched integrated work changes by 0.148% on
the first wider cutoff and 0.021% on the next. Instantaneous
work changes by 6.41% and 3.91%, respectively. Halving the
wide timestep gives only 0.376%, and the fine transverse
comparison gives 0.294%. A second-order state preparation
reduces one cutoff comparison to 4.87%. Keep the 1% phase-sensitive
cutoff criterion open. The exact static edge decomposition
shows that a local gradient subtraction matched at one
harmonic cannot remove the missing-neighbor term across
harmonics. The sharp cutoff static response scales as
$|k|$ at leading order, where the local term scales as $k^2$.

**Transport gate.** In a controlled compact-pulse preparation,
retain the full covariance in a free tail with stated decay and
elastic dephasing. The first allowed volume-averaged boundary
is 18.65 after the pulse with 92.02% of parent energy remaining.
Two later boundaries reproduce the daughter and parent history.
The decay-only pressure envelope fails before half the parent
energy is lost. The explicit control dephaser is not a model
of the H55 collision kernel.

**H55.** The actual archived interface remains active. Its
dressed energy gap is $9.76\times10^{-4}$, and free-axis
dephasing removes pressure coherence and ongoing work. At
the 5,760-bin X5 slice the implied comoving fermion-energy
derivative jump is $-0.10831 H_i$ after normalization by
comoving total energy. The terminal N4 refinement miss also
remains. Keep the original H55 state and its recorded scope.

## Next calculation

Replace or complete the sharp spatial regulator and match the
static response across at least three low harmonics with the
same force and stress prescription. Re-run the helicity,
spatial ODE, momentum, expansion and phase-sensitive cutoff
checks. Only a passing regulator study can justify evolving
a spatial approximation to the actual H55 history. The slab
retains one spatial direction, so test that reduction against
spatial data before interpreting a cosmological output.

Then propagate the energy-basis coherence and the moving-basis
connection through the overlap. Derive or constrain the
collision terms from the H55 interactions. Use the pressure,
momentum, work and energy criteria together, and move the
handoff boundary to test the downstream result before any new
release or relic calculation.
