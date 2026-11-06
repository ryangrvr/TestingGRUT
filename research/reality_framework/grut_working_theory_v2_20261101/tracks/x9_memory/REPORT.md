# X9 · Which parts of the preparation survive?

The tested early H55 sterile phase has a small effect on bulk
counts at the archived diagnostic time, despite a large
difference in instantaneous work. The continuously coupled
bath gives a separate result. The total daughter count can
lose that dependence while the daughter channels retain it.
With live scalar feedback, matched joint preparations can
also partition energy differently. These statements concern
different observables and different levels of the model.
[Equations and scope](MODEL.md)

## Early H55 production

Both branches evolve from the X5 slow-roll state. At
$t=0.16H_i^{-1}$ the alternate branch changes only the
work-carrying Bloch component in the compact momentum
window. The fermion energy fraction there is
$7.95\times10^{-5}$. Scalar and metric data, energy, pressure
and free occupations agree exactly in the checkpoint records.
The first post-change work has opposite sign. At the fine
radial grid it shifts from $0.03146\,\rho H_i$ to
$-0.02552\,\rho H_i$. The largest early work separation is
$2.29\,\rho H_i$. The largest separation in early fermion
energy fraction is $3.84\times10^{-4}$ of the total density.
Both branches keep the Pauli bounds and the reciprocal
pressure-work budget.

The following contrasts are relative percentages averaged
over the final $0.10H_i^{-1}$. The phase comparison is
between matched branches at the same numerical resolution.

| Resolution | Total comoving count change | Largest individual count change | Fermion energy density change |
|---|---:|---:|---:|
| 2,880 radial bins | 0.00704% | 0.0311% | 0.00381% |
| 5,760 radial bins | 0.0111% | 0.105% | 0.00344% |
| 2,880 bins and half timestep | 0.00744% | 0.0313% | 0.00363% |

The individual-count column includes physical and comoving
counts. The late absolute species densities still move under
radial refinement, by up to 0.88% in the last window for the
unperturbed branch. This supports a subpercent bound for the
tested perturbation at each resolution. It gives less
precision for the sign and magnitude of a small residual
difference. The inherited terminal N4 refinement miss is a
separate diagnostic. [Full production](results/production.json)

An independent calculation fixes the actual scalar midpoint
history and maximizes the count and energy differences over
all allowed signs of the initial Bloch component in this
momentum window. It gives an upper envelope of 0.0381% on
the window-averaged total number and 0.00527% on the averaged
total fermion energy. The largest endpoint number envelope
is 0.0660%. The fixed-history envelope uses positivity and
unitary mode rotations. The nonlinear coevolved pair above
tests the feedback separately.
[Admissible phase envelope](results/envelope.json)

## Continuously coupled bath

X9 represents the two signed vacuum lepton–Higgs spectra
as a Hermitian quadratic bath. The real mass is matched on
the zero-drive shell. The resonance widths agree with the
Born widths within 0.36% in the checked examples. A dense
matrix exponential checks the finite propagator. A direct
weak-coupling integral checks the predicted bath amplitudes.
The interaction stays on during the state comparison.

The two joint preparations have the same sterile covariance,
bath energy-bin occupations, total number and energy. They
differ by the sign of a coherence between a sterile-centered
packet and a bath packet in the positive spectral subspace.
For six H55 modes on a common early background, that
difference generates nonzero daughter energy and count
contrasts. The test uses one positive excitation per mode as
a normalization, with the common interacting vacuum
canceling in the difference. It uses a prescribed scalar
trajectory from $0.16$ to $0.42H_i^{-1}$.

The first grid and sparse measurements have substantial
resolution sensitivity. The dense comparison uses the same
physical averaging window at every resolution. In it, the
frequency and timestep refinements change all registered
window contrasts by less than 0.25% of their comparison
scales. The high-frequency cutoff remains relevant. Its
change moves the bath energy contrast by up to 7.36% of
the registered comparison scale. That scale is the larger
of the magnitude of the final-window difference and 10%
of its peak over the sampled history. The strongest bath
energy difference is 0.0834% of the matched initial
excitation energy, with a stable sign under the tested
cutoff, frequency and time refinements. This is a resolved
contrast in the vacuum model for selected modes. The
integrated H55 daughter yield remains open.

The dense Padé runs conserve norm and close the laboratory
energy plus scalar-work and redshift ledger within
$1.73\times10^{-11}$ relative. The separate continuous
torque integral converges with timestep. The frequency
cutoff and state-preparation assumptions remain explicit.
[Uniform bath](results/bath.json) · [Paired bath refinements](results/bath_refined.json) ·
[Dense bath histories](results/bath_dense.json)

## When memory is erased from an observable

The stationary control follows a bath with recurrence time
beyond the integration interval. The two preparations differ
in sterile count by up to 0.614 during the transient. In the
late window both their sterile counts and their total positive
daughter counts agree within $1.13\times10^{-6}$ per initial
excitation. Their daughter chirality counts still differ by
0.1087. Independent frequency and timestep changes reproduce
this channel contrast. The exact stationary energy weights of
the two conjugate states agree. The late outgoing channel
comparison is taken before the recurrence of the finite bath.
It is a chirality diagnostic of the quadratic model.

In the finite live scalar–sterile–bath control the same
matching procedure gives identical initial energy and
counts. The evolving scalar and bath energies differ at
the endpoint by 0.0965 and 0.4556 in the control units.
The total energy ledger closes to roundoff, including
interaction energy. This demonstrates why stationary
loss of dependence in a count does not settle an evolving
production problem. [Stationary and live controls](results/live.json)

The low-scale working model can now distinguish an early
phase variable that matters little for a bulk count from
a joint phase that is visible in decay channels. The
continuously coupled bath still needs a physical interacting
initial state and feedback on the spatial H55 source. A
cosmological relic prediction requires those conditions
and a controlled particle handoff. The inherited strict
ledger contains zero accepted points out of 135.
[Working theory](../../WORKING_THEORY.md) · [Verification](results/verification.json)
