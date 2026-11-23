# Separate critical review of GRUT-D2

23 November 2026. This is review by a separately tasked AI within this Work
run. It is not an external RAI session, external replication, human expert
review, or peer review. The reviewer constructed full physical-spin matrices
independently of the builder's reduced wall-path implementation.

## Verdict

**The scoped finite-window mathematical construction passes this review.**
The physical Pauli Hamiltonian is local, time independent and exactly
equivalent to an orientation bit and a finite XX path of domain walls and
reservoir excitations. It orders a symmetry-neutral family, transfers initial
defect-count information to a reservoir, and permits bounded recovery after
one specified endpoint flip. Both orientations are retained.

**The strict target-free development gate is not met.** Ferromagnetic
alignment, a resonant interface, an initially pure reservoir, the size/time
scaling and the collective count readout are supplied. A passive reservoir
record is established; record-dependent corrective control, an emergent
observer, consciousness, a new GRUT physical restriction and cosmic planning
are not established. The large pure reservoir buys the finite lifetime.

One frozen control forecast was **falsified** during review. A maximally mixed
active system can transfer initially supplied reference correlations into its
bath, even though its active marginal remains stationary. The original
protocol is preserved and the correction is in
[CONTROL_CORRECTION_01.md](CONTROL_CORRECTION_01.md).

## Independent physical-matrix checks

Run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 independent_check.py
```

The script starts from the physical Ising, three-spin wall-hopping,
three-spin interface and XX-reservoir terms in MODEL.md. It builds full
Hilbert-space matrices for all nine combinations `n=2,3,4`, `L=1,2,3`
(dimensions 8 through 128). It explicitly enumerates initially populated
physical spin strings and propagates their mixture using a dense
eigendecomposition. Separate small one-particle exponentials also test the
spectral formula. These finite matrices do not use the free-fermion reduction
to generate their reference probabilities.

Checks include the exact physical-to-wall basis permutation; the factors and
signs in the interface; conserved anchor and global-flip symmetry; complete
initial-count/bath-count joint distributions; unconditional organization;
post-flip expected wall number; phase-echo organization; present-system and
count-record invariance under the bath kick; subsystem entropies; initial
connected/disconnected energy moments; and the echo energy variance.

The physical matrix and encoded XX matrix agree exactly in these cases.
Four readout times per case were checked, plus separate kick and echo
histories. The largest observed discrepancy was `4.09e-14` for echo energy
variance; the largest expected-wall discrepancy was `3.34e-15`. All are below
the frozen `1e-11` small-matrix tolerance. Full printed results can be retained
as INDEPENDENT_CHECK.json by the packet runner.

These checks independently validate the reduction on the listed small
systems. They are not a separate full-matrix simulation of the much larger
production reservoirs, an interval-arithmetic evaluation of floating-point
results, or external replication of the whole packet.

## Review of the uniform bound

For the half-line hopping propagator, the image formula and Bessel recurrence
give, for integer `i>=j>=1`, a phase times

\[
 \sum_{r=0}^{j-1}(-1)^r
 \frac{i-j+2r+1}{t}J_{i-j+2r+1}(2t).
\]

The Fourier integral for an integer Bessel function gives
`|J_k(x)| <= 2/pi` for `k>=1`: the absolute integral of `cos(k theta)` on
`[0,pi]` is 2. The positive indices in the displayed sum add to `ij`.
Consequently `|U_ij(t)| <= 2ij/(pi t)` and the source-block Frobenius norm is
at most `2*S2/(pi*t)`.

A path affected by the reservoir's far boundary has at least `K=2(L+1)`
hops between the scored source sites. The exponential-series tail is bounded
by its first term times a geometric tail. Throughout the scored window,
`x=2*t`, `x/K<1/4`, and `e<3` give the conservative entrywise bound
`epsilon <= 2*(3/4)^K`. Its Frobenius contribution is at most `m*epsilon`.
This is a finite-reservoir bound, not a thermodynamic-limit substitution.

For initial occupation probability one half, the exact Gaussian baseline
formula and an elementary union bound give

\[
 P(N_{\rm wall}=0)=\det(I-\tfrac12 AA^\dagger)
 \ge 1-\tfrac12\|A\|_F^2.
\]

Using the stated family, the rational comparisons with `pi>3.14` certify
the entire declared window above `0.9968`, not merely the 65 sampled times.
This means a single readout at any chosen time in that window has the stated
marginal probability. It does not bound a joint continuously monitored
history, a no-jump event, or success at repeated readouts; measurements can
disturb the subsequent evolution.
Because the bath initially has no excitations, `K=N_wall+N_B` is conserved;
the declared decoder `K_guess=N_B` succeeds exactly when `N_wall=0`.
It records the initial **count**, whose entropy is binomial-count entropy,
not the entire original wall pattern. It is a specified decoder, not an
assertion of optimal guessing at every time.

For the flip, the normal covariance update is valid even though the state is
generally no longer number-conserving Gaussian. With `k=m-1`, its diagonal
sign matrix has entries `+1` before k, `0` at k and `-1` afterwards. At equal
write/recovery intervals, the transferred source block is

\[
 V=-A(2\tau)+A(\tau) B A(\tau),
 \qquad B=\operatorname{diag}(2,\ldots,2,1).
\]

Since `||B||op<=2`, the direct norm bound
`||V||F <= ||A(2tau)||F + 2||A(tau)||F^2` is valid. Together with the endpoint
column bound it yields the stated recovery certificate above `0.9922`.
This argument does not assume that a cross term has a favorable sign; such a
sign assumption would be false in some cases.

The reviewer deliberately tested an invalid alternative: applying the
number-conserving vacuum determinant to the normal covariance after the
physical X flip. For `n=4,L=3,tau=0.731,s=0.413`, the full matrix gives
organization `0.1104317370`, whereas that determinant gives `0.1073053558`.
The difference `0.00312638` demonstrates why the packet correctly uses
expected-wall bounds after the kick.

## Falsified hot-reservoir forecast

Take one wall and one reservoir spin, each initially a fair bit, with an
inert reference labelling the initial wall bit `K`. The active wall/bath state
is `I/4`, but the reference and active system are correlated. At `t=pi/2`
the exchange swaps the wall and bath occupations. The following are all true:

| Quantity | Initially | After exchange |
|---|---:|---:|
| Active state | I/4 | I/4 |
| Aligned probability | 1/2 | 1/2 |
| I(initial K : bath count) | 0 bits | 1 bit |
| Fixed decoder success | 1/2 | 1 |

The full physical matrix reproduces this result. It invalidates the frozen
"no new count record" forecast with the declared past-state tag. It does not
invalidate stationarity of the active marginal or the absence of additional
ordering. A fully maximally mixed state including the reference would have no
correlations to transfer, but would also lack this valid initial tag. These
are different preparations and must not be conflated.

## Interpretation and control limits

1. **An inserted relation remains.** The two orientations are symmetry
   equivalent, but alignment is explicitly favored. The conserved microscopic
   anchor supplies the final orientation. This removes an inserted absolute
   direction, not the preferred organization class.
2. **A sink is not an actuator reading a memory.** Removing the interface
   disables wall transport as well as recording. This cannot isolate a causal
   benefit of stored count information. Its initial mean energy is matched,
   while its variance differs by one half in the declared units. The packet
   reports this difference.
3. **The echo concerns quantum environmental state.** A bath-only diagonal
   phase changes future organization while preserving the current physical
   system and count distribution. This shows that the present system state
   and classical count do not specify the future. It does not show that the
   actuator consults that count. The phase controller and timing are supplied.
   Zero mean work does not remove their preparation or control costs, and
   energy variance can change; both energy diagnostics were checked.
4. **Recovery has a restricted disturbance class.** The endpoint X creates or
   removes a wall and supplies work. One subsequent wall-removal event without
   resetting the bath is demonstrated. No general repeated-repair lifetime,
   distributed damage repair or evolution of a developmental program follows.
5. **Memory readout is supplied.** Measuring reservoir total excitation number
   is collective access to a designated observable. No localized reader,
   self-model or autonomous observer is constructed.
6. **Time and resources remain supplied.** The baseline Hamiltonian is fixed,
   but its Schrodinger time parameter is external. Reservoir size and scored
   lifetime grow as `O(n^3)`; the initial bath has L pure spins. Global entropy
   is conserved and subsystem entropy is redistributed. A fixed finite
   reservoir eventually recurs.

## Why the result is not a cosmic-selection theorem

Finite fixed-H observable probabilities are analytic finite Fourier sums.
A newly attained probability cannot remain exactly one on any open interval
unless it was always one; a convergent full density matrix must have been
stationary. The finite-window approximate claim avoids these obstructions.

Closure alone also permits interacting Hamiltonians with maximally mixed
initial states and no newly generated correlations, provided every reference
is included in that fully mixed state. Additional preparation or dynamical
restrictions are necessary to force developmental histories. This packet
does not provide a new restriction selecting observer-producing universes.

Finally, in the preceding D1 circuit the coherent input subsystem was pure
before and after correction. Its computational-basis outcome entropy fell;
its von Neumann entropy did not. A diagonal mixed input instead produces
classical memory/actuator correlation, rather than their entanglement. Three
scheduled CNOTs alone were not an autonomous clock-and-energy construction.

The defensible result is therefore an exact, reproducible standard-physics
example of autonomous finite-window ordering, a passive internal count record,
and specified disturbance recovery. The stronger GRUT claims remain open.
