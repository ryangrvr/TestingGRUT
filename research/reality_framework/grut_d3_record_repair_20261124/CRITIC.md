# Separate critical review of GRUT-D3

24 November 2026. This review was carried out by a separately tasked AI in
the same Work run, using an independent physical tensor-Hamiltonian
implementation. It is not external RAI execution, external replication,
human expert review, or peer review.

## Verdict

**The scoped mathematical claim passes this review:** a classical pointer
acquired during a writing phase controls which historical label is restored
after one specified erasure, with both restoration and pointer retention
scored at the terminal readout. Both autonomous phases use the same fixed
Hamiltonian. The matched interventions establish a causal dependence on the
record's association with history, rather than merely an improvement from
changing the reservoir or the Hamiltonian.

This advances beyond D2's passive reservoir-count record. The intact and
coin-scrambled cases have identical unconditional active-plus-coin density
matrices and identical complete energy distributions, yet different final
correlations with the inert historical reference. The difference is meaningful
because the historical success observable is fixed before the intervention.

**The result remains an engineered standard-physics mechanism.** The
Hamiltonian supplies the pointer alphabet and conditional writing/repair
relation. The ready memory, chains, fault register, coin, timing and readout
are supplied. The display does not include an autonomous clock, work store,
or physical implementation of the timed pulses. It therefore establishes a
controlled write/fault/repair protocol, not an autonomous closed-universe
history or emergence of an unsupplied developmental program.

## Independent physical reference

Run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 independent_check.py
```

The review script constructs the Hamiltonian directly on the two physical
qutrits and their spin chains, in the physical order `S,M,W...,B...`. It
includes all three levels of each qutrit and all spin configurations. It
uses dense matrix eigendecomposition/exponentials and explicit initial-label
and coin histories, rather than using the builder's reduced grid to generate
the reference probabilities. Inert D and C are represented through those
conditional histories; they do not enter the physical evolution matrix.

The six reservoir-size pairs are `(1,1),(1,2),(2,1),(2,2),(3,2),(2,3)`, with
full active dimensions 36 through 288. Each is checked at three nonspecial
write times and three nonspecial repair durations, giving 54 terminal points
and both historical labels. Checks cover:

- conservation of both charges and the simultaneous color-flip symmetry;
- exact writer and cut-corner grid matrices, including invariance in the
  larger physical Hilbert space;
- the two color sectors' common corner;
- dynamically acquired pointer accuracy and the complete post-fault vector;
- the full conditional energy law, with degenerate eigenspaces grouped;
- equality of unconditional active-plus-coin states and nonmemory joint
  marginals under the controls;
- mean common-fault work;
- every entry of the two-by-two memory/system outcome table;
- intact, inverted and coin-masked historical scores;
- the short-chain pointer-consumption counterexample and an unflagged-fault
  counterexample.

All listed small-system checks pass the frozen `1e-10` tolerance. The largest
state-vector discrepancy is `4.69e-14`; the largest conditional energy-law
discrepancy is `7.64e-15`; the active-plus-coin density discrepancy is
`2.31e-14`. The physical-to-grid matrix, invariant-subspace, charge,
global-symmetry and shared-corner checks agree exactly in these cases.

The script prints its full result for retention as INDEPENDENT_CHECK.json.
This is independent small-system implementation, not a separate full
physical-Hilbert simulation of the large production reservoirs. The large
grid calculations are evaluated numerically; the uniform gate comes from
the analytic proof, not a claim of interval-certified floating-point output.

## Review of the exact reduction and uniform bound

Initially `Q_B=0`, so S retains its initial ordinary label and B remains
vacuum while the memory is written. The `Q_W=1` writer is exactly the open
path described in PROOF.md. The pointer is dynamically acquired; it is not
initialized directly to the historical label.

After the common SWAP, S is exactly erased and D holds the old label. The
conditional repair sector is the product of the writer and repair paths
with the two origin edges removed. The physical isolated corner is common
to the two color sectors: they are overlapping invariant subspaces, not an
orthogonal direct sum containing two physical copies of that corner. The
separate implementation checks this intersection explicitly.

The half-line estimate follows from integration by parts and the stated
finite sine identity. The finite-path walk bound is conservative for both
the origin and its first neighbor. The exponential/Chernoff tail, the
choices of reservoir lengths and time horizons, and the rational use of
`pi>3.14` and `e<3` are valid.

For the cut-corner comparison, the exact perturbation norm is

\[
 \|V\phi\|^2/g^2
 =|c_Wa_B+a_Wc_B|^2+2|a_Wa_B|^2.
\]

The proof retains its interference term and bounds it by absolute values.
On the first unit interval the largest eigenvalue of the stated matrix is
less than seven; on the remaining interval the square-root and logarithm
bounds are conservative. Duhamel's estimate therefore supplies a valid
state-vector error bound. Projecting onto the joint-failure subspace and
using the triangle inequality yields the rational lower certificate
`0.99156823649...`, above the frozen `0.99` gate throughout the full declared
write/read rectangle. The pointer-acquisition certificate is above
`0.99999`.

These uniform statements concern **a single terminal read at any chosen
time** in the specified windows. They are not joint probabilities of a
continuously monitored history, a never-blank pointer trajectory, or success
under repeated reads; measurements can change the evolution. Joint terminal
retention proves more than transferring the label to S while leaving M blank,
but does not assert an uninterrupted perfectly protected record.

## Why the matched memory intervention works

Let `F_M` exchange the two pointer labels and preserve blank, and let
`G=F_S F_M`. The complete Hamiltonian commutes with G. Immediately after
the common fault, S has support only on its erased state, which `F_S`
preserves. Consequently

\[
 F_M\rho_rF_M=G\rho_rG
\]

on every actual checkpoint state, including its coherent blank component.
Every spectral projector of H commutes with G. Thus the control preserves
the **entire** energy law, conditional on the true historical label and the
coin value, rather than merely matching the mean or variance. This is a
state-support statement; `F_M` need not commute with H on arbitrary states.

Fair averaging of the historical labels makes the active checkpoint state
invariant under `F_M`. With the same maximally mixed fair coin present in
every condition, the complete unconditional active-plus-coin state matches
as well. A local memory unitary also preserves every nonmemory joint marginal.
The control changes correlations with the historical label, not the present
system state or the reservoirs' preparation.

The inversion retains information but reverses the controller's encoding.
The coin mask makes the accessible active state independent of the old
label; the information survives in correlations involving the coin. No global
erasure of information is claimed, and H never reads that key or D.

The fixed terminal observable is

\[
 P_{\rm historical}
 =\sum_{r=0}^1 P_r^D\otimes P_r^S\otimes P_r^M.
\]

The wrong-label sector cannot reach S=r, so its historical scores are zero.
The coin arm is an equal mixture of the intact and wrong arms. Its scores
are exactly one half of intact, giving a uniform joint advantage greater
than `0.495`. These equalities hold for every pair of write/repair times,
not just for favorable numerical samples.

Scoring only whether S is ordinary, or only whether S and M agree, would
miss this distinction. Those symmetry-invariant outcomes can agree across
the arms. Likewise the identical active states necessarily have identical
unconditional future active states under the same H. The different
historical correlations do not contradict that fact: the full states
including D are different.

## Negative controls and limits

**Pointer consumption.** At the registered short-chain times, the full
physical calculation gives historical restoration `1.0000000000000013` and
joint restoration/retention `5.24e-31`. M is blank with probability one at
that read. The extra digits are floating-point rounding; the exact result
is restoration one and joint success zero. This supports the need for the
separate retention gate and exposes the resource cost of the longer writer
reservoir.

**An unflagged color error is not repaired.** An additional adversarial
diagnostic replaces the erasure by a color flip. It leaves `Q_B=0`, so S
stays ordinary, B stays vacuum and the repair interface remains inactive.
The physical test gives historical restoration zero while the probability
of an ordinary state is one. This diagnostic is scored against the known
initial preparation label and is outside the common SWAP-erasure sequence.
It rules out interpreting the result as correction of arbitrary damage.

**Only a classical pointer is copied.** The preferred basis and the erased
symbol are supplied. Agreement of these classical labels does not imply
recovery or duplication of an arbitrary unknown quantum state. No
noncommuting-observable or consciousness claim is tested.

**The rule is supplied.** The two labels are treated equally, and the one
restored is selected by the acquired record. Nevertheless the Hamiltonian
explicitly implements conditional copying and repair. It does not derive a
developmental program, a complex morphology, an emergent controller, a
preferred observer-producing history, or a new GRUT restriction on nature.

**Resources and timing remain explicit.** The blank pointer, both vacuum
chains, pure erased fault register and fair coin are prepared resources.
The SWAP transfers the old S state into D and injects mean work `Delta`.
The displayed Hilbert space does not contain a dynamical clock or work store
implementing that pulse. Pre/post energy-law preservation under the memory
controls is not by itself a constructed finite-time, energy-conserving pulse
controller. The undisturbed phases are autonomous; the full intervention
history is externally timed.

**The demonstrated task is one historical bit and one erasure.** Increasing
reservoir length supplies a finite useful window. It does not establish
scaling to complex patterns, indefinite pointer protection, a repeated-repair
lifetime, or a thermodynamic attractor. A finite closed Hamiltonian remains
recurrent. Reusing the apparatus would require a separate account of ready
memory and reservoir resources.

The earned result is therefore causal use of a dynamically acquired local
classical record in a specified finite repair task, under unusually strong
state and energy matching. Its mathematical success does not upgrade the
program to new fundamental physics, cosmic predetermination or consciousness.
