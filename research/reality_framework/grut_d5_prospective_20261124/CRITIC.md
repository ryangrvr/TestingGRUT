# Separate critical review of GRUT-D5

25 November 2026. This review uses an independent implementation of the
reversible local head/tape rules and their reachable physical Hamiltonian.
It is separate internal AI review, not external RAI execution, independent
experimental replication, human expert review or peer review.

## Verdict

**The scoped mathematical construction passes this review.** A supplied
finite local reversible program acquires a classical record, grows an
initially incomplete pattern, reconstructs it with authenticated damaged
symbols in its archives, and maintains the completed word for a bounded
clock interval. Its head is a physical dynamical stage variable; no runtime
gate sequencer drives its internal scheduled operations. Laboratory tape
flips are also treated without measuring the head or discarding clock tails.

The record intervention is causal in the main implementation. The matched
pilot null is equally important: it acquires the same record and produces
the same intact growth and repair while never using that record for writes.
Thus appearance plus an accurate correlated record is insufficient evidence
of functional record use. The intervention distinguishes the two laws.

**This is programmed development, not an independently emerging goal.** The
target family, source phase, rule bit, route, blank archives, visit schedule
and transport weights are supplied. The original pilot h remains physically
available, including when M is inverted or disabled. The experiments show
which input the main law uses; they do not establish information-theoretic
necessity of M, anticipation of an unavailable future choice, consciousness,
or new fundamental physics.

## Independent physical-graph oracle

Run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 oracle.py
```

The oracle imports none of the builder's implementation. Its basis states
contain the head position, all nine head registers and every displayed and
archive qutrit. It implements the actual reversible local hops and their
inverses from MODEL.md. For laboratory disturbances it closes the reachable
physical graph under both Hamiltonian hops and the tape pulse, rather than
restricting propagation to the original legal history.

The implementation tests N=8 with short routes L=64 and 128. These are
implementation references, not members of the frozen L=256N success family.
The L=64 pulse-closed graphs contain 780 physical basis states. Direct sparse
exponential propagation is compared to the grouped endpoint-amplitude
formula in 27 external main/null cases. Eighteen additional direct graph
propagations test the internal controls, including inactive-coin branches.

All tests pass the frozen `1e-10` implementation tolerance. The worst
endpoint-amplitude discrepancy is `7.87e-15`; target-plus-archive probability
discrepancy is `2.56e-15`; canonical stage probability discrepancy is
`6.67e-16`; normalization error is `2.89e-15`. The internal endpoint/control
error is below `1.34e-15`. The main/null external success discrepancy is zero.
The clock spectrum and initial spectral-law discrepancies are below
`2.00e-15` and `5.24e-15`, respectively.

The oracle also verifies the assay, both unfinished-prefix stages, the last
rebuild as the first completed ordered word, reverse-pass parity alignment,
archive authentication, every main control for both coins and every seed,
and the null's invariance under M controls. Complete intact data histories
of main and null agree exactly, including the acquired record and archives.

The deliberately short L=64 third-pulse test has canonical completed-stage
probability about `0.10954` and joint target/archive success about `0.55143`.
These values are retained. They demonstrate that the oracle actually includes
early/late clock components and endpoint misses. They are not misses of the
frozen large-route threshold, which uses a much longer padded clock.

This is an independent finite reachable-graph calculation. It does not
materialize the complete exponential many-cell Hilbert space, test fabrication
of the graph, or provide interval arithmetic for floating-point solvers.

## Review of the clock and pulse proof

Each local update is a reversible permutation. The swap-based rebuild and
maintenance moves retain the overwritten symbol rather than erasing it.
Reversing the order of operations on backward hops is essential; the oracle
checks these inverses. The maintenance pre-toggle gives the correct parity
at the first reverse site and avoids an alternating-pattern phase error.

The supplied nearest-neighbor head exchange has operator norm at most J
per edge. The graph geometry and head alphabet are fixed local resources,
with sizes independent of N per site. On the full one-head work sector,
the history-prefix gauge is unitary and gives

\[
 H=D(\omega J_x\otimes I)D^\dagger.
\]

The binomial head distribution, perfect endpoint transfer, energy spectrum,
input spectral law and recurrence follow from this identity. It remains valid
for work states outside the intended legal history, which is necessary after
an unscheduled tape pulse. It is not an argument confined to one preferred
classical trajectory.

For a data pulse at tau, the exact endpoint vector has one common phase and
the contributions

\[
 \sum_j p_j(\tau)
 V_LV_j^\dagger X V_j|s\rangle.
\]

Therefore identical **complete final basis states must be grouped before
squaring**. A group of total binomial weight q contributes q squared to its
unconditional endpoint probability. Treating those weights as a classical
mixture would give the wrong head-arrival and success probabilities. The
independent physical propagation tests this coherence directly.

The constant-work-band comparison is also valid: the state-vector error is
at most `2*sqrt(1-q)`, giving failure at most `4(1-q)`. The three declared
bands, phase windows and padded lengths satisfy the stated mean-to-boundary
inequalities. Hoeffding's bound and the conservative rational exponential
envelope justify the strict large-route thresholds without rounding a tiny
failure probability to zero. The subsequent data-persistence bound follows
by adding the late-stage failure amplitude once; recurrence still ends the
useful interval.

The initial whole-energy law is the same across internal selectors, data
labels and the null. This follows from the common clock spectral measure at
head position zero. It is stronger than matching a mean energy. A laboratory
pulse is different: its mean work is zero for the real-permutation/chiral
preparation, but its energy distribution and variance can change. In the
short oracle the third-pulse variance changes from approximately `0.06059`
to `0.36794`. Zero mean work does not make that intervention resource free.

## What prospective indistinguishability does and does not prohibit

Let F label a proposed future target and let `rho_f` describe the complete
present physical information available before the scored action. Include
the seed, clock, reservoirs, archives, controllers and correlations, not
just the displayed tape. If all later processing before that action is a
common channel Lambda receiving no new F-dependent input, and mutually
exclusive whole-target projectors `Q_f` satisfy `sum Q_f<=I`, then

\[
 P_{\rm success}
 =\sum_f\pi_f\operatorname{Tr}[Q_f\Lambda(\rho_f)]
 \le P_{\rm guess}(F\mid\text{present information}).
\]

Pulling the outcome effects back through Lambda produces a POVM on the
present state, proving the inequality. For two equally likely targets the
right side is `1/2 + ||rho_0-rho_1||_1/4`. Identical conditional states give
only the largest prior. Four independent equally likely whole words therefore
give at most `1/4` success before their label becomes available.

This does not prohibit predicting a future event correlated with current
records, or responding after a new condition arrives. It also does not apply
the `1/d` baseline to overlapping rewards. For identical present states and
score effects `G_f`, the independent-target baseline is

\[
 \lambda_{\max}\!\left(\sum_f\pi_fG_f\right).
\]

For example, with phase known, a fixed uniform word already has average
per-cell agreement `3/4` against equally likely uniform/alternating targets.
Its exact whole-word success is only `1/2`. The model correctly scores
whole words.

A later trace-preserving operation cannot change the probability of an
earlier instrument outcome. Conditioning instead on a selected later
outcome can change that earlier conditional distribution; it is postselection
and must include its selection probability. A future timestamp alone is
not evidence of information arriving from the future.

D5's final family member is already specified by the initial physical h,G
and the supplied rule. Its successful unfolding is fully compatible with
these bounds. The data do not test an independent inaccessible future target
or a law selecting observer-producing histories.

## The matched null and the limits of the memory claim

The null is a distinct, explicitly specified Hamiltonian with the same
route, clock, registers, assay and perturbations, whose deposits read h
instead of M. Because intact histories have M=h at every deposit, their
entire states agree, not merely their visible shape. A tape-only pulse
does not break equality of those two head values, so the intact null also
matches the damaged-tape trials.

An M intervention changes the main model's subsequent colors and leaves
the null's colors unchanged. Inverting G affects both. These are valid
functional tests of which input the transition law consults. They do not
show that the acquired M contains information absent elsewhere: h remains,
and disabling by a reversible swap moves the old value into a backup.
The control coin also retains information in correlations. No globally
irreversible forgetting or unique internal self-model has been demonstrated.

The comparison must therefore preserve two conclusions simultaneously:
the acquired record is causally used in the main law, and the same growth,
record acquisition and repair appearance occurs in a law that does not use
that record. Appearance alone does not establish a mechanism of recognition.

## Supplied resources and interpretation boundaries

1. **The route is a program.** Its assay, growth tranches, fault edges, return
   route, rebuild and maintenance passes are selected physical interactions.
   The deposition rule is compact and local, but the route layout and
   engineered weights scale with N. No emergent selection of this program
   or of the uniform/alternating family is proved.
2. **Regeneration is not discrepancy diagnosis.** The sweeps overwrite correct
   cells as well as damaged ones. They demonstrate regeneration using stored
   input and fresh archives, rather than D4's discrepancy-triggered response.
   Internally scheduled diagnostic faults are part of the supplied circuit;
   the unmeasured laboratory pulses are the separate unscheduled tests.
3. **Clock position is coherent.** Coordinate snapshots describe conditional
   history states, not an observed head moving monotonically. Backward
   amplitude can undo and redo scheduled gates. The endpoint net program is
   exact, and one-time laboratory stage probabilities are bounded, but no
   joint observed milestone history or strict absence of all early completed
   amplitude follows.
4. **The word expands a short description.** Four possible length-N patterns
   carry only two bits of family choice. The acquired M is one correlated
   bit. Increasing N does not create N independent target bits or derive
   complex biological development.
5. **Prepared storage is substantial.** The head's internal alphabet, localized
   one-head preparation, blank tape and two blank archive banks are supplied.
   Arbitrary reversible overwrite of an N-qutrit input while preparing one
   fixed output requires at least `3^N` distinguishable environment records.
   The archival banks supply capacity explicitly. Their reset or indefinite
   reuse is not part of this result.
6. **The clock has a cost.** Its initial spectral law, spatial weights and
   precision are assumptions. The chosen normalization bounds bond strengths
   while transfer time grows linearly with N. A clock and route are now
   physical model components, rather than an external gate sequencer, but
   laboratory disturbance/readout timing and their actuators remain supplied.
7. **Closure does not derive purpose.** The internal dynamics is unitary and
   finite, with entropy transferred into correlations and archives and a
   full recurrence. It supplies neither consciousness, unknown-quantum-state
   copying, a new physical law, nor a consciously planned Big Bang.

The earned result is an exact finite example of programmed development and
regeneration, with a physical stage head, causal record controls, authenticated
damage and a null that rules out interpreting appearance as evidence of
record use. Its limitations remain part of the result.
