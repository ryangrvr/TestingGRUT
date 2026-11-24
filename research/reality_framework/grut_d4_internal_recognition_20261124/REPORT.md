# GRUT-D4 internal recognition and historical identity maintenance

24 November 2026. One time-independent Hamiltonian writes two copies of
a classical historical bit. It then responds to discrepancies in the
three current copies and restores the historical configuration after
unflagged substitutions. When up to four single faults are separated by
the stated waits, the original system configuration and both records are
jointly correct with probability above **99.51%** at the terminal read.
The probability of all three being correct with the expected cumulative
counts in both reservoirs exceeds **99.19%**. Both bounds apply to the
actual Hamiltonian with the writer still present.

D3's fixed repair dynamics required a designated erasure symbol. D4 puts
an ordinary wrong label in the system or in either memory. The repair
interaction receives only the current physical degrees of freedom. The
two parity comparisons distinguish a single system error from either
single memory error. They can detect more general disagreement, but
historical identification depends on the fault class. [D3 report](https://github.com/ryangrvr/TestingGRUT/blob/0f18c546ab0b0ad2f0f5968c52bc02349f2c01e4/research/reality_framework/grut_d3_record_repair_20261124/REPORT.md)

| Criterion | D4 verdict |
|---|---|
| Historical restoration | Proved for either historical label and any sequence of up to four separated single-site faults. The score requires the historical bit in S and both memories. |
| Internal fault detection | Two present-state parity comparisons activate the corresponding exchanges without a fault flag. The interpretation of a syndrome as a culprit holds for the single-fault class. |
| Repeated self-maintenance | Proved within the declared finite horizon with both reservoirs initialized once. Same-site reuse and memory faults are included. |
| Autonomous timing | Repair responds under the always-on H. Preparation and read windows and the label-blind fault pulses use laboratory time. An internal clock or autonomous source of faults remains open. |
| Physical novelty | A finite realization of classical repetition-code correction under standard quantum mechanics. The parity and majority interactions are supplied. |

## Model and bound

Beginning with blank memories, a resonant writing exchange copies the
current classical source bit into both memories and excites a writer spin.
The writing success exceeds 99.9996% in the declared phase window. On a
written triple, let `r` be its current majority and let the two syndromes
be `x_i=S xor M_i`. The Hamiltonian exchanges a syndrome excitation
with the corresponding spin chain, holding the current majority fixed.
The interaction treats the two historical colors equally. The physical
site geometry specifies the comparisons with S.
[MODEL.md](MODEL.md)

Each drain has 2304 initially vacuum spins. The encoded readout is proved
for any inter-fault and terminal wait in `[128/g,144/g]`, and for a writer
phase within `1/512` of `pi/2`. Writing is slow because the supplied
coupling is `g/65536`. The waiting bound comes from finite-path return
amplitudes and a norm estimate for any residual excitation seen by the
next fault. The proof also bounds the writer's reverse exchange. There
are no resets or intermediate reads. Each bound concerns one read after
the chosen number of faults. [PROOF.md](PROOF.md)

| Separated faults | Historical codeword lower bound | Codeword and exact counts lower bound |
|---:|---:|---:|
| 1 | 0.99986510 | 0.99986510 |
| 2 | 0.99924396 | 0.99911457 |
| 3 | 0.99778066 | 0.99680031 |
| 4 | 0.99518102 | 0.99194364 |

The numerical comparator retains the actual physical fault maps and
omits only the writer exchange after acquisition. It uses exact finite
path modes and coherent logical-sector recombination. It tests 243
four-fault timing rows covering all 81 ordered site sequences at each
of three waits. Its lowest joint score is 0.99999917797. The conservative
physical bounds above include the effect of the writer and apply
between sampled waits. The numerical score measures the comparator.
[RESULTS.json](RESULTS.json)

## What the system can infer

The central obstruction is an exact collision of present states. History
0 with one S flip and history 1 with two memory flips give the same
physical triple, writer and reservoirs. With identical accessible states
and equal priors, any physical recovery operation has at most 50% success
across those histories. The two-copy obstruction has the same structure
after just one allowed substitution. A mismatch shows that copies differ.
Knowing which of them represents the past requires redundancy, a fault
assumption or additional evidence.

For d equiprobable historical alternatives, repair success is bounded by
the optimal probability of guessing the historical label from the full
accessible state. If all available history information lies in a
D-dimensional register, success is at most `min(1,D/d)`. For example,
one bit alone supports at most 25% for eight alternatives. Exact repair of
the no-fault and all single-substitution inputs for a fixed three-copy
codeword requires `1+3(d-1)` distinguishable environment outputs. Those
are separate constraints on historical information and error storage.
[INFORMATION_BOUND.md](INFORMATION_BOUND.md)

## Failed controls and resource costs

The two-site simultaneous faults activate a discrepancy and produce the
incorrect unanimous word with probability above 99.9865% at a registered
read, by the one-fault theorem applied to that majority. All three flips
leave the opposite unanimous word
and no parity flag. Both types have exactly zero joint historical success
after nominal writing. With vacuum reservoirs the undamaged control has
zero repair-bath activation and zero source color change. It retains a
small writer reversion probability, about 0.00772% over the longest tested
wait. A prepared excited bath can instead create a discrepancy, and the
one-spin reservoir recurrence returns a fault with probability one at
about `6.2831853/g`. [CONTROLS.md](CONTROLS.md)

The first S fault requires work 8g and either first memory fault requires
4g. The larger cost for S reflects two parity excitations, independent
of historical color. Four repeated S faults leave between 7.9355 and 8
absorbed quanta in expectation by the physical bound. The total net pulse
work is bounded conservatively by about 15.8711g to 48.00002g. The pure
reservoirs and blank memories provide information storage. No dynamical
work source and clock implement the test pulses in this model.
[CONTROLS.md](CONTROLS.md)

## Review and reproduction

The design was frozen before production in [commit 5c59a2d](https://github.com/ryangrvr/TestingGRUT/commit/5c59a2dce5858752a785f563dd779fa1034a9371).
The protocol and model retain their original hashes. The separate AI
review constructed four small full physical Hamiltonians with the writer
present. It also compared nine scenarios to a dense ready-sector oracle.
The maximum joint probability discrepancy in that comparison was
`4.22e-15`. All 243 timing rows and 4,648 stored numbers replayed exactly
in the recorded runtime. The separate review ran within Work.
Standalone external RAI and human review remain open.
[CRITIC.md](CRITIC.md) · [VERIFICATION.json](VERIFICATION.json)

The earned result is an internal discrepancy response coupled to bounded
historical maintenance, with its reliance on redundancy, vacuum and fault
spacing explicit. In D4 the comparisons and recovery rule are part of the
Hamiltonian. Autonomous quantum correction by engineered interactions
has established precedents. Deriving a new physical law or an emergent
recognition rule remains a separate research problem.
[Cohen and Mirrahimi](https://arxiv.org/abs/1409.6759)
