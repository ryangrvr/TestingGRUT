# Separate critical review of GRUT-D4

24 November 2026. This is review by a separately tasked AI in this Work run.
The reviewer implemented independent physical tensor Hamiltonians. It is not
external RAI execution, external replication, human expert review or peer
review.

## Verdict

**The registered finite mathematical claim passes this review.** One fixed
Hamiltonian writes two classical replicas, responds to unflagged discrepancies
on any of the three sites, and preserves the historical codeword with the
stated terminal probability through up to four separated faults. No bath
reset, intermediate read, fault flag, or historical-label input to the
Hamiltonian is used. The analytical bound includes the writer that remains
present during correction.

This is stronger than D3's specified erasure repair. It is nevertheless a
supplied majority/parity correction rule under a restricted fault class.
Three replicas do not determine an arbitrary past from their current state.
The successful score depends on the historical majority surviving each
separated fault, with residual-discrepancy errors explicitly bounded. The
identical-state and simultaneous-fault counterexamples establish the limit.

The model does not derive its correction rule, pointer basis, preparation,
geometry, timing or reservoir resources from closure. It establishes neither
conscious recognition, a new GRUT law, nor a cosmological selection principle.

## Independent physical implementation

Run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 independent_check.py
```

The review code constructs S, both complete qutrit memories, W and every
reservoir spin. Its four full physical cases have chain lengths
`(1,1),(1,2),(2,1),(2,2)` and dimensions 144, 288, 288 and 576. Both the
actual writer-on H and the comparator omitting only its exchange are built.
Dense matrix eigendecomposition evolves repeated actual physical pulses.

A separate three-physical-bit construction with lengths `(3,3)` has dimension
512 and tests the ready-sector Gram/Fock routine. It uses direct physical
basis amplitudes rather than the builder's time-mode reduction. Nine pulse
scenarios include four repeated/alternating four-fault sequences, all three
single faults, a simultaneous pair and a simultaneous triple. The complete
joint distributions include logical majority, both parities and both bath
counts. Different count-axis sizes are padded and compared in full.

The largest ready joint-probability discrepancy is `4.22e-15`; the largest
normalization discrepancy is `1.34e-15`. The ready-coordinate permutation,
physical fault maps, charge commutators, color symmetry, logical conservation,
writer subspace, ready embedding and identical-state collision all agree
exactly in the tested matrices. All are below the frozen `1e-10` tolerance.

The actual writer-on and exchange-omitted vectors differ by at most
`2.34e-5` across the small repeated-fault tests, below the corresponding
`kappa*T=6.8817138671875e-5` bound. An acquisition-window residual is also
tested. Exact writing is evaluated with the direct two-state Rabi formula
after independently checking its invariant physical subspace. This avoids
accumulating large-time eigensolver phase roundoff at the very slow nominal
writing time.

The script prints a full JSON result for archival as INDEPENDENT_CHECK.json.
It validates the implementation independently on these small systems. It is
not a full physical-matrix simulation of the large reservoirs or an
interval-arithmetic certificate for floating-point results.

## Audit of the actual writer-on proof

The supplied charges commute with the complete H. The initial writing
subspace is exactly two dimensional for each label; its on-site energy is
constant and its exchange is `kappa*sigma_x`. The phase-window acquisition
error and the `1/512` state-distance bound are valid.

The ready basis permutation is bijective. Its current-majority coordinate
is not a hidden historical register. Substitution into the inverse map
gives the displayed physical fault maps, including the possibility of
logical-majority changes on residual discrepancies. That possibility must
not be replaced by a logical coordinate frozen through every fault.

The half-line return estimate, the conservative finite-boundary walk count,
and the exponential tail bound are valid. The on-site Delta phase does not
change the return modulus, but is retained in the exact multiple-pulse
calculation. Moving an annihilation operator through the endpoint Majoranas
using the CAR yields the stated harmonic waiting bound without resetting a
previously used chain.

For the syndrome-only K comparator the logical value is a fixed basis
label. On the component where the true F changes that label, the F and K
outputs have orthogonal logical labels; their difference has norm
`sqrt(2)` times the component norm. Thus the universal `2*a_ell` pulse
comparison bound is valid. The writer exchange has operator norm kappa,
and its interval errors accumulate to at most `kappa*T` under the same
unitary pulses. Adding acquisition, pulse and terminal-source errors gives
the registered rational physical-codeword bound above `0.995`.

The expected-charge projectors commute with propagation and source
occupations. The charge-error hybrid estimate and the final commuting-event
union bound therefore also justify the joint codeword/count certificate.
The proof displays a conservative lower bound `0.9919236501...`; the
executable rational square-root envelope can tighten this slightly. Neither
uses rounded numerical samples for acceptance.

In the exact ready solver, historical scoring coherently recombines the two
logical-X sectors. In particular an S pulse in the negative sector includes
the phase of `(X1 Z1)(X2 Z2)=-Y1 Y2`. The physical matrix comparisons test
this phase-sensitive rule, repeated source reuse, and simultaneous pulses.

These are single-terminal-read guarantees for every allowed wait and fault
sequence, not a continuously monitored no-error history, a sequence of
successful measured corrections, or an unlimited lifetime. The large ready
numbers omit the writer exchange as a mathematical comparator. They must
not be quoted as direct probabilities of the actual writer-on model; the
analytical physical error budget is essential.

## Information-theoretic limit: a response cannot identify an absent fact

Let `sigma_h` be the complete accessible present state for historical label
h, including all reservoirs and any clock used by the repair interaction.
Any fixed-H evolution and terminal readout is a quantum channel followed by
a measurement. Historical restoration cannot outperform optimal
discrimination of that ensemble:

\[
 P_{\rm restore}\le P_{\rm guess}(\{\pi_h,\sigma_h\}).
\]

For two equally likely labels the latter is

\[
 \tfrac12\bigl(1+\tfrac12\|\sigma_0-\sigma_1\|_1\bigr).
\]

Identical present states therefore permit average success at most one half.
No derived response, probe, susceptibility or subsequent internal dynamics
can recover a distinction absent from its complete input state. If a model
does distinguish such histories, it contains additional accessible
information, for example a fault-dependent reservoir, clock, controller or
reference. That information may be physically legitimate, but must be
identified and counted. Equality only of local marginals or energy
distributions would not suffice for this no-go premise.

### Two replicas

A pair `(a,b)` with distinct labels can result from historical a with the
second replica changed, or historical b with the first changed. With common
accessible ancillary preparation these histories have identical current
inputs. A deterministic decoder fails on at least one; a randomized or
quantum decoder has worst-case success at most `1/2`. An asymmetric fault
prior can favor trusting one replica, but does not identify the true culprit
from the pair alone.

This obstruction concerns **unflagged substitutions**. Two replicas can
repair a single located erasure because its flag supplies the missing
information. D3's erasure setting is therefore not contradicted.

### Three or more replicas and simultaneous faults

For an n-fold repetition code over d labels, let `f_r` be the number of
positions equal to r. Under at most t substitutions, label r is compatible
with the current string precisely when `f_r>=n-t`. Unique correction for
every allowed history is guaranteed only when `n>2t`.

If k labels are compatible with the same observed state, worst-case
historical success is at most `1/k`. For `t<n`, the maximum ambiguity is

\[
 k_{\max}=\min\left(d,\left\lfloor\frac n{n-t}\right\rfloor\right).
\]

For three binary replicas, `001` can mean historical zero with one error or
historical one with two. For `d>=3`, a three-distinct-symbol string has three
possible historical labels under two errors, yielding a `1/3` worst-case
bound. A majority identifies a culprit only under its stated fault model.
Correlated simultaneous faults are not excluded by a single-fault theorem.

D4 provides an exact stronger physical witness: nominally written historical
zero with `X_S`, and historical one with `X_M1 X_M2`, both give the same full
state `100` with W ready and vacuum reservoirs. The independent code checks
this collision and its energy law exactly. The subsequent histories cannot
both be restored correctly. Three flips produce the opposite unanimous
codeword with zero parity discrepancy. The chosen H preserves that wrong
majority rather than recovering the inaccessible past.

## Syndrome storage and alphabet scaling

For each fixed historical r, a three-replica d-label single-substitution
task has `K=1+3(d-1)=3d-2` orthogonal no-fault/single-fault inputs. If one
common unitary maps all to the same pure codeword from the same pure
environment preparation, the final environment must retain K orthogonal
records. Thus `dim(E)>=K` and its necessary capacity is `log2(K)` bits.

For an initial environmental state of rank `r_E`, the exact support bound
strengthens to `dim(E)>=K*r_E`. Uniformly mixing those K inputs also gives
the entropy-capacity condition

\[
 \log_2\dim(E)-S(E_0)\ge\log_2K.
\]

These are necessary exact-correction bounds, not sharp resource estimates
for the approximate D4 certificate.

| Alphabet d | Physical ready dimension d^3 | Required exact syndrome dimension at least |
|---:|---:|---:|
| 2 | 8 | 4 |
| 3 | 27 | 7 |
| 4 | 64 | 10 |
| 8 | 512 | 22 |
| 16 | 4096 | 46 |
| 64 | 262144 | 190 |

The two parity drains have capacity for four orthogonal no-fault/single-fault
records after transfer: bath counts `(0,0),(1,0),(0,1),(1,1)`. This satisfies
that basic binary requirement. Their cumulative counts do not generally
retain fault order. No larger-
alphabet Hamiltonian or probability certificate is established by the table.

This accounting also exposes a possible false construction: several
degenerate orthogonal errors coupled identically to one common bath endpoint
leave dark source combinations. Only their bright combination drains. A
three-error rank-one interface would permit at most `1/3` correction from
each equal-weight basis error unless further mixing or channels are supplied.
D4 avoids that particular failure through its explicit two-syndrome
structure; it does not erase syndrome information by a many-to-one unitary.

## Physical negative controls and resource limits

1. **Finite recurrence is exact.** The one-spin-bath M1-fault sector is the
   three-state path with couplings g and kappa. Its matrix and invariant
   physical subspace agree exactly. The defect returns with probability
   `1.0000000000000004` numerically at `2*pi/sqrt(g^2+kappa^2)`; the exact
   probability is one.
2. **A loaded drain can generate a discrepancy.** In that same physical
   sector, an occupied bath gives a fault probability
   `0.999999999767...` at the tested time, versus approximately zero for the
   vacuum bath. These are different energy/preparation resources, not a
   matched causal comparison. Reusing the large drains within the finite
   theorem is justified by its return bound, not by assuming they remain
   empty or harmless indefinitely.
3. **Fault work is not uniform.** From a perfect initial codeword, the tested
   work is 8 for S and 4 for either memory in the registered units. S changes
   two parities. Later pulses may encounter occupied sources and have a
   different work law. Cross-location experiments are not claimed to have
   matched energy costs.
4. **Prepared capacity is substantial.** The 4608 pure drain spins, both blank
   memories, W and the slow `kappa/g=1/65536` writer are supplied. The long
   nominal preparation time and its phase precision are also resources.
   This is not a minimality result or derivation of those scales.
5. **The history is not fully autonomous.** H writes and responds without a
   correction trigger, but the tested faults and final read are externally
   timed ideal operations. W is not a demonstrated internal clock. No
   dynamical pulse source or work store is included, and global unitary
   entropy preservation alone does not close that implementation gap.
6. **The alphabet is classical and supplied.** Three copies of the selected
   basis labels and their correction do not constitute three clones of an
   arbitrary unknown quantum state, protection against all quantum errors,
   or self-generated semantic identity.

The earned result is an explicit standard-physics realization of acquired
redundancy, internal discrepancy response and finite repeated historical
repair, with honest ambiguity and resource bounds. A new physical law,
emergent recognition rule, independent experimental validation, consciousness
and cosmic predetermination remain separate, unestablished claims.
