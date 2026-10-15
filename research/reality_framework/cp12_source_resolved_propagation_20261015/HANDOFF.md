# GRUT-RAI handoff: source-resolved locality challenge

Parent: `4e7d6c5a91fe96737ab9a8a140ec2b50eb140ebc`.
All changes are isolated to this research packet. Canonical theory,
Stage3, PR12 and the frozen CP-12 law are unchanged.

## Result requiring a decision

**The recorded-ring candidate metric fails the specified six-site test.**
Under the complete frozen Hamiltonian, an initial direct chord carries
a much larger local density response across ring distance 3 than the
response to a receiver at ring distance 2. At waiting time 5, lag 1/4:

- `chi_03` is `[0.00201985496,0.00202060216]`.
- `chi_02` is `[2.85250e-6,3.65290e-6]`.
- The absolute response ratio is `[552.946,708.362]`.
- The same original ring remains with probability
  `[0.904057487,0.904058719]`.
- The gauge-invariant density-phase source has variance
  `[0.230599074,0.230599117]`, so it is demonstrably active.

The lag-1/8 ordering also fails, while its smaller receiver response is
unresolved about zero. The five-site equal-ring-distance control gives
a chord-endpoint/nonendpoint response ratio `[1683.50,1701.37]` at 1/4.
These are signed susceptibilities, not transfer probabilities. Pointwise
nonmonotonicity is not by itself a universal no-go theorem for local
quantum theories; the conclusion concerns this recorded metric and its
active bypass, not every possible CP-12 phase.

## Formation and scope

The matched sequence uses **Q=N+1**, unchanged coefficients and original
equal coherent neutral preparation. It is auxiliary, not the primary
Q=2N scaling sequence. Exact spanning-ring formation at time 5 is
8.025916%, 1.753043%, 0.344661% for N=4,5,6. The N=6 longest-chord class
occurs with probability 0.115830%; with the fixed source on its chord,
0.038610%. One labelled endpoint core-host outcome has probability
`6.434999e-6`. No free favorable state was supplied.

The original `N=5,Q=10` sector is also run from its neutral ray at waiting
time 2. Its exact-ring probability is approximately 0.126724%; at lag 1/4
the ring-distance-2 response is resolved and negative, while the near
response is positive. Maximum distance there is only two. This packet
does not settle large-size Q=2N propagation.

The six-site complete core-host instrument includes all 540 successful
records plus failure, has energy change about 0.00151714 and record
entropy about 0.06444334 bits. Full labelled N=4 controls confirm that
nonselective recording itself changes later response. Final-core joint
and conditional responses are reported separately; only raw local
receiver responses form the primary comparison.

## New computational capability

The N=6 neutral matrix has 1,980,546 invariant orbits and 133,348 complete
electric-boundary rows. A real 128-column Krylov reduction, with explicit
boundary and relation residuals, gives full-state error below
`2.162e-9` at waiting time 5. A larger 144-column calculation agrees
with its state to about `2.01e-14`, within both bounds. The source basis
has 513,066 orbits, every outgoing hopping/rewiring row retained for
error control, and positive source variance. Full labelled N=4 response
and instrument checks pass. Ordinary floating point is disclosed;
these are not outward-rounded formal certificates.

## What has and has not been earned

The tests do **not** demonstrate Milestone A: a viable emergent spatial
phase. Both prevalence and the candidate transport metric are problems.
They do not establish Milestone B, a relativistic/gravitational continuum,
or Milestone C, a distinctive empirical prediction. Quantum dynamics,
external time, fermionic statistics, gauge structure and coefficients
remain inputs. Ring ball growth is selected by the event definition;
it is not a derivation of spatial dimension.

The recommended next constructive assignment is a separately specified
metric/transport architecture with prospectively declared receiver tests
and held-out sizes, or a fresh full-active-host geometry hypothesis.
Do not protect the ring metric by relabelling the observed chord after
the response is known. `ARCHITECTURE_CHARTER.md` states the required
formation, invariance and universality checks, including why simply
blocking all neutral-state hopping prevents geometry formation.

The five research gates in `PROGRAM_TESTS.md` adopt the reviewer's latest
acceptance standard. Consciousness/observer mechanics supplies no
exception to a failed propagation calculation. Preserve the mathematics
and error-control methods; require the architecture to earn locality.
