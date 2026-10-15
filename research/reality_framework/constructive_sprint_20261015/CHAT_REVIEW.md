# Review handoff: constructive GRUT working theory

This packet continues a2f447d through published parent1801afe without modifying
the frozen CP-12 Hamiltonian. Read MANUSCRIPT.md and EXECUTIVE_REPORT.md first.
The objective is to assess a constructive, conditional matter–geometry mechanism,
not certify a completed theory of reality.

## New claims worth independently challenging

1. **Actual graph-choice access.** The original neutral ray has nonzero exact
   charged-circulation coefficients at N6,N7. At N7 the two exhaustive four-regular
   host families have unequal per-host coefficients, ratio0.5979550201916068.
   Check the principal-minor determinant reduction, signed path/cycle covers,
   multiplicity4 for distinct charged vertices and neutral normalization.
   `charged_access/run.py` includes16 independent integer half-leg checks.
2. **Induced interactions without added laws.** Virtual matter generates
   `+J_l(W_l+W_l†)`, with
   `J_l=l Catalan_(l-1) t^(2l)/(2m)^(2l-1)`, and a positive coordination-variance
   term in the controlled static heavy-matter limit. Separate supplied zero
   incidence angles from optimized holonomies. At frozen parameters use the
   energy-dependent Feshbach result, including the kappa*Q energy shift and full
   residual; do not treat sixth-order perturbation at t=m=1 as controlled.
3. **Exact conditional locality.** For a single graph sector,
   `||P_out H P_g v||²=(Gamma/M)²(M-Q) sum_e ||r_e v||²`.
   The state error is bounded by `Gamma sqrt(Q(M-Q))/M * |T|`.
   Keep the escape diagonal in `P_g H P_g`. The proof does not permit replacing
   a coherent multigraph neutral state by an average of these bounds. Unbounded
   onsite rotors are handled before restricting the local even observable algebra.
4. **A finite neutral-origin recorded chain.** N4,Q4, wait5, lag0.5 produces
   the charged-circulating event with weight0.0550221, recorded host persistence
   near1 and a positive density susceptibility about0.000131169. Verify the
   complete graph/event channel, physical-basis expansion after recording, fixed
   initial-record receiver, conditional denominator, and record energy/entropy
   costs. The reported final circulation allows either orientation and any final
   vertex occupation. This is external diagnostic recording, not autonomous selection.

## Exact and numerical controls

Run `python verify.py --rerun` with BLAS/OMP threads set to1. Review the nested
derivations, source hashes, two electric cutoffs, omitted-space estimates and
independent integer controls. The previous1801afe reproduction has now completed
fresh, without a cache: all65 scientific comparisons pass. Its old immutable
pending status is preserved as history.

The source pulse convention is `exp(-i theta B)`, so
`chi=-i <[O_H,B]>=2 Im <U psi,O U B psi>`. Pinching and selective analysis must
not silently discard the event's weight or energy cost. The code uses ordinary
floating point plus stated analytical bounds and numerical allowances; it is not
a formal interval-arithmetic certificate. The worker review is internal.

## The gap and next experiment

None of these claims selects three dimensions, proves a macroscopic causal cone,
derives Einstein equations or chiral Standard Model physics, explains the Born
rule, or eliminates external Hamiltonian time. The neutral formation state is
coherent and externally prepared. A fixed bounded-degree host is supplied in the
locality theorem. The finite recorded prototype is rare and explicitly measured.

The next discriminating test is growing-size **unfiltered neutral evolution**,
measuring charged-network probability, relational persistence, intrinsic ball
growth and spatial density response in the same state, compared with random
and frozen-host controls. The new formation kernel and locality bound should
reduce that calculation. Do not count another isolated small-system nonzero
response as proof of autonomous spacetime formation.
