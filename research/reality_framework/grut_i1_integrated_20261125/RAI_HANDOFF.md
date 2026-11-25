# External RAI replication packet for GRUT-I1

**Status.** Prepared for an external reviewer. The internal Work review and
its reproduction are recorded separately. No work has been dispatched to
a separate RAI session.

Use TestingGRUT, parent draft PR #20 at
cecc5f99202f77c55255b528372bf3d53e12a149. Retrieve the I1 freeze from
cee2d2197b700b1932cdba3db758bcbecd3b864e. Check that the supplied model
and grid match the preproduction files. Compute independent results before
reading the archived result files when possible.
[Protocol](REALITY_TEST_PROTOCOL.md) · [Frozen grid](PARAMETERS.json)

## Required checks

1. Reconstruct D7's fixed Bell measurement and access interface from its
   proof. Reproduce weak shallow \((1+\sqrt{5/8})/2\) and weak deep \(0.75\)
   alongside the optimized blind comparators. Identify every subsystem
   unavailable to the downstream recovery and distinguish the optimized
   comparisons from fixed-memory interventions.
   [D7 parent proof](../grut_d7_functional_recovery_20261125/PROOF.md)
2. Independently exponentiate the two-site pulses. Evolve \(B,M,Q,R\),
   with explicit conditional states for the two initial source signs.
   Trace out the correct bath cells and preserve real joint evolution on
   repeated contacts. Compare Bell score, population, phase, record
   contrast and both angular-error conditions with exact idle. Distinguish
   contrast to the initial and terminal source.
3. Check binding hypotheses for the accessible-programme optimum.
   In particular, try an added source bypass, early encoder or access to R
   only as an explicitly changed interface. Compare the original interface
   against ordinary optimal quantum recovery.
4. Recompute the hot and cold energy and entropy ledgers. Verify the
   state dimensions and preparation for every contact, including unused
   cells. In the reused case check the fault in replacing coherent reuse
   with a new thermal state. Check when the continuous generator is valid
   and what caps remain on the runtime and contact resolution.
5. Compare the independent and common two-site repair arms under both
   registered matches. Separate the density score of large independent
   arrays from an interacting-system phase transition. Recompute one
   network family outside the ring and exhibit what geometry can be
   recovered from the selected source and receiver alone.
6. Challenge identifiability and novelty. Show where the resource-limited
   value-of-information comparison requires a broader set of possible
   policies, where a useful score follows established physics, and what
   observation would actually exclude the conventional comparator. Check
   that the broad statements about fundamental principles respect the
   preserved locality failure, the separate Candidate 1 EFT and the
   absence of a cosmological construction.

## Expected report

Return identifiers for immutable inputs, independent executable steps,
numerical precision, convergence and stopping conditions, deviations
from the frozen gates, equations establishing the strongest upheld
claims, and one A–D verdict for I1. Support novelty with primary
precedents or with a quantitatively testable distinction. Label an
independent script replay, an independent code implementation,
an independently derived proof and hardware replication separately.

The laboratory prediction is a proposed transmon test. Report a circuit
result and its access and calibration uncertainty before drawing any
physical conclusion. The preregistered interval includes both contrasts,
and frozen systematic allowances remain separate from shot noise.

[Review scope](CRITICAL_REVIEW.md) · [Reproduction](REPRODUCTION.md)
