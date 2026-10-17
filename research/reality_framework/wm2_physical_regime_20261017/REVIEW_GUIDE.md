# Central mathematical review — no research-history reconstruction needed

The S9 claims are Builder derivations with execution checks, not external
validation. The smallest useful review consists of the following steps.

1. In FORMATION_ANALYSIS.md (2)-(4), check the complete-word algebra and
   ambient charge-sector extension. Cross-type append/pop products vanish,
   while two append operators generally do not commute. This is not the
   failed radial or subset reduction.
2. Check the endpoint integration-by-parts estimate (6) and finite closed-
   walk tail (7) on the **physical** L=64N^2 chain through 8N^2.
3. Check the Heisenberg identity (9) and differentiation of
   `[A_(u-v)](t-v)`. Both creation and absorption commutators must remain.
   Contact occupation bounds the latter and the vacuum +1 in the former
   must be paid. Verify the history-retardation error (11)-(14).
4. Check that the instantaneous product S^dagger A_u is the weighted free-
   event projector sum. Then the current comparison is exact up to its
   bounded remainder; it is not an assumed Markov master equation.
5. Check the negative-damping Gronwall step and Markov inequality for the
   quantum free-count observable. Rational tests verify the constants and
   N_0=2^28 endpoint; they do not replace steps 1-4.
6. In BALLISTIC_TRANSPORT.md, check the common-Hilbert-space comparison
   and equality of all short interaction-picture words. Both complete
   tails are paid. The half-line is only a reference, not a physical
   truncation or a supplied operational metric for the original system.
7. Check the boundary-acceleration, moment and Paley-Zygmund bounds, then
   subtract escaped-prefix and storage-error events in the integrated
   original-trial probability.
8. Separately check certificate.py's direct residual/unknown-Gram bound
   and structural batch guards. No square sums across reversals are allowed.
   Its primary numerical pass/fail cannot validate the analytic all-N proof.

Disproof of any of these steps requires withdrawing the corresponding
claim; a passing saved verifier or an AI-authored audit cannot overrule it.
All source, initial freeze and failed numerical endpoints remain archived.
