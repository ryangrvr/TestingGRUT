# Short review route

1. Check the source and comparator pins against B2-SS ac4bf51 and B2-SO
   3700986. Do not treat their inherited proofs as externally validated.
2. Re-derive DECOMPOSITION.md (1)-(6), especially the unnormalized squared
   joint probability, high-measure versus high-state cost, and cut identity.
3. Check the exact energy/type cross terms in (9)-(12).
4. Prove or break the orthogonality-aware capture lemma (13). Verify that
   its conditional 24/25 capture threshold implies the original 17/50
   condition, and that capture itself has NOT been established.
5. Check dimensionless Fourier scaling and all BV endpoint terms in
   (17)-(20). None of the regularity or singular-weight hypotheses is
   claimed for the physical h.
6. Independently sum SCALAR_LIMITATION.md's rational weights and variance.
   Verify its continuous early bound and exact recurrence at 128*pi.
   Confirm that it is NOT the actual source or a physical counterexample.
7. Run controls.py, then verify.py after restoring the archived
   CONTROLS.json. Their agreement checks implementation consistency only.

The all-time inference limitation has an elementary exact proof. The
actual source-to-receiver and late fuel estimates are still unresolved.
One assistant wrote proofs and code; no independent review, external
validation, empirical test or new physical regime is claimed.
