# Independent reviewer and architecture handoff

The next reviewer should be a different researcher, not another run of
this same author's verifier. All results in this packet remain authored
within the Builder program. No external reviewer has been contacted.

## Smallest proof review

1. **Charge blocks and vacuum:** reconstruct F1 from H_0 and H_-1 without
   imposing a photon cutoff. Try to find a missing sector, boundary term,
   or misuse of a creator norm. Check the conserved-charge energy offset.
2. **Kernel:** rederive the boundary spectral density, two endpoint
   estimates, closed-walk order and finite tail through 8N^2. A failure
   of the O(N^2)-horizon comparison breaks the claimed formation theorem.
3. **Memory/current:** retain both [S,A]a^* and [S^*,A]a. Check the
   equal-time operator ordering and the powers of N in all three errors.
   Inspect the damping integral and macro-failure threshold independently.
4. **Transport:** construct the common space explicitly as in
   TRANSPORT_REVIEW.md. Seek a physical short word that differs before
   q, including storage inconsistency, field motion and return paths.
   Check both tails, the boundary force sign, position domains, moments
   and the discrete Paley-Zygmund event.
5. **Composition:** use the actual conditional receiver, every accepted
   word, original roots, failed trials and frozen time sample. Recording
   and injection must not quietly replace the prepared state.

A counterexample or failed estimate should identify exactly which S9
claims must be withdrawn. Passing CHECKS.json cannot overrule one.

## Bounded new coupling result

Review M1-M9 against C as frozen at `a0c53ef`. Covariance means simultaneous
relabeling of disorder, root and events, not invariance of each fixed
inhomogeneous realization. The root flag is an existing distinguished
register; M1 does not require it. Check the norm, fixed-time event,
Duhamel coefficient, prefix-loss count and acceptance multiplier.

If a further constraint excludes M1/M2, identify it explicitly and say
whether WM2 already contains it for b or it is a new primitive. A smaller
class with a newly supplied bond selection rule does not refute the
counterexample to C. Likewise, reproducing a's occupied-bond hopping for
b does not establish universality by a new mechanism.

The conclusion needed for the architecture decision is already bounded:
current symmetries alone do not select shared matter geometry. Preserve
the construction; do not open another solver or gravity campaign to
avoid that result. A successor would need an actual common-coupling
principle and tests of its permitted direct and virtual violations.

## Receiver lemma review

Check R1's factors and all vanishing errors before using R2. Check that
the bandwidth/center family has uniform constants and positive overlap.
The flat-receiver obstruction uses energy expectation, not a bath
approximation. Neither example reduces the supplied mode-count scaling.

## Reproduction and provenance

From this directory:

```sh
OPENBLAS_NUM_THREADS=1 python verify.py --run-checks
```

Default `python verify.py` checks hashes, frozen source pins and stored
claim/status consistency only. `--run-checks` executes fresh finite
diagnostics and exact rational checks, while keeping CHECKS.json archived.
No option launches the closed S9 workers or a physical-L formation solver.
Float calculations are ordinary double precision with asserted tolerances;
they are not interval enclosures. Integer/rational checks are exact.
