# SCOUT_0 W2 P-23a RESULT — prediction gate on the amplitude/width-free variance relation

**Candidate (spawned by P-17 §8):** at `O(β)` the odd-in-`a` retained-mean shift is a fixed linear
functional of the measured variance, `(S − D)(t) = −12β ∫₀ᵗ g(t−s) μ(s) v(s) ds`, with no bath
amplitude, width or temperature left once `v` and `g` are measured.
**Gate:** the auditor's four questions, plus the canonical gate's R1 (derived), R4 (non-absorbable)
and R5 (reach) from `GRUT_PREDICTION_GATE_GAMMA_T.md` §8. Script: `p23_invariant_checks.py` (A1).

| Gate question | Answer |
|---|---|
| Does it follow from GRUT-specific certified structure? | **No.** It is the generic small-noise expansion. For **any** drift `f` with additive noise, `Δc₂ = T f″(a)` (A1, symbolic in `f`), so `S − D = ½∫ g f″(μ) v ds` at first order. The quartic `f″ = −24βx` gives the record's `−24βTa`. The only GRUT input is the declared drift (S-9, supplied). |
| Or from generic perturbation theory for any additively forced nonlinear system? | **Yes:** standard second-order (van Kampen / moment-closure) expansion. |
| Can a non-GRUT model satisfy it? | **Yes, every additively forced smooth model does**, with its own `f″` (A1: double well `x − x³` gives `−3∫gμv`). By P-17/P-18 it is also satisfied identically by exogenous noise, a harmonic bath and a chaotic driver. |
| What observation could violate GRUT while satisfying the alternatives? | **None.** A violation would falsify "additive forcing + smooth drift", which is shared by every alternative in the class, not GRUT. |
| R1 / R4 / R5 | R1 fails (rests on the supplied drift, not on derived structure); R4 fails (absorbable into any additive-noise model); R5 not reached |

**Verdict: REDISCOVERED-KNOWN / NON-DISTINCTIVE — NOT A GRUT PREDICTION. Retired.**
