# Handoff: q-theory vacuum construction (state on 22 October 2026)

Claude's usage limit stopped this pass partway through. The four build tracks and three of their verifiers finished.
These four agents failed at the limit and never ran:
- the dissipation verifier (`tracks/v_dissipation/` may be partial);
- the relation agent;
- the physics review;
- the scope review.

Nothing below has had the final adversarial review. Treat it as **draft, unreviewed**. The records are
`tracks/*/NOTES.md` and `tracks/v_*/VERIFY.md`; read BRIEF.md and BRIEF_STAGE2.md first.

## The construction (QV1, borrowed from Klinkhamer–Volovik and credited to them)

**Action.**
- S = ∫√−g [½K(q)R − ε(q)] + S_m, with q the 4-form field strength.
- Field equations: (q) ε′ − ½K′R = μ, a global constant; (g) KG = T − (ε − μq)g − (g□ − ∇∇)K.

**State condition.** P_q: μ = μ₀[ε_∞], the equilibrium of the final, fully renormalized vacuum. This is a final-state
condition on the state. The dynamics do not enforce it. The tuning moves into μ (dμ₀/dc = 1/q₀), which Klinkhamer and
Volovik concede.

**Verification of the equations.**
- The KV papers, read in full from TeX.
- A blind derivation.
- An exact Bianchi-identity check.

## Results so far

1. **Solved example** (`tracks/scenario`, verified in `tracks/v_scenario`: 25 claims confirmed, 6 corrected, none of
   which changes a verdict).
   - **What the example shows:**
     - Under P_q every constant cancels, including Standard Model loops, thresholds and a heavy constant c.
     - The q-mode tracks adiabatically through the electroweak and QCD crossovers and leaves no constant remnant.
     - A Planck-scale vacuum relaxes and, with the scalaron decay source, reheats to about 5×10¹⁶ GeV.
   - **Checks:** the reduction error is 0.03(H/m)², validated at m/H = 10–10⁴, and the constraint holds to ≤10⁻¹³.
   - **QV1 fails retrospectively.** It predicts essentially zero vacuum energy today (4.7×10⁻¹⁶⁹ GeV⁴), giving
     100θ* = 0.943, which is −327σ; it also fails SN and BAO.
   - **QV1-cal** inserts ρ_Λ through a calibrated offset of μ and reproduces ΛCDM. The calibration supplies the
     residual.
2. **Bound** (`tracks/bound`, verified).
   - For an unscreened q-mode, the non-calibrated dark energy today is Ω_dyn = (9/4)Ω_m²(H₀/m_q)², and it is stiff.
   - Eöt-Wash at α = 1/3 gives m_q ≥ 3.65 meV, so Ω_q ≤ 3.5×10⁻⁶².
   - A light q also overproduces oscillation energy (Ω_q h² from 10² to 10⁶; scenario §4).
   - The screened escape is the published disappearing-Λ f(R) class (Starobinsky; Hu–Sawicki). Its final state is
     de Sitter, which conflicts with P_q's final-state reading.
3. **Multiple vacua** (`tracks/vacua`, verified: 20 confirmed, 9 corrected).
   - "The final vacuum" is ambiguous once several vacua exist.
   - Under the natural readings (the global minimum, or our vacuum with a nonzero decay rate), P_q requires the
     electroweak vacuum to be **absolutely stable** (or degenerate within ρ_Λ). That turns into a Standard Model
     relation, the stability bound on m_t(m_h, α_s).
   - Measured masses lie on the wrong side at about 3.1σ (verifier value; the track had 3.3σ). This is retrospective.
   - The relation assumes the SM holds up to the Planck scale; new physics below the instability scale is the escape.
   - Prior art on vacuum stability is credited. The reading is an added assumption. The prediction is **not yet
     frozen**, because the relation agent did not run.
4. **Dissipation** (`tracks/dissipation`, **unverified**).
   - KV's dissipative variants form μ only in an empty, vacuum-dominated universe, after about 10¹²³ e-folds.
   - With a radiation era and BBN, μ still has to be prepared, to about 10⁻⁵⁵ or better.

## Status in Codex's terms (proposed, not yet reviewed)

QV1 is the strongest partially specified construction, with a solved example. It is consistent, but it does not
produce the observed residual and fails late-universe data retrospectively. It has no GRUT-native element. GRUT's
contribution is organizational: the P_q formulation and the protect-or-form bookkeeping.

**Single next mathematical dependency (proposed).** Derive the present offset δμ = μ − μ₀[ε_∞] (required: −ρ_Λ/q₀)
from dynamics or from a stated state rule, in a universe with a radiation era. Without that derivation the residual
remains a calibration.

## To do after the limit resets

1. Rerun the failed agents: resume workflow `wf_fd579cde-b38` from its script; completed agents replay from cache.
2. Freeze the vacua m_t relation in `tracks/relate/results/frozen_predictions.json` before any new m_t data (FCC-ee
   threshold scan) is used.
3. Write CONSTRUCTION.md with the reviewer corrections, then regenerate MANIFEST.
4. Add an AGENT_COORDINATION entry and reply to Codex comments 6324148604 and 6324495415 on issue #13.

Candidate 1 status is unchanged (Codex): "Candidate 1 is not independently test-ready; universal equivalence to cold
dark matter remains unproved."
