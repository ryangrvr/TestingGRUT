# Skeptic review of track cf_vacuum (P3: closure as vacuum equilibrium)

27 October 2026. Claude, adversarial skeptic `cf_vacuum_skeptic`, closure-forcing parallel pass (issue #13 context).
Reviewed: `NOTES.md`, `results/*.json` and `code/` of this track, against `../../BRIEF.md`, FORMAL_CORE §3–4, the vacuum
SYNTHESIS, the q-construction HANDOFF and track `vacua`, Codex's FOUNDATION_GATE (2908e91) and exact architecture
(04c53f5: THEORY.md, ACTION.md, THEORY_DATA.json).

**Labels.** [D] derived here; [C] computed by `skeptic/code/` (output in `skeptic/results/`); [Q] quoted from a source
read by me (read level in `skeptic/results/sk_sources.json`); [E] estimate; [A] argued; [R] recalled.

**Scope.** No validation, no exclusion of the real world, no new theory. Every comparison with an existing
measurement is **retrospective**. This is same-session AI review, not outside replication.

**Reproduction.** `cd skeptic/code && python3 -B sk_line.py && python3 -B sk_checks.py && python3 -B sk_summary.py`
(about 3 min). The code imports nothing from `../../code`. The RGEs were transcribed again by hand from the
Buttazzo et al. TeX. The weak-scale matching is either Buttazzo 2013 ("B13") or Bednyakov et al. 2015 ("B15"); B15
has not been used anywhere in the program before. The track's line uses SMDR.

---------------------------------------------------------------------------------------------------------------------

## 0. Outcome

The track's computations survive. Every number I recomputed agrees within the stated errors, by a different method,
code or matching. Its physics conclusions also survive: P3a + SM is a strict, testable, non-GRUT-native restriction,
GRUT has no internal escape at the reference point, and the Σ corner is tuned and invisible.

Two verdict labels do not survive:
- **P3a + GRUT reference** is labelled INCONSISTENT WITH GRUT BRANCH. This is refuted.
- **MIXED** is used for P3a + SM and for P3d. It is not one of the brief's five verdicts.

I also found five issues the track missed or under-stated. The two that bear most on the claims are:
- **A one-sided UV systematic: the sign of the Planck-suppressed |H|⁶ coefficient.** GRUT's own L_higher contains this
  coefficient, and it is unstated.
- **The α_s category.** α_s = 0.1221, the z = 0 value, is PDG's own electroweak-subfield average.

| Postulate | Refuted? | What I recomputed | Track verdict → should be |
|---|---|---|---|
| P3a-SM | no (substance) | line 170.96 (B15) / 170.90 (B13) vs 171.00 (SMDR); z = 3.54 / 3.67 vs 3.45; slopes, μ*, Λ_λ agree | MIXED → **STANDARD** (relabelled published absolute-stability hypothesis; testable; retrospectively disfavoured 1.6–3.7σ) |
| P3a-GRUTref | **yes (verdict)** | portal needs 0.118 / 0.176 (ratios 23.5 / 24.8 to the invisible-width bounds); hidden bound 9.05e-4 | INCONSISTENT → **STANDARD** (inherits P3a-SM; adds nothing measurable; the conflict is with retrospective data, not with the branch) |
| P3a-GRUT-Σ-corner | no | λ_ΣY needed 0.0469 / 0.0569; Y-vacuum 2.2e98 ρ_Λ deeper; CW −1.369e-5; θM_S identity; valley identity | UNTESTABLE (agree) |
| P3d-MPP | no (substance) | same line (two-sided) | MIXED → **STANDARD** (relabelled Froggatt–Nielsen / Volovik); the equality is UV-sensitive downward (§3.3) |
| P3v-volume | no | — | STANDARD (agree) |

## 1. Data re-checked at the primary source [Q]

| Input | Track | Primary source today | Status |
|---|---|---|---|
| M_h | 125.13 ± 0.11 | PDG 2026 S126M: 125.13 ± 0.11 (S = 1.5) | confirmed |
| M_t direct | 172.60 ± 0.27 | PDG 2026 Q007TP: 172.60 ± 0.27, OUR AVERAGE (S = 1.5) | confirmed |
| M_t σ(tt̄) | 172.1 ± 0.6 | PDG 2026 Q007TP4: 172.1 ± 0.6 (8 inputs) | confirmed |
| α_s(M_Z) | 0.1180 ± 0.0009 | PDG 2026 QCD review (dated 1 June 2026): 0.1180 ± 0.0009; FLAG 2024 0.1183 ± 0.0007 (eq. 9.23) | confirmed |
| α_s, electroweak subfield | not used | Table 9.1: **0.1221 ± 0.0027** | see §5, N2 |
| arXiv IDs credited | 20 | all 20 resolve to the stated authors and titles (arXiv API) | confirmed |

## 2. Independent recomputation of the SM line [C]

**Method.**
- **RGEs.** Full 3-loop for g1, g2, g3, y_t and λ (4-loop pure QCD in g3), transcribed by hand from 1307.3536 App. A.
- **Matching.**
  - B13: Buttazzo's fits at μ = M_t.
  - B15: Bednyakov et al. Table 1 at μ = 173.21 GeV. This includes 2-loop EW, O(αα_s) and 4-loop QCD.
- **Criterion.** min_{μ ≤ M_Pl} λ(μ) = 0, as in the track.

**Validation.**
- λ(M_Pl) = −0.01436 and y_t(M_Pl) = 0.3825 at Buttazzo's central point; the paper gives −0.0143 and 0.3825.
- B15 at Bednyakov's inputs gives **171.441**; the paper gives 171.44.
- B13 at Buttazzo's inputs gives 171.07. The paper's formula gives 171.115 and its quoted value is 171.53. The quoted
  value does not follow from the paper's own matching and RGEs either. This confirms the track's (and vacua's)
  observation.

**PDG 2026 results.**

| Quantity | SMDR (track) | B15 (this review) | B13 (this review) |
|---|---|---|---|
| M_t,crit (125.13, 0.1180), loop 3 | 170.997 | **170.964** | 170.902 |
| dM_t/dM_h | 0.494 | 0.4955 | 0.4950 |
| per 0.0009 in α_s | 0.355 | 0.352 | 0.353 |
| loop 2 − loop 3 | −0.110 | −0.108 | −0.109 |
| μ* (GeV) | 4.96e17 | 4.66e17 | 4.59e17 |
| Λ_λ at M_t = 172.60 | 4.9e10 | 4.5e10 | 3.5e10 |
| Λ_λ at ±1σ (172.33 / 172.87) | 1.4e11 / 2.0e10 | 1.2e11 / 1.8e10 | 9.5e10 / 1.5e10 |
| Λ_λ at σ(tt̄) mass 172.1 | 3.9e11 | 3.4e11 | 2.5e11 |
| λ(M̄_Pl) at 172.60 | V = −7.5e70 GeV⁴ (λ_eff, e^{4Γ}) | −0.0105 → RG-improved tree V = −9.2e70 GeV⁴ (3.7e117 ρ_Λ) | −0.0109 |

The three matchings span 0.095 GeV, inside the track's ±0.113 GeV. Bednyakov's own scale-variation error on
M_t^cri is asymmetric, +0.17/−0.36 GeV [Q]; using the upward 0.17 gives z = 3.42 instead of 3.54. **SMDR gives the
highest line of the three, so the track's z values are the most conservative.**

**Tension (track's σ recipe), B15 line, retrospective.**

| Top-mass input | z |
|---|---|
| direct | 3.54 |
| direct, +0.3 GeV MC | 2.97 |
| direct, +0.6 GeV MC | 2.16 |
| ATLAS+CMS | 3.12 |
| σ(tt̄) | 1.61 |
| direct, FLAG α_s | 3.75 |
| direct, **electroweak-subfield α_s 0.1221 ± 0.0027** | **0.03** |

With the B13 line every z is about 0.1 higher. The track's numbers are confirmed.

## 3. Is the consequence forced by the stated postulate?

### 3.1 Λ_ours = F(V_ours − V_true) [D]
With K′(q) = 0 the q-equation fixes the same q in every vacuum, so Λ_i − Λ_j = V_i − V_j exactly. F is then the
identity, and the track's "≥" is an equality. With K′ ≠ 0, F is only monotone with F(0) = 0. That is enough for
ΔV ≈ 3e117 ρ_Λ. The step holds.

### 3.2 The strict restriction rests entirely on the clause "the residual offset may not absorb vacuum-energy differences" [A]
QV1 without calibration predicts Λ_ours = 0, which is retrospectively excluded (HANDOFF). QV1-cal adds an offset δμ.
- If δμ were a free calibration, it would absorb ΔV, and S_C = S_base: no restriction.
- P3a therefore needs δμ to be of order ρ_Λ, i.e. P_q holds to ρ_Λ precision. That is the natural reading of
  "approximately P_q", and the track states it in the postulate, so nothing is smuggled.
- But the residual's origin is the program's open "single next dependency" (HANDOFF). The dissipation track finds
  that μ must be prepared to about 1e-55.

So the restriction is conditional on a residual mechanism that is not yet derived. **It is not derived from P_q
alone.** The track should say this next to the gate-item-3 derivation.

### 3.3 "SM valid to M_Pl" contains a sign condition on Planck-suppressed operators (new, computed) [C][E]
GRUT's action explicitly contains L_higher, "higher-dimensional operators, with a stated cutoff and matching
coefficients" (THEORY.md §2). THEORY_DATA lists "higher-dimensional coefficients" as **missing**. I added
c₆|H|⁶/M̄_Pl² to the RG-improved potential and redid the P3a criterion, V ≥ 0 up to the cut. Line shifts relative to
c₆ = 0, B15 matching:

| c₆ | −1 | −0.1 | −0.01 | −0.001 | +0.01 | +0.1 | +1 | +4π |
|---|---|---|---|---|---|---|---|---|
| cut at M̄_Pl = 2.435e18 | **no M_t in [150, 180] satisfies P3a** | −8.27 GeV | −0.78 | −0.07 | +0.008 | +0.020 | +0.041 | **+0.075** |
| cut at M_Pl = 1.22e19 | none | none | none | −1.96 | +0.008 | +0.020 | +0.041 | +0.075 |

Consequences:
- **For falsification the line is robust.** No perturbative positive c₆ raises it by more than 0.075 GeV. To stabilise
  the PDG central point you would need c₆ ≈ 1.8e12. The track's statements "z > 5 ⇒ P3a + SM fails" and "the escape
  must sit below Λ_λ" are therefore confirmed.
- **For confirmation the line is not robust.** A negative c₆ of order −0.01 (or −0.001 with the cut at M_Pl) lowers
  the line by about 1–2 GeV or more. So P3a also forces **c₆ ≳ 0** in GRUT's L_higher, and the result depends on the
  cutoff convention. A future M_t below 171 GeV would therefore **not** confirm P3a. The MPP equality (P3d) is defined
  only with c₆ = 0.
- **The error budget should say so.** The frozen relation's ±0.11 ± 0.03 GeV should state that it assumes c₆ = 0
  (more precisely c₆ ≥ 0) and a cut at M_Pl. "Still supplied" (NOTES §7) should list the sign of GRUT's Planck-suppressed
  Higgs operators explicitly.
- **Credit.** Branchina–Messina 1307.5193 make the analogous point for the lifetime [Q abstract].

### 3.4 Other forcing steps checked
- **Worldline reading reduces to (a) when Γ > 0.** Confirmed as a fixed-point argument (vacua §3.3). It is a measure
  choice; the track labels it [A].
- **GRUT's fields cannot stabilise at reference.** My own RG extension uses 1-loop portal pieces derived symbolically
  from ½Tr M⁴. It gives b ≥ 0.118 and λ_HX ≥ 0.176 (track: 0.117 and 0.175), with invisible-width bounds 5.03e-3 and
  7.11e-3, i.e. factors 23.5 and 24.8. Σ at reference lies above the instability scale. Y has no vev. **Confirmed.**
- **The KV citation for "global".** KV 1102.3152 (TeX, lines 653–667) does say μ₀ "corresponds to the true quantum
  vacuum of a global equilibrium (the final Minkowski state…)". Their contrasted "false vacuum" is the **thermal**
  intermediate phase above a transition, not a zero-temperature metastable vacuum. Applying it to EW metastability is
  the program's extension. The credit stands, labelled [A].

## 4. Standard physics or a relabelled result?

- **The relation and the hypothesis are published.**
  - The relation is the standard absolute-stability bound: Degrassi 1205.6497, Buttazzo 1307.3536, Bednyakov
    1507.08833, Hiller 2401.08811.
  - "Our vacuum is the global minimum, with the true vacuum at Λ = 0" is published: AHHKM astro-ph/0005111, KV.
  - Hiller et al. already state the precision needed for a 5σ verdict (M_t error of 200–300 MeV) [Q].
  - So the **testable content of P3a + SM is identical to the published absolute-stability test**. A z > 5 outcome
    would test that hypothesis and says nothing about closure.
- **The combination is not in the literature.** My arXiv queries ("q-theory" + Higgs; Klinkhamer/Volovik + Higgs
  vacuum) found no paper that combines q-theory's P_q with Higgs absolute stability. This agrees with vacua §4. The
  combination is bookkeeping, as the track says.
- **The verdict word.** Under the brief, "relabelled published result" maps to STANDARD. "MIXED" is outside the
  five-verdict set. The track's prose ("not GRUT-native", "best *testable* candidate") is accurate, but the label
  should be **STANDARD (relabel; testable; retrospectively disfavoured)**. "Best testable candidate" is acceptable
  only with the qualifier that its test is the published stability test.

## 5. Verdict strength

**P3a-GRUTref is refuted as a verdict.** INCONSISTENT WITH GRUT BRANCH would need the GRUT branch to require
something P3a forbids. It does not:
- GRUT does not fix (M_t, M_h, α_s).
- The hidden sector meets its own P3a bound at reference (margin about 1100, reproduced).
- Candidate 1's tree potential has X = 0 as its minimum.
- If the SM is on the stable side, the heavy window's Higgs-instability obligation (track §5.4) **disappears**.

P3a + GRUT reference is therefore internally consistent for every M_t ≤ M_t,crit. The only conflict is with measured
inputs, at 1.6–3.7σ retrospectively. That is neither an inconsistency with the branch nor an exclusion, which the
brief forbids claiming. Correct label: **STANDARD (inherits P3a-SM; GRUT adds no measurable content and no internal
escape at reference)**. Keep the conditional: "if metastability is established, P3a + GRUT reference fails."

**P3a-Σ-corner: UNTESTABLE stands.** All its numbers were reproduced by a different method.
- **Hidden-sector bound.** A brute-force 2-D grid in (|Σ|, |Y|) at cos = −1, with 3-variable refinement, gives
  λ_Σ ≥ 9.051e-4.
- **λ_ΣY needed** in the window: 0.0469 at λ_Σ = 2.8e-5 and 0.0569 at 1.8e-10.
- **Depth.** At λ_Σ = 2.8e-5 and λ_ΣY = 0, the Y-condensed vacuum is 2.2e98 ρ_Λ deeper.
- **CW coefficient.** Exactly −μ⁴/(24π²M_Y⁴) = −1.369e-5 (sympy series).
- **Mixing identity.** θ·M_S = √(2δ)v is reproduced by exact diagonalisation to ≤ 1e-4 relative.
- **Valley-running identity.** Verified symbolically.

The corner's mass edge inherits the matching dependence of Λ_λ, 3.5–4.9e10 GeV. That is a soft factor and does not
change the verdict.

## 6. New issues and missed items

1. **N1. Planck-suppressed operator sign** (§3.3). P3a forces c₆ ≳ 0 in GRUT's L_higher, and P3a cannot be confirmed
   by masses. The upward systematic is ≤ +0.075 GeV and the downward one is unbounded.
2. **N2. α_s category (retrospective).**
   - The z = 0 value α_s = 0.1221 equals PDG 2026's electroweak-subfield average, 0.1221 ± 0.0027 (Table 9.1). With
     it, z ≈ 0.03 (B15) or 0.08 (B13).
   - The headline should keep the world average. The note should say that the tension is carried by the non-EW α_s
     categories together with the MC-template top mass, and that z spans 0 to 3.75 across published α_s categories.
     A decisive test needs both inputs improved, as the frozen relation already says.
3. **N3. The GRUT reference hidden quartics are non-perturbative before M_Pl.**
   - With one-loop 16π²β = 20λ² for λ_Σ = 1 at h and λ_Y = 1 at M_Y, λ reaches 4π at 7.6e15 GeV and 5.7e15 GeV. The
     Landau poles are at 1.4e16 and 1.1e16 GeV [C][E]. This is a 1-loop estimate; λ_ΣY stays about 0 because only
     e⁴ ≈ 1e-28 terms feed it.
   - So "GRUT valid to M_Pl", the premise of P3a + GRUT, is false at the reference.
   - The track's "hidden-sector P3a satisfied, margin 1100" holds only for field values ≲ 1e16 GeV. Beyond that, P3a
     needs a hidden-sector UV completion with no deeper vacua.
   - The Higgs-direction conclusions are unaffected, because the sectors are decoupled at reference.
4. **N4. The ξ_H obligation has an upper end.**
   - The track gives ξ_H ≳ 0.05–0.07 (2011.03763).
   - Curvature-induced Higgs production at the end of inflation and during preheating destabilises the vacuum for
     ξ ≳ 1 (Herranen–Markkanen–Nurmi–Rajantie 1506.04065; Ema–Mukaida–Nakayama 1602.00483) [Q abstracts].
   - The heavy-window obligation is therefore a window, roughly 0.06 ≲ ξ_H ≲ O(1), and depends on the reheating
     model. It is independent of P_q, as the track says.
   - Supporting numbers: the per-patch Gaussian probability of |h| > h_max is about 0.10 against e^180 patches, which
     confirms that the obligation is real.
5. **N5. Theory error and conservativeness.** The three matchings span 0.095 GeV, and Bednyakov's asymmetric
   +0.17/−0.36 GeV is larger upward than ±0.113. The track's SMDR line is the highest of the three, so z = 3.45 is the
   conservative choice. No change is needed, but the asymmetric literature error should be quoted.
6. **N6. Attribution of the "brief's formula".**
   - "171.53 + 0.15ΔM_h − 0.34Δα_s" does not appear in `BRIEF.md` on disk. It must have come from the computed task
     text, and the track should cite it that way.
   - The formula's M_h slope is also wrong, by a factor of 3.3 (0.15 against 0.495), not only the sign of its α_s term.
7. **N7. Missing row: reading (b) with uncalibrated QV1** (vacua §3.2). It gives the same inequality through Γ = 0, but
   QV1 with Λ_ours = 0 already fails retrospectively (−327σ; HANDOFF). It has no separate SM content and needs a
   one-line entry for completeness.
8. **N8. Held-out FCC-ee reach (illustrative).**
   - z = 12.1 assumes that the conversion from the threshold mass to M_t(pole) cancels the pole-mass renormalon
     ambiguity, about 110 MeV (1605.03609 [Q abstract]).
   - If it does not cancel, z ≈ 9.3. With a 0.25 GeV ambiguity z ≈ 5.7 [E]. With the σ(tt̄)-like central value and
     0.11 GeV, z ≈ 6.4.
   - Recommendation: also freeze the line in a short-distance quantity, such as m_t(m_t) or y_t(M_t).

## 7. Files

- `skeptic/code/sk_rge.py`: own 3-loop SM RGE and the B13/B15 matchings.
- `skeptic/code/sk_line.py` → `skeptic/results/sk_line.json`: validation, line, tension, Λ_λ, c₆ scan.
- `skeptic/code/sk_checks.py` → `skeptic/results/sk_checks.json`: hidden-sector bound, CW, mixing and valley
  identities, portal loop stabilisation, invisible width, Landau poles, inflation numbers.
- `skeptic/code/sk_summary.py` → `skeptic/results/sk_summary.json`: key numbers and verdict checks.
- `skeptic/results/sk_sources.json`: primary sources, read levels, searches.
- Downloads (the PDG 2026 QCD review PDF and the KV 1102.3152 TeX) are kept outside the repository, in the scratchpad
  folder `closure/cf_vacuum_skeptic/dl/`.
