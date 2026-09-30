# DP-R0-01 HOSTILE AUDIT

> Independent hostile review of the DP-R0-01 execution (parameter-space reconstruction).
> The audit attacks the result in BOTH directions per the dispatch: falsify the lower
> bound, falsify the upper boundary, attack the window splice, attack the selector
> conclusion, attack the prediction firewall. Canonical GRUT-RAI untouched.

---

## 0. Audit scope and inputs

Inputs (read-only): `DP_R0_01_EXECUTION.json`, `DP_R0_01_EXECUTION.md`,
`DP_R0_01_PROTOCOL.md` (frozen), `DP_R0_01_MODEL_REGISTRY.json`, primary sources
re-inspected for contested claims. The five malformed comparative-foundations JSON
cases were NOT touched (separate repair process).

## 1. Attack 1 — the lower bound (Gran Sasso, 4 Å)

| Question | Finding |
|---|---|
| Same DP variant? | The Donadi 2021 experiment targets the **Markovian white-noise DP with nucleon mass-density source**, i.e. the DP-2007/finite-R₀ family. The execution's applicability field ("all finite-R0 Markovian DP variants") is correct with one refinement: the bound is derived under the **standard nucleon mass-density composition assumption**; a variant with genuinely different primitive mass content (e.g. effective density spreading into extended degrees of freedom beyond the stated model) is not directly bounded by this channel. |
| Radiation channel → R₀ mapping model-specific? | YES — the bound comes from the collapse-induced diffusion rate of charged nucleons in Germanium, translated through the DP diffusion coefficient. The mapping is variant-specific: it constrains the **charge-associated diffusion** which in DP tracks the mass-density localization rate. A variant in which collapse acts on different operator content (e.g. correlates with different degrees of freedom) would need its own translation. **Classification of the bound as variant-specific HARD_EMPIRICAL is CORRECT, but "applies to all finite-R0 Markovian variants" should carry the explicit qualifier: "within the standard mass-density-coupled formulation."** |
| "Hard" justified? | YES for the tested formulation: the measurement is a dedicated underground radiation search with independently characterized backgrounds; the bound is data-driven, not criterion-based. The pre-2021 bounds (heating, GW, neutron-star) were correctly superseded. **No downgrade.** |
| Assumption transfer between variants? | One flag: DP-1989 (unregularized) is marked "effectively excluded by any finite lower bound" — this is correct in the sense that its collapse rate diverges as R₀ → 0, but the execution should note explicitly that the Donadi bound excludes it **by implication through the R₀ family**, not by direct measurement of the 1989 form. Recorded as a scope note, not a downgrade. |

## 2. Attack 2 — the upper boundary (effectiveness criterion)

| Question | Finding |
|---|---|
| Is 10⁻⁴ m a constraint? | **NO — correctly classified as EFFECTIVENESS_CRITERION.** The source's own language confirms: "Clearly there is a lot of ambiguity in this request, given that the macroscopic domain cannot be defined precisely." The number follows from choosing (a) a specific object (graphene plate of particular L, d), (b) a perception timescale τ_obs ~ 0.01 s, (c) a definition of "smallest visible object." Each choice shifts the bound. |
| Does the cited literature establish a numerical upper bound? | It establishes a **regime motivation**, not a bound. Different object choices in the same paper yield different values. The execution's warning ("to some degree arbitrary") is accurate and **should be retained verbatim** in any downstream use. |
| Applies to every finite-R₀ implementation? | The criterion is formulation-level (it asks whether the model guarantees collapse of macroscopic superpositions at all), so it applies broadly to the DP family — but as a **requirement on the model's purpose**, not as data. Confirmed classification. |

**Audit conclusion on upper boundary: the execution correctly refused to treat 10⁻⁴ m
as data. The conditional character is not cosmetic; the upper end of the window is
a modeling preference.**

## 3. Attack 3 — the six-order window as a single class

The window splices:

- a hard empirical **lower** bound from one formulation (mass-density-coupled Markovian DP);
- a criterion-based **upper** bound from a classicality-requirement analysis;
- across variants that share the 2007 Markovian core but differ in the 2013 conjectural repair.

**Verdict: the splice is legitimate but only as a "family-level" statement.** The window
is coherent **iff** (a) all members share the mass-density collapse operator (required by
the lower bound's translation) and (b) the 2013 conjecture is accepted as the repair
mechanism. If the 2013 conjecture fails an internal-consistency test, the upper
portion of the family loses its repair and the window collapses to the subset where
the Markovian form is well-defined without it. **The execution's `conditional_on`
list captures this; no change required, but the audit adds: the window is a
CONDITIONAL FAMILY region, not a single-model region.**

## 4. Attack 4 — the selector conclusion

| Route re-attacked | Result |
|---|---|
| Hidden microphysical scale? | Re-checked: no source derives R₀ from an independent scale (e.g. nucleon radius was a naturalness choice, now excluded; no gravitational derivation exists). **ABSENT confirmed.** |
| Consistency condition? | The 2013 repair exists *because* the unregularized form has an energy defect; R₀ is the repair's knob, not an output. **ABSENT confirmed.** |
| Symmetry/dimensional? | Dimensional analysis constrains ranges; no symmetry principle selects a value. **ABSENT confirmed.** |
| Dynamical selection? | No mechanism in any searched source. **ABSENT confirmed.** |
| Post-hoc classification check | The natural parameter-free choice predates the experiment — correctly NOT classified as post-hoc fit. The surviving window is defined by the bound, not by a fitted value. **CONFIRMED.** |

**Audit conclusion: `FREE_PHENOMENOLOGICAL_PARAMETER` survives.** No hidden selector
was found. The strongest hostile counter-consideration — that some deeper theory
might derive R₀ — remains a speculative future possibility, not current evidence,
and is correctly absent from the execution's positive claims.

## 5. Attack 5 — the prediction firewall

| Question | Finding |
|---|---|
| Distinctive observable? | The radiation channel (which produced the Gran Sasso bound) is the strongest lever. Within the surviving window, **stronger versions of the same experiment probe the upper region**; this is a genuine, falsifiable, quantitative channel. |
| Distinguishable from environmental effects? | The execution correctly notes the environmental-confound firewall is unresolved for much of the window. The discriminator exists **in principle** via the differing d/mass dependence (DP ~ E_G(d)/ħ vs. linear-in-d environmental), but no current experiment has demonstrated the separation across the full window. |
| Does the window constitute a prediction? | **NO — correctly not claimed.** The window is a surviving parameter territory, not a prediction. The radiation channel is a *test vehicle*, not a *verified prediction*. The execution's firewall held. |

## 6. Discrepancy summary

| Finding | Severity | Action |
|---|---|---|
| Lower-bound applicability should carry explicit "standard mass-density-coupled formulation" qualifier | Minor scope note | Record in audit; no execution change required |
| DP-1989 exclusion is by implication through the R₀ family, not direct measurement | Minor scope note | Record in audit |
| Upper-bound arbitrariness correctly flagged | — | Retained |
| Window is family-conditional, not single-model | Clarification | Retained in audit; no execution change |
| Selector conclusion survives | — | Confirmed |
| Prediction firewall held | — | Confirmed |

**No substantive scientific error found.** Two minor scope notes recorded above.

## 7. Audit verdict

**CONFIRMED_WITH_MINOR_SCOPE_NOTES**

The DP-R₀-01 result is defensible as executed:

1. Lower bound: **HARD_EMPIRICAL**, variant-scoped (mass-density-coupled Markovian family).
2. Upper bound: **EFFECTIVENESS_CRITERION**, correctly refused as data.
3. Window: **conditional family region**, splice legitimate under stated conditions.
4. Selector: **FREE_PHENOMENOLOGICAL_PARAMETER** — survives hostile attack.
5. Prediction firewall: held — no inflated prediction claimed.
6. The parameter-free DP version remains **empirically excluded**.

The result stands as: *a conditional DP finite-R₀ model retains a phenomenological
viability region, but no independently demonstrated principle selects R₀ within it,
and the parameter-free version that supplied much of DP's compression advantage is
empirically excluded.*

## 8. Authorization status

- Execution: CONFIRMED, ready to commit (as a separate commit from the freeze package).
- Compression rescoring: **AUTHORIZED AS A SEPARATE PHASE** following this audit, using
  the frozen scoring rules — DP's effective compression has degraded (parameter-free
  version excluded; six-order window carries three burdens). Do not manually lower;
  apply the frozen rules and report the result whatever it is.
- Born audit (E&C-BORN-01): **NOT YET AUTHORIZED** — should be decided after the
  compression rescoring updates the escape-door map.

---

Canonical GRUT-RAI: UNMODIFIED. All artifacts: TestingGRUT laboratory record.
