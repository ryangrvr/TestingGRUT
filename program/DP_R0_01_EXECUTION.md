# DP-R₀-01 EXECUTION REPORT

> Diósi–Penrose R₀ parameter-space reconstruction. Executed under the frozen protocol
> (`DP_R0_01_PROTOCOL.md`). Audit artifact only: no rescoring of the compression matrix,
> no model repair, no viability declaration. Canonical GRUT-RAI untouched.

---

## 1. Verdict

**CONDITIONAL_SURVIVING_REGION**

$$
4\times10^{-10}\,\text{m} \;\lesssim\; R_0 \;\lesssim\; 10^{-4}\,\text{m}
$$

with three explicit conditional burdens:

1. the finite-R₀ construction itself remains a **2013 CONJECTURE** (mass-density resolution);
2. the surviving dynamics is the **Markovian (2007) form**, not the original 1989 dynamics;
3. the upper boundary is an **effectiveness criterion**, not a law of physics — its own source flags it as ambiguous and notes DP "strictly speaking" does not satisfy the classicality requirement at face value.

The selector question was answered **after** parameter-space reconstruction, per protocol order.

## 2. Variant registry (excerpt)

| Variant | R₀ role | Status | Fate |
|---|---|---|---|
| DP-1989 (unregularized) | none (δ-like density) | NOT_APPLICABLE | known energy defect; effectively excluded |
| DP-2007 Markovian | Gaussian smearing width | PHENOMENOLOGICAL, conditionally constrained | open |
| DP-2013 finite-R₀ | mass-density resolution | CONJECTURAL_PHYSICAL_SCALE + free parameter | open, conditional on 2013 conjecture |
| DP natural (parameter-free) | R₀ ~ nuclear size | naturalness choice (post-hoc-flavored, pre-2021 good-faith) | **EXCLUDED** by Gran Sasso |

## 3. Bound registry (types never merged)

| Bound | Type | Parameter | Result | Status |
|---|---|---|---|---|
| Ghirardi et al. 1990 (heating) | CONDITIONAL_EMPIRICAL | R₀ lower | weak | superseded |
| GW-detector noise budgets | CONDITIONAL_EMPIRICAL | R₀ lower | intermediate | superseded |
| Neutron-star heating | CONDITIONAL_EMPIRICAL | R₀ lower | intermediate | superseded |
| **Gran Sasso 2021** ([Donadi et al., Nat. Phys. 17, 74](https://arxiv.org/abs/2111.13490)) | **HARD_EMPIRICAL** | R₀ lower | **R₀ ≳ 4 Å** (~3 orders above previous) | active |
| **Effectiveness criterion 2024** ([arXiv:2406.18494](https://arxiv.org/html/2406.18494v2)) | **EFFECTIVENESS_CRITERION** | R₀ upper | R₀ ≲ 10⁻⁴ m (for the graphene-plate construction used) | active but criterion-dependent |

Key provenance facts:

- The Gran Sasso result is a hard empirical lower bound measured via collapse-induced radiation from charged particles, with dedicated underground background characterization. Its quoted conclusion: the bound "rules out the natural parameter-free version of the Diósi–Penrose model."
- The 2024 upper bound is **not** an experimental result. It is the requirement that collapse make the smallest visibly distinguishable object classical within perception time; the source states plainly: *"Clearly there is a lot of ambiguity in this request, given that the macroscopic domain cannot be defined precisely."* It also reports DP "strictly speaking" fails that requirement and only becomes effective for larger objects.

## 4. Surviving window — boundary status

| Boundary | Status |
|---|---|
| Lower: 4 Å | **HARD_EMPIRICAL_BOUND** (Gran Sasso) |
| Upper: ~10⁻⁴ m (10⁶ Å) | **EFFECTIVENESS_CRITERION** — to some degree arbitrary; depends on chosen object, τ_obs, d; not data |
| Entire window | conditional on 2013 conjectural repair + Markovian form |

$$
\text{surviving window} \neq \text{physically grounded candidate};\qquad \text{falsifiable} \neq \text{derived}.
$$

## 5. Distinguishing observables within the window

| Channel | Distinguishing? |
|---|---|
| Collapse-induced radiation (charged particles) | **YES** with adequate background — strongest existing lever |
| Optomechanical/oscillator heating | CONDITIONAL on environmental noise budget |
| Interferometric visibility loss (mesoscopic superpositions) | CONDITIONAL — must beat Γ_env by design; DP's dependence on superposition distance d (through E_G) differs from linear-in-d environmental decoherence |
| BEC superposition tests | proposed, not yet decisive |

**Environmental-confound firewall:** unresolved for much of the window at current sensitivity. Any future discriminator claim must show Γ_DP(d, R₀) separation from the standard decoherence budget, using the differing d/mass dependence.

## 6. R₀ selector investigation (run only after §4)

| Route | Result |
|---|---|
| Derivation from gravitational field equation | ABSENT — τ ~ ħ/E_G never derived |
| Derivation from QFT/UV structure | ABSENT |
| Consistency requirement (energy conservation) | ABSENT — finite R₀ **is** the conjectural repair; the scale is the repair's free ingredient |
| Symmetry/dimensional argument | ABSENT — dimensional analysis fixes ranges, not values |
| Independent experimental constraint | PARTIAL — hard lower bound (Gran Sasso); upper bound only criterion-based |

**Verdict: FREE_PHENOMENOLOGICAL_PARAMETER** (with CONDITIONALLY_CONSTRAINED status).

Post-hoc-fit check: the natural parameter-free choice (R₀ ~ nuclear size) was made in good faith **before** the 2021 experiment and was excluded by it. The surviving window is therefore defined by experiment, not by an independent selection principle. No mechanism selects R₀ anywhere in the 6-order interval.

## 7. Compression-matrix impact

**Not rescored** (separate authorized phase). Recorded observation only: DP's compression advantage in the comparative matrix derived substantially from the parameter-free version; that version is empirically excluded. DP's effective compression has degraded to a phenomenological window with three attached burdens (conjectural repair + free scale + environmental discrimination).

## 8. Limitations

- Bounds quoted from source abstracts/figures; exact confidence intervals not re-derived from raw data.
- Upper bound criterion-dependent and flagged as ambiguous by its own source.
- Non-Markovian DP variants and relativistic formulations: open, not fully bounded here.
- B4's applicability assumes standard nucleon mass-density composition.

## 9. Sources

- [Donadi et al., Nat. Phys. 17, 74–78 (2021), arXiv:2111.13490](https://arxiv.org/abs/2111.13490) — EXPERIMENTAL_RESULT (primary): R₀ ≳ 4 Å; parameter-free DP excluded.
- [arXiv:2406.18494 (2024), "On the effectiveness of the collapse in the Diósi–Penrose model"](https://arxiv.org/html/2406.18494v2) — ANALYSIS: R₀ ≲ 10⁻⁴ m under the stated (ambiguous) effectiveness criterion.
- [Diósi, J. Phys. Conf. Ser. 442, 012001 (2013)](https://www.nature.com/articles/s41567-020-1008-4) — CONJECTURE: mass-density resolution; R₀ free.
- Bassi et al., Rev. Mod. Phys. 85, 471 (2013) — REVIEW (corroboration).

---

Canonical GRUT-RAI: UNMODIFIED. All artifacts: TestingGRUT laboratory record.
Machine-readable companion: `program/DP_R0_01_EXECUTION.json`.
