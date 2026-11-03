# Track e2e_card: integrator of the end-to-end pass, and skeptic for e2e_sm and e2e_higgs

**Deliverables:**
- [END_TO_END.md](END_TO_END.md): the Lagrangian, the history, the ledger, supplied vs computed, the open items and the
  fine-tuning.
- [CANDIDATE_M2.json](CANDIDATE_M2.json): the full parameter card, written by `code/make_card.py`.

**Ground rules.**
- Every number rests on model choices.
- Nothing here is a validation, an exclusion, a new theory or a GRUT prediction.
- Data comparisons are retrospective.
- Budget was tight, so the earlier tracks' results are reused read-only and only the decisive numbers were recomputed.

## 1. Skeptic checks of e2e_sm

| Item | Verdict | Evidence |
|---|---|---|
| Σm_ν = 0.0588 eV (NO, m₁ = 0, NuFIT 6.0) | **holds** [D] | `code/sm_spot.py` |
| CP-even Casas–Ibarra sum rule Σ η_j m̃_j = Σ m_i | **holds** [D] | Own derivation: R_ij² = η_j\|R_ij\|², and the rows of R R^T = 1 give the rule. Numerically 0.05878 for S1 |
| m̃_j, R R^T = 1, y²_jj, Γ_j = (0.533, 2.49, 0.341, 0.123) H_i | **hold** [C] | Rebuilt from sm_card.json |
| ε_j per unit δ = (2.07, 0.32, 0.96, 0.67)×10⁻² | **holds** [C], to 3×10⁻⁴ | Own CRV code, with e2e_sm's complex CP-odd direction A = where(η_jη_k > 0, i, 1)·A_ODD |
| The card stores the CP-odd direction | **incomplete** | sm_card.json stores only A_ODD. The real direction gives ε smaller by ×0.47–0.72 and O(δ) m̃ shifts. CANDIDATE_M2.json states the full direction |
| **η_B/δ = 1.74×10⁻⁴; δ_needed = 3.5×10⁻⁶; ×540 non-thermal enhancement** | **refuted** [C] | See below |
| Spontaneous asymmetry erased (survival ≤ 10⁻⁴⁵ from t ≤ 30) | holds | Asymmetry made that early meets the full later washout. No θWW̃ in the mass-phase basis [A], consistent with 2006.03148 (verified) |
| UV completion: θ̇/f max 7.4 (σ), 11.1 (χ) | **σ underestimated** [C] | The bench's homogeneous σ̇ peaks at θ̇/f₁ = 11.6 at t = 0.204 (energy-based 7.5 at t = 0.19). Weak coupling then needs λ ≥ 270, and θ̇/m_ρ = 0.92 at m_ρ = 4πf. The negative verdict holds and is stronger |
| ΔN_eff = 0.0009–0.0043 (N_s dark radiation) | **holds on its history**; larger at the M2 points | ΔN_eff = 2.86 × f_s,dec [D]. With e2e_chain's f_s at the s_* floors: 0.004 (M2-X) and 0.017 (M2-Q); at the card's s_* = 1.5×10¹⁶ GeV: 0.008 and 0.021. All lie far below the Planck/ACT errors |
| Contradiction with e2e_chain, which lets s decay via N_s → LH | **settled: dark radiation adopted** | The chain's +0.07% / +0.26% H_req entropy shift is dropped |

### The baryon-asymmetry error (decisive)

**What sm_lepto does.**
- It sets X_NT,j = f_Nj ρ a³ / E_j, read from the e2e_kinetic export, and adds ε_j Γ_j (M/E) X_NT as a source. That
  term has no inverse-decay counterpart.
- The thermal N_j are evolved a second time through X_th → X_eq.

**What the export actually contains.** I compared f_N4 with the Boltzmann equilibrium energy density at T_thermal:
- 1.00–1.02 at t = 100–300;
- 2.4–3.6 at t = 10³–5×10³, where χ → N₄N₄ feeds it;
- 1.23 at 7×10³;
- 1.00 at t ≥ 8×10³.

**Why that is an error.** At late times the "non-thermal" source is the equilibrium population. An equilibrium
population makes no asymmetry, and its tail at z₄ = 6–10 is only weakly washed out (survival 0.05–0.9). That tail made
the ×540. It is history-independent, which is why η_B/δ came out at 1.737×10⁻⁴ to four digits on the kinetic, GM and
P0 histories alike.

**Correction** (`code/lepto_fix.py`). e2e_sm's code is imported unchanged; only the history is patched, so that f_Nj
becomes the excess over equilibrium.

| History | η_B/δ, washout ×1 | δ_needed | ×0.5 | ×2 |
|---|---|---|---|---|
| kinetic HFa (1.3×10¹¹) | 8.9×10⁻⁷ | 6.9×10⁻⁴ | 8.6×10⁻⁶ (7×10⁻⁵) | 3.8×10⁻⁷ (1.6×10⁻³) |
| final GM 1.3×10¹¹ / 8.4×10¹⁰ | 8.5 / 5.7×10⁻⁷ | 7.2×10⁻⁴ / 1.1×10⁻³ | | |
| final P0 1.3×10¹¹ / 8.4×10¹⁰ | 8.6 / 5.7×10⁻⁷ | 7.1×10⁻⁴ / 1.1×10⁻³ | | |
| **M2-X** (ln H interpolation) | 7.5×10⁻⁷ | **8.2×10⁻⁴** | 9.5×10⁻⁵ | 1.9×10⁻³ |
| **M2-Q** | 6.7×10⁻⁷ | **9.2×10⁻⁴** | 1.1×10⁻⁴ | 2.1×10⁻³ |

**Independent minimal model** (`code/lepto_check.py`, own code).
- Source: condensate energy loss → N₄ pairs (Y_N4 = 0.17–0.18).
- Washout: inverse decays of all four N, plus ΔL = 2.
- Result: a non-thermal part of 1.1–1.4×10⁻⁷ δ (survival 4×10⁻⁵, source-weighted T = 3.1–3.3×10¹² GeV, z₄ ≈ 5). Added
  to e2e_sm's thermal-only 3.2×10⁻⁷ δ, that gives about 4.5×10⁻⁷ δ, within ×2 of the corrected value.
- Washout ×0.5 / ×2 changes the non-thermal part by ×60 / ×10⁻³: it is exponentially sensitive.

**Verdict.**
- Baryogenesis is possible. η_B = 6.12×10⁻¹⁰ needs δ ≈ 10⁻³ (range 10⁻⁴–2×10⁻³), not 3.5×10⁻⁶. O(1) phases still
  overproduce, by ×10²–10³.
- The tuning is ~10⁻³ against the χ vacuum phase 1.31 rad [A].
- The result is unflavoured and has no ΔL = 1 scatterings [E]. T ≈ 1–3×10¹² GeV sits at the τ-flavour edge.

## 2. Skeptic checks of e2e_higgs

| Item | Verdict | Evidence |
|---|---|---|
| ξ_min during inflation: 0.046–0.055 | **consistent** [E] | Own constant-H Hawking–Moss criterion (8π²V_b/3H⁴ > 3N₀) on e2e_higgs's cached λ_eff tables (read-only) gives ξ(μ) = 0.014–0.022, i.e. ξ(M_t) = 0.039 / 0.044 at 172.57 / 172.86. That is within 20% of the track's tree-level 0.049 / 0.053. ξ_H = 0.1 keeps a margin of about ×2 |
| N-gas medium mass ≥ 10³ H_i² by t = 0.2 | **holds** [E] | Own Σ y_j² ρ_Nj/E_j² on the kinetic history: 1.6×10⁴ at t = 0.2, 8.8×10⁵ at 0.3, 2.9×10⁵ at 3 (H_i²) |
| M_t,crit at the final H_i | **resolved** (e2e_higgs had to extrapolate to 3.7–6.5×10¹¹) | The final roots, 1.123 / 0.989×10¹¹ GeV, lie inside the computed range. Quadratic in ln H_i: 172.94–172.96 (M2-X), 172.98–173.00 (M2-Q); 172.63 / 172.67 at α_s = 0.1171 |
| −0.30 GeV per e-fold of H_i | holds | 5×10¹⁰ → 2.4×10¹¹: −0.48 GeV over 1.57 e-folds |

## 3. Contradictions settled as integrator

1. **K: lattice continuation vs kinetic vs chain.** The chain's energy-conserving lattice → kinetic histories are
   adopted (GM 15.3, P0 8.8, GM with spread 12.6). T_RH and t_SM,50% are identical for P0 and GM (3.15×10¹² GeV and
   4.0×10³ at 1.3×10¹¹); the K difference is early radiation [C, read from the chain's fluid exports].
2. **M2-Q "only on P0".** The skeptic's pocket on GM with the fragment spread is recorded as an alternative, and the
   wider phase window is adopted.
3. **μ_eff components.** The m2_points.json values are used; the NOTES/digest values were wrong.
4. **N_s.** Dark radiation (see §1).
5. **e2e_higgs H_i range.** Resolved by interpolation (§2).

## 4. Files

**code/**
- `lepto_fix.py`: corrected rerun of e2e_sm's sm_lepto, imported unchanged.
- `lepto_check.py`: own minimal leptogenesis model.
- `sm_spot.py`: neutrino sector, CP asymmetries and the σ θ̇/f check.
- `higgs_spot.py`: HM ξ_min, M_t,crit interpolation and N-gas mass.
- `make_card.py`: writes the card.

**results/**
- `lepto_fix.json`, `lepto_own_model.json`, `lepto_summary.json`;
- `sm_spot.json`, `higgs_spot.json`, `higgs_xi_running.json`.

**Read-only inputs, besides the tracks' results:** e2e_chain's scratch `kin/fluid_ho*.json` (the final histories) and
e2e_higgs's scratch `smcache`. Temporary patched histories are in this track's scratch.

## 5. Limits

- The leptogenesis correction reuses e2e_sm's simplified Boltzmann system and inherits its omissions. The "excess over
  equilibrium" split is a model choice [A]. The full fix is a momentum-resolved N₄ source in e2e_kinetic, exported as a
  separate non-thermal column.
- The Higgs checks are estimates. No mode or danger runs were made.
- Not re-checked: e2e_sm's one-loop light-mass correction, the S2 scenario, the N_s soft-mass variant, and the lattice
  and kinetic internals (checked only through their exports).
