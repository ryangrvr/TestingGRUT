# M2 end to end: the Lagrangian, the history, the observables ledger

*Integrator: e2e_card, end-to-end pass, 3 November 2026.*

M2 is M1r's hidden sector together with the HFa multi-field inflaton sector. HFb is an alternative that changes only
n_s. There are two benchmarks:
- **M2-X:** free X at e = 10⁻⁶, with e_crit ≈ 2×10⁻⁶.
- **M2-Q:** Q-matter at e = 10⁻⁶.

**Ground rules.**
- Every number rests on model choices: the potentials, couplings, closures, splice rules, relic rule and seed rule.
- Closure fixes none of the supplied numbers. Nothing here is a GRUT prediction, and nothing validates or excludes
  anything.
- All data comparisons are retrospective.
- Every parameter, with its status, source and window, is in [CANDIDATE_M2.json](CANDIDATE_M2.json).
- Claim labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## 0. Verdict in one paragraph

**Every link from CMB horizon exit to today now has a computed number, given a set of supplied closures.** No link
breaks the model:
- no domain walls form;
- reheating completes;
- baryogenesis is possible;
- the Higgs survives.

**Four links are weaker than the tracks reported:**
1. **Pauli blocking.** Which blocking closure is physical is still open. That closure sets K (8.8 or 15.3) and decides
   how robust M2-Q is.
2. **The axial couplings' UV completion.** It is in control only at strong coupling, and only at the edge of control.
3. **The baryon asymmetry was overcounted by e2e_sm** [C, this track]. The corrected η_B needs a CP phase of about
   10⁻³, not 3.5×10⁻⁶. The result is exponentially sensitive to the washout normalisation.
4. **The relic branch is a tuned choice.** It is set by the holder strength c0 to a few percent, by the phase to a few
   degrees, and by the spectator's initial value s_*, which no equilibrium supplies.

So M2 is complete as a model. It is not a consequence of GRUT, and it is not natural.

## 1. Lagrangian (signature −+++, reduced Planck mass M_P)

```
L/√−g = M_P² R/2
   − ½(∂σ)² − ½(∂χ)² − ½(∂s)² − V(σ,χ) − V_s(s;χ,R) − λ_s s⁴
   + L_SM − ξ_H R |H|²
   + Σ_{j=1..4} [ ½ N̄_j i∂̸ N_j − (M_j/2) N̄_j e^{iγ5 φ_j/f_j} N_j ] − (ỹ_{αj} L̄_α H̃ N_j + h.c.)      φ_1 = σ, φ_{2,3,4} = χ
   + ½ N̄_s i∂̸ N_s − (y_s/2) s N̄_s N_s
   + L_hidden(X, φ, Y, Σ, A_μ) − (g_h e^{iα} s Σ* Y² + h.c.)

V(σ,χ) = (λ/4)(χ² − v²)² + ½ g² σ² χ² − c χ + c₁(σ − σ_c) − ½ m₂ (σ − σ_c)² + V_1-loop(μ_R² = κ g²(2σ_c)²) + V_const
V_s    = ½ g_s² χ² s²                      (HFa; HFb: ½ (m_s0² + ξ_s R) s², m_s0 = 3H_i, ξ_s = −3/4)
L_hidden = −|D X|² − ½(∂φ)² − |D Y|² − |D Σ|² − ¼ F² − λ_X|X|⁴ − g² φ²|X|² − λ_φ/4 (φ² − v_φ²)² + ε v_φ (φ − v_φ)³
         − M_Y²|Y|² − λ_Y|Y|⁴ − λ_Σ(|Σ|² − h²/2)² − (μ Σ* Y² + κ Y* X³ + h.c.) − ξ_X R |X|²
```

Charges are X 1, Y 3, Σ 6. Choices made in this Lagrangian:
- **Portals.** No tree-level portals between {σ, χ, s} and the Higgs, and none between {σ, χ} and the hidden
  sector. The induced holder portal vanishes at M_Y (c_L = 0).
- **N_s.** It is massless, with N_s → iγ5 N_s exact, so it has no seesaw Yukawa.
- **UV completion [A].** Φ_a = (ρ_a/√2) e^{iθ_a}, with σ = f₁θ₁ and χ = f₂θ₂. The hybrid potential is a monodromy on
  these angles (Kaloper–Sorbo 0811.1989; McAllister–Silverstein–Westphal 0808.0706). The radial modes are taken at
  m_ρ = 4πf.

Values are in CANDIDATE_M2.json, under `common`.

## 2. The history, stage by stage

Units: t is measured from the slice (σ = σ_c, a = 1) in units of H_i⁻¹.
- M2-X: H_i = 1.123×10¹¹ GeV, so H_i⁻¹ = 5.9×10⁻³⁶ s. Closure GM.
- M2-Q: H_i = 9.89×10¹⁰ GeV, so H_i⁻¹ = 6.7×10⁻³⁶ s. Closure P0.

Sources are given in brackets.

| # | Stage | M2-X | M2-Q | Label, source |
|---|---|---|---|---|
| 0 | **CMB exit** | N_* = 51.55 before the slice; σ_* = 0.021 M_P on the attractor; s = s_* (frozen) | N_* = 51.76 | [C] e2e_chain + skeptic, mf_hybrid |
| 1 | **Slow roll.** σ falls to σ_c. χ sits at χ_min = c/M_χ² > 0. s holds the hidden phase: g_h s_* = 19.07 H_i (M2-X) / 17.84 H_i (M2-Q), i.e. μ_eff = −3.94 + 15.50i / −3.22 + 14.50i | V0^{1/4} = 6.88×10¹⁴ GeV | 6.46×10¹⁴ GeV | [C/D] e2e_chain (corrected μ_eff, skeptic) |
| 2 | **Waterfall** | ε_H = 1 at t = 0.19 (H = 0.997 H_i). χ → +v everywhere (bias). The spectator mass turns on in ≈ 0.2 | same | [C] mf_hybrid, e2e_lattice |
| 3 | **Preheating and fragmentation** | N production by the axial phases up to 36% of the energy at t = 0.3–0.45. f_coh < 0.1 at t = 0.47 (GM). Gradient energy peaks at 0.22. The −v fraction peaks at 3% and is transient; no cluster percolates; **no walls** | the same at t = 0.49 (P0); −v peak 1.8% | [C] e2e_lattice (64³/96³), e2e_kinetic |
| 4 | **Hidden release** (t ≲ 3) | ε = 0.74, E/E_c = 1.49, N_YL = 35, N_YH = 23 | ε = 0.61, E/E_c = 1.91, N_YL = 6.4, N_YH = 9.2 | [C] e2e_chain on the final histories |
| 5 | **Hand-over at t_sw = 3** (lattice → kinetic, energy-conserving) | H continuous to 5×10⁻⁴ | same | [C] e2e_chain (fixed a 19%/75% N-energy loss) |
| 6 | **SM bath** | Thermal by t ≈ 7.4. T_max ≈ 4.3×10¹³ GeV at t ≈ 2 | thermal by t ≈ 6.5; T_max ≈ 4.7×10¹³ GeV | [C] e2e_kinetic network, interpolated in ln H [E] |
| 7 | **χ-remnant era.** χ → N₄N₄, Pauli-blocked by its own products; the remnant dominates | SM reaches 50% at t ≈ 4.7×10³; T_RH = 2.7×10¹² GeV; K = 15.7 | t ≈ 5.3×10³; T_RH = 2.4×10¹² GeV; K = 9.05 | [C] e2e_chain final histories |
| 8 | **Leptogenesis** at T ≈ 1–3×10¹² GeV, z₄ ≈ 5–10. Strong washout (K_j = 46–85). Sphalerons convert at T ≲ 10¹² GeV | η_B/δ = 7.5×10⁻⁷ ⇒ δ = 8.2×10⁻⁴ | 6.7×10⁻⁷ ⇒ δ = 9.2×10⁻⁴ | [C on Q/E] e2e_card (corrected e2e_sm, §5) |
| 9 | **Heavy remnants and the spectator.** Y_L and Y_H decay deep in radiation (Γ_L = 10⁻⁷ H_i); their heat is in the ledger. s decays at a ≈ 4.8×10³ into N_s (dark radiation) | s energy fraction 0.27% (at s_* = 1.5×10¹⁶ GeV) | 0.73% | [C/E] e2e_chain, e2e_sm |
| 10 | **Hidden relaxation and the branch** at e = 10⁻⁶ | S/Q = 4.03 > 3.47: **free X**; e_crit = 1.72×10⁻⁶ (bulk) | S/Q = 2.43: **Q-matter**; e_crit = 0.43×10⁻⁶ | [C] e2e_chain + skeptic |
| 11 | **Relic** | Asymmetric free X: M = 1.0003 m_X, X̄/X = 1.6×10⁻⁴ | Q-matter: ⟨Q⟩ = 7.5×10¹⁹ by ripening; M = 0.802 m_X; vapour f_v = 0.094–0.138 | [C] e2e_relic relic() (⟨Q⟩ transferred [Q]) |
| 12 | **BBN** | The only new light state is N_s (dark radiation); every other new state has decayed by about 10⁻²⁹ s | same | [D] e2e_sm |
| 13 | **Today** | Ω_DM fitted through H_i; n_X = 3.2 cm⁻³ locally | 9.6 μg (0.44 M_Pl) objects of radius 0.73 nm; n = 4.9×10⁻¹⁴ m⁻³ | [C] e2e_relic, e2e_chain |

## 3. Observables ledger (retrospective)

| Observable | M2-X | M2-Q | Bound or measurement (verified arXiv) | Label |
|---|---|---|---|---|
| n_s | 0.9649 (HFa target; HFb 0.9743) | same | Planck 0.9649 ± 0.0042 (1807.06211); P-ACT 0.974 ± 0.003 (2503.14452), so HFa sits at −3.0σ and HFb at −2.2σ vs Planck | **fitted** (m₂), not an output |
| r | 2.05×10⁻⁷ | 1.59×10⁻⁷ | < 0.036 at 95% (BICEP/Keck, 2110.00483) | [C] MS modes; ±4–5% (H_i rule) |
| α_s (running) | +9.5×10⁻⁵ | +7.4×10⁻⁵ | Planck −0.0045 ± 0.0067 [R, 1807.06211] | [C] |
| CDM isocurvature β_iso | 0.022 at s_* = 1.5×10¹⁶ GeV (floor 1.13×10¹⁶) | 0.031 (floor 1.36×10¹⁶) | < 0.038 (Planck 2018 X, 1807.06211) | [C/D]; D = 5.8 / 7.9, ±10% |
| f_NL (local) | ≈ (5/12)(1 − n_s) ≈ 0.015 from the single clock; spectator curvature contribution ~10⁻⁹ of the power | same | −0.9 ± 5.1 (Planck 2018 IX, 1905.05697) | **[E] only, not computed** |
| Ω_c h² | 0.120 | 0.120 | 0.120 ± 0.001 (1807.06209) | **fitted** (H_i, D target 4.4×10⁻¹⁰ GeV) |
| Relic type | asymmetric free X (92.4 MeV) | Q-matter | — | [C] |
| Object mass, size | — | 9.6 μg = 0.44 M_Pl; R = 0.73 nm | — | [C] e2e_relic, skeptic |
| Flux | n = 3.2 cm⁻³, 1.1×10⁸ cm⁻² s⁻¹ | 0.55 m⁻² yr⁻¹ (n⟨v⟩A) or 0.27 (isotropic flat plane); impulse 3.6×10⁻²⁷ kg m/s on 1 mg at 1 mm | Proposed sensors 1903.00492 and Windchime 2203.07242: purely gravitational Planck-mass detection needs about 80 dB noise reduction, and 0.44 M_Pl lies below that target | [C]; sensitivity not assessed |
| Self-interaction σ/m | 6.0×10⁻⁸ (XX), 7.1×10⁻⁷ (XX̄) cm²/g | 7×10⁻⁹ cm²/g (contact) | < 0.7 cm²/g (Bullet Cluster, 0704.0261) | [D/C] |
| P(k) | CDM (free-streaming length 0.3–0.5 pc) | CDM-like (ζ ≈ 10⁻⁷ vapour) | Lyman-α, satellites | [Q] qm_transfer fit |
| Σm_ν | 0.0588 eV (NO, m₁ = 0) | same | < 0.0642 eV (DESI DR2 + CMB, 2503.14744); < 0.12 eV (1807.06209) | [D]; splittings fitted (NuFIT 6.0, 2410.05380) |
| η_B | 6.12×10⁻¹⁰ at δ = 8.2×10⁻⁴ (window 9.5×10⁻⁵–1.9×10⁻³) | at δ = 9.2×10⁻⁴ (1.1×10⁻⁴–2.1×10⁻³) | 6.12×10⁻¹⁰ (Planck ω_b h² = 0.0224, 1807.06209) | **fitted** (δ); [C on Q/E] |
| ΔN_eff | 0.008 (N_s dark radiation; 0.004 at the s_* floor) | 0.021 (0.017 at the floor) | Planck N_eff = 2.99 ± 0.17 (1807.06209); ACT DR6 2.86 ± 0.13 (2503.14454); BBN N_ν = 2.86 ± 0.15 (1912.01132) | [E] 2.86 × the s energy fraction |
| BBN | no new state except N_s; N lifetimes ~10⁻³⁵ s, s lifetime 4×10⁻³⁰ s | same | — | [D] |
| Gravitational waves | Ω_GW h² ~ 10⁻¹²–10⁻¹¹ at ~0.1–3 GHz | same (K⁻² scaling) | far below the ΔN_eff ceiling (~1.7×10⁻⁶); outside every detector band | [C/E] e2e_lattice; amplitude unconverged (factor few–10) |
| Walls | none (−v ≤ 3%, transient, never percolates) | none | — | [C] lattice |
| PBHs | none (waterfall ≈ 0.2 H_i⁻¹; ζ_wf ≲ 6×10⁻³ on sub-Hubble scales) | none | — | [C/E] mf_hybrid |
| Strings | U(1)_a never thermally restored (T_max/f₁ ≈ 0.3); hidden U(1) Higgsed during inflation (h = 40 H_i) | same | — | [D]/[A]; non-thermal restoration and hidden strings not computed |
| Higgs stability | Gravitational-only M_t,crit = 172.96 GeV (172.63 at α_s = 0.1171); with N-gas and SM media, safe to ≥ 174.6 GeV; ξ_H lower edge 0.046–0.055 | 172.99 (172.67) GeV; same | M_t = 172.57 ± 0.29 GeV (PDG 2024 [R]) | [C on E] e2e_higgs, interpolated in ln H_i [E] |

## 4. What is supplied and what is computed

| Supplied (model choices) | Fitted (to a named datum) | Derived or computed |
|---|---|---|
| Architecture (HFa/HFb, axial singlets, spectator holder, hidden M1r) | c₁, m₂ → A_s and n_s (Planck or ACT) | N_*, r, α_s (MS modes) |
| λ, v, g, c, κ (μ_R), and the V0 rescaling rule | H_i → Ω_DM (D = 4.4×10⁻¹⁰ GeV) | K, T_RH, f_rad(t), Γ_th(t) (lattice and kinetic) |
| M_j, f_j, m̃_j (S1), the CP parities, the CP-odd direction A | light-ν splittings and PMNS → NuFIT 6.0 | Absence of walls; fragmentation; GWs |
| g_s (HFa) or m_s0, ξ_s (HFb); y_s; λ_s (≥ 7.48 g_h⁴) | δ (lepton CP phase) → η_B | Release ε, E/E_c, N_YL, N_YH; the S/Q ledger |
| **c0 = 1.55 / 1.45, phase δ_h = 0°**, α + sign(s_*) (the DM sign) | — | S/Q(e), e_crit, relic M, ⟨Q⟩ ripening, object observables |
| **s_* = 1.5×10¹⁶ GeV** (initial condition; not an equilibrium) | — | Isocurvature D and the s_* floor; ΔN_eff |
| **e = 10⁻⁶** | — | Σm_ν; ε_j; washout K_j |
| Hidden: g, λ_X, λ_φ, v_φ (m_X), ε, M_Y, h, μ, κ, λ_Y, λ_Σ, ξ_X; c_L = 0; sequestering | — | M_t,crit, ξ_H window, medium deadlines |
| ξ_H = 0.1; zero portals; N_s massless | — | — |
| Closures: Pauli GM (M2-X) / P0 (M2-Q); hand-over at t_sw = 3; relic rule; seed rule; leptogenesis source | — | — |
| ρ_vac (measured, underived) | — | — |

## 5. Settled contradictions and corrections (this pass)

1. **K.** The lattice's own fluid continuation (K = 3.6–6.1) lacks the χ → N₄N₄ self-blocking. The kinetic K = 16.35 has
   a homogeneous parent. Both are superseded by e2e_chain's energy-conserving lattice → kinetic histories:
   - GM: K = 15.3;
   - GM with a fragment spread σ = 0.3: K = 12.6;
   - P0: K = 8.8.

   T_RH (3.15×10¹² GeV at 1.3×10¹¹) and t_SM,50% are the same for P0 and GM. K differs because P0 puts more energy into
   radiation early [C, read from the final histories].
2. **Does M2-Q exist only on P0?** Not quite (e2e_chain skeptic).
   - It is robust on P0.
   - On GM with the fragment spread there is a narrow, bulk-only pocket: c0 = 1.485–1.52, δ_h = +2…+3.5°.
   - On plain GM it is absent.
   - The phase window on P0 is at least −4…+4°. c0 = 1.40 on P0 is bistable.
3. **The μ_eff components** in e2e_chain's NOTES and digest were wrong. The values in m2_points.json are correct
   (−3.94 + 15.50i and −3.22 + 14.50i).
4. **N_s.**
   - e2e_chain let s decay to SM entropy (N_s → LH); e2e_sm fixes N_s as massless dark radiation. **Adopted: dark
     radiation.**
   - The chain's +0.07% / +0.26% H_req entropy shift is dropped (it is below the ±2% H_i systematic).
   - ΔN_eff = 2.86 × the s energy fraction at decay [D]: 0.004 / 0.017 at the s_* floors, 0.008 / 0.021 at the card's
     s_*.
5. **Baryogenesis: e2e_sm overcounted η_B** [C, e2e_card/code/lepto_fix.py].
   - sm_lepto reads the "non-thermal" N population from the kinetic f_Nj columns. But f_N4 equals its
     thermal-equilibrium value to 1–2% at t ≤ 300 and t ≥ 8000. It is ×1.2–3.6 above equilibrium only at
     t = 10³–7×10³, where χ → N₄N₄ feeds it.
   - Adding εΓX_NT for the equilibrium part, with no inverse-decay counterpart, creates asymmetry from an equilibrium
     population. The late, weakly washed-out tail (survival 0.05–0.9 from t ≥ 10⁴) gave the ×540 "non-thermal"
     enhancement. This explains why η_B/δ came out identical, to four digits, on every history.
   - With only the excess over equilibrium as the source, e2e_sm's own code (unchanged, imported) gives η_B/δ:

     | washout | kinetic history | M2-X | M2-Q |
     |---|---|---|---|
     | ×1 | 8.9×10⁻⁷ | 7.5×10⁻⁷ | 6.7×10⁻⁷ |
     | ×0.5 | 8.6×10⁻⁶ | — | — |
     | ×2 | 3.8×10⁻⁷ | — | — |

   - So δ ≈ 8–9×10⁻⁴, within 10⁻⁴–2×10⁻³.
   - An independent minimal model (code/lepto_check.py: condensate loss → N₄ pairs, inverse decay plus ΔL = 2 washout)
     gives a non-thermal part of 1.1–1.4×10⁻⁷ δ (survival 4×10⁻⁵). Added to the thermal 3.2×10⁻⁷ δ, that is
     4.5×10⁻⁷ δ, within ×2 of the corrected value.
   - η_B is still overproduced at O(0.1–1) phases (×100–1000). The χ vacuum phase v/f = 1.31 rad must be cancelled to
     about 10⁻³, not 10⁻⁶ [A].
   - The result depends exponentially on the washout normalisation: the non-thermal part changes ×60 for washout ×0.5
     and ×10⁻³ for ×2.
6. **The weakly coupled UV completion is worse than e2e_sm found** [C, code/sm_spot.py].
   - The bench's homogeneous σ reaches θ̇/f₁ = 11.6 at t = 0.204, against e2e_sm's 7.4, which was estimated from the
     energy fraction at t ≥ 0.19.
   - Freezing ρ then needs λ ≥ 270, and θ̇/m_ρ = 0.92 at m_ρ = 4πf.
   - The conclusion is unchanged, only stronger: only the strongly coupled edge, or a raised f (not rerun), holds.
7. **e2e_sm numbers reproduced by independent code** [C]: Σm_ν, the CP-even sum rule, m̃_j, R R^T = 1, Γ_j and ε_j
   (to 3×10⁻⁴). ε_j needs e2e_sm's complex CP-odd direction A = where(η_jη_k > 0, i, 1)·A_ODD. sm_card.json stores only
   A_ODD; the real-A direction gives ε smaller by ×0.5–0.7. The card now states the direction.
8. **e2e_higgs at the final scales.** e2e_higgs extrapolated to the chain's in-progress H_i (3.7–6.5×10¹¹ GeV). The
   final roots are 1.12×10¹¹ / 0.99×10¹¹ GeV, inside its computed range.
   - The interpolated gravitational-only M_t,crit is 172.94–172.96 (M2-X) and 172.98–173.00 (M2-Q), above M_t + 1σ.
     At α_s = 0.1171 it is 172.63 / 172.67, so the central M_t is safe by 0.06–0.10 GeV.
   - Own Hawking–Moss estimate: ξ_min(μ) = 0.014–0.022 at constant H, which is ξ(M_t) = 0.039 / 0.044 at
     M_t = 172.57 / 172.86. That is within 20% of e2e_higgs's 0.049 / 0.053.
   - Own N-gas medium mass: m_H² ≈ 1.6×10⁴ H_i² at t = 0.2 and 9×10⁵ at t = 0.3, so e2e_higgs's "≥ 10³ H_i² by
     t = 0.2" is confirmed.

## 6. What remains uncomputed, or held only by a closure (honest list)

1. **Pauli-blocking closure (GM vs P0).** It sets K (15.3 vs 8.8) and whether M2-Q is robust. No first-principles
   criterion was computed.
2. **Axial UV completion.**
   - Weak coupling fails the centrifugal test; strong coupling sits at the edge (θ̇/m_ρ = 0.6–0.9).
   - The f ×3–4.4 fix would slow reheating and was not run.
   - Parametric radial excitation and non-thermal U(1)_a restoration were not computed.
   - The θWW̃ cancellation is argued [A]. It is consistent with the coupling-redundancy bookkeeping of Domcke–Ema–
     Mukaida–Yamada (2006.03148); even if a spontaneous asymmetry were made, washout erases it (≤ 10⁻⁴⁵ from
     t ≤ 30).
3. **η_B** is unflavoured, has no ΔL = 1 scatterings or thermal masses, and is exponentially washout-sensitive. The
   kinetic history used m̃₂ = 0.05 eV, against S1's 0.0912 eV (K shifts by about −2%).
4. **Relic microphysics.**
   - ⟨Q⟩ and f_v are transferred from qm_micro and its skeptic, not recomputed on the hybrid history.
   - A seed or nucleation path is assumed for 2.93 < S/Q < 3.47 (the GM pocket and M2-X at 3×10⁻⁶).
   - A mixed-phase relic function at the 3.47 threshold is missing, so the branch map is discontinuous.
5. **No history was built at the roots.** H_i is interpolated in ln H between 8.4×10¹⁰ and 1.3×10¹¹ (±2% rule
   systematic). D is computed at a single node (±10%).
6. **GW amplitude** is not converged in lattice spacing (factor few–10). f_NL was not computed (estimate only). Hidden
   Z6 strings and non-thermal defects were not re-examined.
7. **Higgs.** Thermal tunnelling was not computed. The upper ξ edge was not recomputed. The media are hard-gas
   estimates, and the SM sign is gauge-sensitive.
8. **No link failed outright.** If the axial completion is rejected, M2 is an effective theory of axial couplings
   below an unspecified strongly coupled UV scale. If the GM closure is right, M2-Q exists only as a bulk-only pocket
   that needs a seed path. If neither closure is right, the self-consistent H_i and c0 move. Γ_th = 10 moved them to
   c0 ≈ 1.25–1.30 at H_i ≈ 5×10¹⁰ GeV in the earlier bracket.

## 7. Fine-tuning

| Quantity | Value | Window | Comment |
|---|---|---|---|
| c0 (M2-X) | 1.55 | 1.45–1.60 for self-consistent free X; e_crit ≈ 2×10⁻⁶ only near 1.55 | e_crit moves to 9.6 / 2.3×10⁻⁶ for a ±3% shift in the lattice-stage m_s² |
| c0 (M2-Q) | 1.45 | 1.40 (bistable) – 1.50 on P0; GM pocket 1.485–1.52 | ±3% |
| Holder phase | 0° | M2-X −2…+6°; M2-Q ≥ −4…+4° | N_YL varies 0.04 → 38 within a few degrees near arg μ ≈ 106° |
| s_* | 1.5×10¹⁶ GeV | ≥ 1.13 / 1.36×10¹⁶ (isocurvature) | Not an equilibrium value: V_min(\|u(s)\|) is tachyonic (m²_eff = −6.3×10⁻⁴ H²); P_SY(window) ≲ exp(−6×10⁷); s_* is set before the observable e-folds |
| λ_s | 10⁻¹⁴ | ≥ 7.48 g_h⁴ (3.1 / 1.4×10⁻¹⁵) | Needed for boundedness; the one-loop value is too small |
| ξ_s (HFb) / g_s (HFa) | −3/4 / g_s v = 3H_i | not scanned | Sets m_s = 3 H_i |
| e | 10⁻⁶ | M2-X: < 1.7×10⁻⁶ (×1.7); M2-Q: > 0.43×10⁻⁶ (×2.3) | Supplied |
| Lepton CP phase δ | 8–9×10⁻⁴ | 10⁻⁴–2×10⁻³ | Must cancel an O(1) vacuum phase (1.31 rad) to about 10⁻³ |
| m̃_j | S1 | ×4.1 cancellation in m_ν | — |
| n_s | via m₂ | — | 58% of η is radiative, so n_s depends on the renormalisation choice κ |
| c_L | 0 at M_Y | ≲ 1/16π² | S/Q +0.05–0.1 at 1/16π² |
| f_j | 1.4–1.5 f_edge | — | At the axial EFT edge |

## 8. Literature verified (WebSearch, this or a previous e2e session)

- **CMB.**
  - Planck 2018 VI: 1807.06209 (verified here).
  - Planck 2018 X: 1807.06211 (e2e_chain skeptic).
  - Planck 2018 IX: 1905.05697 (verified here).
  - BICEP/Keck: 2110.00483 (verified here).
  - ACT DR6: 2503.14452 (skeptic) and 2503.14454 (e2e_sm).
- **Neutrinos.** NuFIT 6.0: 2410.05380. DESI DR2: 2503.14744 (e2e_sm).
- **BBN.** Fields et al.: 1912.01132 (verified here).
- **Dark-matter self-interaction.** Bullet Cluster: 0704.0261 (verified here).
- **Spontaneous baryogenesis.** Domcke et al.: 2006.03148 (verified here).
- **Detectors.** 1903.00492 and 2203.07242 (skeptic).
- **Thermalisation.** Kurkela–Lu 1405.6318; Kurkela–Moore 1107.5050; Harigaya–Mukaida 1312.3097; Besak–Bödeker 1202.1288
  (e2e_kinetic).
- **Preheating and gravitational waves.** 0707.0875, 1110.0337, hep-ph/0012142 (e2e_lattice).
- **Stochastic spectator.** Starobinsky–Yokoyama: astro-ph/9407016.
- **Recalled only [R].** PDG 2024 M_t; the Maldacena consistency relation for f_NL; Covi–Roulet–Vissani
  hep-ph/9605319 and Buchmüller–Di Bari–Plümacher hep-ph/0401240 (cited by e2e_sm).
