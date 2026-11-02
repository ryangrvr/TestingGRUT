# SKEPTIC of track inf_chain (key inf_chain_skeptic)

3 November 2026 · Claude, adversarial pass on `tracks/inf_chain` (BRIEF_INFLATON.md; GitHub issue #13 context).
- This is same-session AI work on one supplied model, not outside replication.
- Nothing here is a validation, an exclusion of the real world, or evidence for any sector.
- **Every potential, slice, holder calibration point, width, portal condition, ledger and survival criterion is a model
  choice.** None follows from closure.
- Data comparisons (Planck 2018, BICEP/Keck 2021, ACT DR6) are retrospective.

Labels: [D] derived, [C] computed with my own code (`skeptic/code/`, output in `skeptic/results/`), [E] estimate,
[A] argued, [Q] quoted from a reviewed track, [R] recalled or literature.

**Method policy.** Wherever possible I recomputed with a different code path from the track's:
- my own background, Mukhanov–Sasaki (MS) modes, decay fluid and e-fold bookkeeping;
- my own Python longitudinal-vector integrator;
- my own deposits, yields and an analytic limit for the S/Q closed form.

Two exceptions:
- **The track's engine was re-run only to test one convention.** This used a patched copy built in my scratch area;
  the default build reproduces the track's numbers bit for bit.
- **Two read-only imports.** The Higgs check imports the shared SM running (`higgs_reheat/code/smpot.py`). The
  convention test imports the track's driver (`ch_common.py`), with its scratch path redirected to mine.

Git was read-only and nothing was committed. Every run was in the foreground, with one process. Scratch was
`$SP/wt2/inf_chain_skeptic/`.

**Process note.** A mistyped shell `cp`, with an unset variable, left a 2.7 KB copy of my own
`sk_gamma_switch.json` at `/sk_gamma_switch.json`, at the filesystem root. Removing it was blocked by the harness's
safety check, so it is still there. It is harmless, and the owner can delete it.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

**The track's central conclusions hold under independent recomputation.** On every benchmark the declared e = 10⁻⁶ is
deep in free X (S/Q = 155–987), with no e_crit. The relic is mostly X X̄, and HT6t and T1 miss their design H_i.

Where I could recompute with my own code, the numbers agree:
- TB CMB and reheating: n_s, r, N_*, H_i, K, x₀ and T_RH;
- K for HT6t and T1;
- the TB1 matter-era vector count;
- every ledger sum, S/Q, yield and H_i,req;
- the large-W limit of the S/Q closed form.

Three statements are wrong, and five issues are new.

| Result | Verdict | Recomputed / reason | Label |
|---|---|---|---|
| C1, C2 HT6t / T1 release | **stand** | My own field-and-fluid ODE gives K = 12.340784 / 15.200283 (track 12.340784 / 15.200284; rel. ≤ 3×10⁻⁸). My engine re-run reproduces the dense late means bit for bit. The NOTES' convention remark is corrected in §2 | C |
| C3 TB release (inputs) | **stands** | x₀ = 21.85060 / 17.06279 (rel. ≤ 1.1×10⁻⁷); m_I = 2.207749 / 1.541287; ε_H(S1) = 0.13470 / 0.22714; K Γ^{1/3}/x₀^{2/3} = 0.501281. N_X and N_Y were not rechecked | C |
| C4 Send sign flip | **stands** | It also survives the width-switch convention: N_X = +3.84×10⁶ → +3.32×10⁶ | C |
| C5 vectors | **stands** (TB checked) | The TB1 static count N_A(k ≤ 1) from my Python Stueckelberg integrator agrees to 0.2% at e = 10⁻⁶ and 10⁻⁵ (§1.3). The deposit arithmetic is reproduced. f_eff was not rechecked | C |
| C6 S/Q, no e_crit | **stands**, with a TB margin | Ledger sums and S/Q reproduced (rel. ≤ 2.4×10⁻⁷; e → ∞ ≤ 1.5×10⁻⁴). TB1's "no e_crit at any e" rests on R_Y = 2.08 against a threshold of 0.85, i.e. a factor 2.4 (§3, N6) | C |
| C7 closed form beyond its fit range | **stands; independent support** | Its W → ∞ limit is the ideal massless X + φ gas at fixed charge. That limit gives 0.190463 W^{3/4} against the fit's 0.190255 (0.11%). closed form/ideal = 0.9990–1.0007 for W = 10³–10⁶. N_c/E^{3/4} = 6.9775676 exactly | D |
| C8 p(S/Q) | not rechecked; consistent | The entropy per final particle, S/Q/(1 + 2p), is 8.5–9.9 in the map and 8.96–8.99 on the benchmarks. That is a frozen Maxwell–Boltzmann gas with (m − μ)/T_f ≈ 6–7 | D |
| C9 relic and H_i,req | **stands** | My yield from N_X and K, and my H_i,req, agree to ≤ 6×10⁻⁸ | C |
| C10 iteration | **numbers stand; "diverging" refuted** | The log mismatch ln(H_req/H_pot) goes −1.565 → −0.989, a contraction of 0.63. The secant fixed point is ≈ 2.0×10⁹ GeV (r ≈ 7×10⁻¹¹). T1's is 1.18×10¹⁰ (§4) | C |
| C11 TB width calibration and CMB | **stands** | Independent MS and N_* (§1.1): ΔN_* = +0.07 / +0.08, H_i −0.13% / −0.16%, n_s +5×10⁻⁵ / +6×10⁻⁵, r −0.26% / −0.31%, T_RH −0.06% / −0.08%, K +2×10⁻⁶ / +3×10⁻⁶, y = 4.14 / 2.62×10⁻⁴ | C |
| C12 Born φ_I → YY | not rechecked | It dominates T1, HT6t_cal1 and T1_cal1, but S/Q without it is still 73–150 | — |
| E1 induced-mixing leak | **structure stands; M_P bracket ×0.45–0.49 for TB** | The M_P bracket uses M1r's ln(M_P/M_Y) = 13.31 for every benchmark. With each benchmark's own M_Y = 30 H_i it is 9.29 (TB1) and 8.94 (TB2). R_leak becomes 8.7×10⁵ / 1.9×10⁶ (track 1.8×10⁶ / 4.2×10⁶). HT6t's near-resonant value moves ×0.57–2.2 for m_σ ± 2% (§5) | C/E |

---------------------------------------------------------------------------------------------------------------------

## 1. What was recomputed independently, and how

### 1.1 CMB numbers and reheating of the quiet T-models [C] (`sk_cmb_fluid.py`)

**Method.**
- Background in e-folds from deep slow roll; scalar and tensor MS modes (Bunch–Davies at k/aH = 300, read at 10⁻³).
- Pivot exit N_* e-folds before ε_H = 1.
- V0 is set by A_s = 2.100×10⁻⁹, and H_i := H(S1).
- **N_* from explicit e-fold bookkeeping:** a₀ in S1 units from the exact late K (ρ_r a⁴ = 3K²M_P²H_i²), with
  T₀ = 2.7255 K, g = 106.75, g_s0 = 3.931 and k_* = 0.05 Mpc⁻¹. This differs from the track's sudden ('wT') history.
- x₀ is the cycle average of ρa³ on my own oscillating T-model from S1. K comes from my own matter → radiation fluid.

| | N_* | H_i (GeV) | H_* (GeV) | n_s | r | α_s | K | T_RH (GeV) |
|---|---|---|---|---|---|---|---|---|
| TB1 mine | 52.413 | 7.457×10¹² | 8.978×10¹² | 0.96190 | 1.312×10⁻³ | −7.21×10⁻⁴ | 2000.48 | 1.997×10¹¹ |
| TB1 track | 52.345 | 7.466×10¹² | 8.989×10¹² | 0.96185 | 1.315×10⁻³ | −7.23×10⁻⁴ | 2000.47 | 1.998×10¹¹ |
| TB2 mine | 52.692 | 1.0588×10¹³ | 1.609×10¹³ | 0.96221 | 4.211×10⁻³ | −7.10×10⁻⁴ | 2589.23 | 1.262×10¹¹ |
| TB2 track | 52.607 | 1.0605×10¹³ | 1.611×10¹³ | 0.96215 | 4.224×10⁻³ | −7.12×10⁻⁴ | 2589.22 | 1.263×10¹¹ |

**Comparison with the data [R].**
- Planck 2018 (arXiv:1807.06211): n_s = 0.9649 ± 0.0042 and A_s = 2.10×10⁻⁹ are met (−0.7σ).
- BICEP/Keck 2021 (arXiv:2110.00483): r < 0.036 is met.
- ACT DR6 P-ACT-LB (n_s = 0.9743 ± 0.0034, arXiv:2503.14452): TB1 and TB2 sit at −3.6σ, as the track says.
- r ≈ 10⁻³ is within the reach of proposed σ(r) ~ 10⁻³ missions. That is an observation, not a claim.

### 1.2 K with Γ_I = 0.5 [C] (`sk_cmb_fluid.py`)

My own field ODE with the Γφ̇ friction and the radiation source runs from S1. It gives K = 12.340784 (HT6t) and
15.200283 (T1); the engine gives 12.340784 and 15.200284. With the width switched on at ε_H = 1 instead, it gives
1.5163631 and 1.8538703 in the a = 1 at ε_H = 1 normalisation, against inf_potential's 1.5163631 and 1.8538703. So
both the engine and the potential track's Python history are reproduced. The ratio is 1.101 / 1.110 (NOTES item 6).

### 1.3 Matter-era vector count for TB1 [C] (`sk_vmodes.py`)

**Method.**
- Python Stueckelberg-angle action (p = a³F²k²/(k² + a²e²F²), q = aF²k²).
- Exact massless de Sitter Goldstone with F_s at the slice; my own T-model background from S1, then an exact matter
  era.
- F² = F_v² + 2χ²/a², with the early F transient ignored, because the angle is frozen outside the horizon.
- Adiabatic read-out at ω/H = 200 (and 100).

| e | N_A(k ≤ 1, static), mine | track | ratio |
|---|---|---|---|
| 10⁻⁶ | 113.85 | 114.07 | 0.998 |
| 10⁻⁵ | 11.37 | 11.40 | 0.998 |

d ln N_A/d ln e = −1.0006. This is the matter-era crossing law I derive: N_A ∝ m_A a_*³ with a_* ∝ m_A^{−2/3}.
Combined with the e^{−1/2} deposit, the static R_vec ∝ e^{−3/2}. The track's e^{−1.42} additionally includes the e
dependence of f_eff. This validates vmodes.c plus the fluid continuation in the matter-era regime, which its M1r
(radiation) validation did not cover.

### 1.4 The ledger arithmetic [C] (`sk_ledger.py`)

From the stored counts I rebuilt R_Y, R_casc, R_YY, R_vec, W, S/Q, S/Q(e → ∞), Y/Y_Q and H_i,req for all seven rows,
with my own ⟨a_L⟩:
- radiation (πK/2Γ_L)^{1/2} for the Γ_I = 0.5 rows;
- my own matter-to-radiation fluid integral for TB.

Every S/Q agrees to ≤ 2.4×10⁻⁷, every S/Q(e → ∞) to ≤ 1.5×10⁻⁴, and Y/Y_Q and H_i,req to ≤ 6×10⁻⁸. The TB ⟨a_L⟩
differ by 2.4×10⁻⁴ because my fluid starts at S1.

### 1.5 The S/Q closed form at large W [D]

At W → ∞ the circular condensate has evaporated (the track finds this at W ≈ 5400). The hidden plasma is then an
ideal massless gas: X is complex (2 dof) and φ real (1 dof), with q = μT²/3 and s = (2π²/15)T³. The circular orbit of
λ_X|X|⁴ gives E_c = 3λA⁴ and Q = 2√(2λ)A³. Together:

S/Q → (π²/15)(30/π²)^{3/4} λ_X^{1/4} W^{3/4}/√2 = 0.190463 W^{3/4},  N_c/E^{3/4} = 2√2/(3^{3/4}λ^{1/4}) = 6.9775676.

The fit's 0.190255 W^{0.75013} matches this to 0.1% from W = 10³ to 10⁶ (ratio 0.9990–1.0007). The closed form is
therefore physically anchored where the track uses it, up to 10⁵, independently of the Hartree solver.

---------------------------------------------------------------------------------------------------------------------

## 2. The width-switch "convention" is a re-normalisation, not a yield change [C] (`sk_gamma_switch.py`)

**The claim.** NOTES §5 item 6 calls the engine's application of Γ from S1 (during the last e-fold of slow roll) "a
convention", and states "the yield moves by K^{3/2}". A perturbative width acts on oscillating quanta, not on a
slowly rolling condensate (cf. warm inflation), so switching it on at ε_H = 1 is the more physical reading.

**The test.** I re-ran the engine (a patched copy, Γ_eff = Γ for ln a ≥ 1) and compared the dense late means:

| switch at ε_H = 1 vs S1 | N_X | E_X | ε | N_YL | N_YH | N_σ | K (asympt.) | yield ∝ N_X/K^{3/2} | R_casc ∝ N_YH K^{1/2}/E_X |
|---|---|---|---|---|---|---|---|---|---|
| HT6t | −13.1% | −18.1% | +0.9% | −16.7% | −13.8% | −11.8% | −9.2% | **+0.5%** | +0.3% |
| T1 | −13.7% | −19.1% | +1.2% | −8.6% | −8.9% | −6.2% | −9.9% | **+0.9%** | +7.0% |

**Reading.**
- The extra friction delays the end of inflation by about 0.05 e-folds. Every comoving quantity therefore shifts by
  its power of a, and the physical combinations do not move.
- The yield changes by < 1%, not by K^{3/2} (−13% to −15%). S/Q changes by ≈ −1% (through ε).
- **Refutation of the NOTES remark:** the convention moves K, N_X and E_X together, so H_i,req is stable to < 1%.
- The default build reproduces the track's dense run bit for bit.

---------------------------------------------------------------------------------------------------------------------

## 3. Margins that the headline "no e_crit at any e" hides [C]

With vectors and A_L deposits removed (e → ∞), the Y_L deposit alone decides the branch:

| | R_Y (e → ∞) | R_Y needed for S/Q < 3.4716 | margin |
|---|---|---|---|
| HT6t | 248 | 4.97 | ×50 |
| T1 | 290 (+ Born) | 7.94 | ≫ |
| **TB1** | **2.08** | **0.85** | **×2.4** |
| TB2 | 4.23 | 0.77 | ×5.5 |

**For TB the margin is small.**
- In the matter era R_Y ∝ ⟨a_L⟩ ∝ Γ_L^{−2/3}. A Y_L width 3.8× larger than dq_heat's vacuum rate (1.006×10⁻⁷ H_i)
  would open a Q-matter branch for TB1 at large e (e ≳ 10⁻⁴–10⁻³). So would a 2.4× smaller deposit, for example
  partial deposition.
- "No e_crit at any e" for TB1 is therefore conditional on the Y_L ledger [Q dq_heat, E].
- The declared e = 10⁻⁶ is unaffected (S/Q = 615).

---------------------------------------------------------------------------------------------------------------------

## 4. The HT6t iteration contracts; "diverging" is not supported [C] (`sk_ledger.py`)

I used the track's own two points per family.

| | H_pot → H_req (GeV) | ln(H_req/H_pot) | contraction | secant fixed point | r there |
|---|---|---|---|---|---|
| HT6t | 1.29×10¹¹ → 2.70×10¹⁰; 2.805×10¹⁰ → 1.04×10¹⁰ | −1.565 → −0.989 | 0.63 | **≈ 2.0×10⁹ GeV** | 6.8×10⁻¹¹ |
| T1 | 1.29×10¹¹ → 9.86×10⁹; 1.043×10¹⁰ → 1.19×10¹⁰ | −2.571 → +0.136 | 0.05 | 1.18×10¹⁰ GeV | 2.3×10⁻⁹ |

- **The HT6t mismatch shrank, in both linear and log measure.** The digest's "There is no consistent point near the
  design value" is true. C10's "diverging" and the NOTES' "the single step moves away from consistency" are not.
- A consistent HT6t member may exist at about 2×10⁹ GeV. There m_I/H_i would be several hundred and
  φ_I → Y_LY_L wide open. That needs a real third step [E].

---------------------------------------------------------------------------------------------------------------------

## 5. The leak estimate: benchmark-specific logarithm [C/E] (`sk_ledger.py`)

`leak_estimate.py` brackets c with 13.31/4π² for "portals vanish at M_P". That number is ln(M_P/M_Y) at M1r's H_i.
For each benchmark's own M_Y = 30 H_i:

| | ln(M_P/M_Y) | R_leak factor | R_leak (M_P), own | track | S/Q(10⁻⁶) with it |
|---|---|---|---|---|---|
| TB1 | 9.29 | 0.488 | 8.7×10⁵ | 1.8×10⁶ | 6.4×10⁴ (track 1.1×10⁵) |
| TB2 | 8.94 | 0.451 | 1.9×10⁶ | 4.2×10⁶ | 1.2×10⁵ (2.1×10⁵) |
| HT6t | 13.35 | 1.006 | 1.32×10³ | 1.31×10³ | 302 |
| HT6t_cal1 / T1_cal1 | 14.88 / 15.87 | 1.25 / 1.42 | 19 / 11 | 15 / 7.6 | unchanged |

- The qualitative conclusion stands: the leak dominates the TB ledger unless |A_IΣ| ≲ 10⁻⁶ H_i.
- **HT6t's near-resonance** (m_I = 42.37 against m_σ = 40) makes its M_P value move ×0.57–2.2 for m_σ ± 2%.
- **Wording.** NOTES §5 says each A_L deposits "about 10³ times its production energy". For TB1 at e = 10⁻⁶ the ratio
  of m_A⟨a_A⟩ to the comoving production energy (m_I/2)a(t ≈ 1/Γ_I) is ≈ 16–20. For HT6t it is ≈ 10⁴. The
  computation uses N_I·BR·deposit directly and is not affected.

---------------------------------------------------------------------------------------------------------------------

## 6. New issues

**N1. Higgs survival during inflation is marginal at the TB Hubble scale [C on a criterion choice, E]
(`sk_higgs_TB.py`).**
- **Method.** inf_higgs's bubble-count criterion, re-implemented: tree-level RG-improved de Sitter barrier and
  ⟨N⟩ = (4π/3)∫e^{3(N₀−N_b)}(H/H_p)³e^{−B}, with the benchmark's own H(N). SM running is from smpot.py (read-only),
  with PDG 2024 M_h = 125.20.
- **Regression.** HT6t at 1.29×10¹¹ reproduces inf_higgs's ξ_min to 10⁻⁶.

| ξ_EW,min (tree) | M_t = 172.28 | 172.57 | 172.86 (+1σ) | 173.15 (+2σ) |
|---|---|---|---|---|
| HT6t (H_i = 1.29×10¹¹, N₀ = 55.76) | 0.0596 | 0.0702 | 0.0790 | 0.0866 |
| TB1 (H_i = 7.46×10¹², N₀ = 58.82) | 0.0781 | 0.0865 | **0.0938** | **0.1006** |
| TB2 (H_i = 1.06×10¹³, N₀ = 59.10) | 0.0788 | 0.0871 | **0.0944** | **0.1011** |

- **Why it is worse.** At H ~ 10¹³ GeV the barrier sits where λ is about twice as negative (λ(10H_p) ≈ −0.008 against
  −0.004), and N₀ is three e-folds larger.
- **Result.** ξ_H = 0.1 passes M_t + 1σ with only a 6–7% margin (log₁₀⟨N⟩ = −13.5 / −12.1). It **fails at +2σ**
  (log₁₀⟨N⟩ = +1.1 / +2.1), where HT6t still passes.
- The one-loop de Sitter potential lowered ξ_min by ≈ 0.013 at 1.3×10¹¹ (inf_higgs). This is not recomputed here.
- The post-inflation danger in TB's 12-e-fold matter era at H ~ 10¹³ GeV is not computed. The NOTES' "Higgs survival
  at the new H_i not recomputed" understates this, because the inflationary stage itself has lost most of its margin.

**N2. The dark-matter abundance on these benchmarks is set by hidden entropy, not by the CP-odd charge
[D on the track's p].**
- With p ≈ 0.05 S/Q, n_X̄/s_h = p/(S/Q) = 0.052–0.055 on every benchmark, so ρ_DM ≈ 0.11 m_X s_h.
- The asymmetric share of Ω_DM is 1/M: 5.8% (HT6t), 1.4% (T1), 1.5% (TB1), 0.9% (TB2), 3.0% (Send).
- So α_I and sign(φ_I), the CP bookkeeping the construction carefully preserves, control only 1–6% of the relic. The
  relic is effectively a symmetric freeze-out of a 92 MeV scalar whose amount is fixed by the hidden heat.
- The holder's role as "the source of the dark matter" (WORKING_THEORY §0) does not survive on these benchmarks. The
  H_i calibration still closes (both s_h a³ and N_X scale as H_i³).

**N3. The bottom line's "frozen pre-slice hold" condition has the wrong sign [A, with computed numbers].**
- WORKING_THEORY §4 and the dq_heat skeptic found that rolling before the slice lowers the vector heat. Q-matter at
  10⁻⁶ needs f_eff < 0.47, not a frozen (f_eff = 1) hold.
- TB2 already has f_eff = 0.468.
- What the benchmarks lack is fast reheating (K ≈ 2) with a quiet inflaton (m_I ≈ 2–4 H_i): two conditions, not three.
  TB1 at a fast width (Γ_I ≥ 10⁻⁴ H_i, K ≤ 85) already reaches S/Q(e → ∞) = 3.28 (track §4).

**N4. The width-switch convention (§2):** it changes K, N_X and E_X by 9–19% together, but the yield and S/Q by
≲ 1%. NOTES item 6 should read "re-normalisation, yield invariant".

**N5. The TB "no e_crit at any e" has a ×2.4 margin on the Y_L deposit (§3).**

**N6. Omitted effects, estimated and negligible for the chain [E] (`sk_estimates.py`).**
- **Fermion preheating at the zero crossings while φ_I → NN is blocked.**
  - The Greene–Kofman count (hep-ph/9807339) transfers 1.2×10⁻⁸ (TB1) and 3.9×10⁻⁹ (TB2) of the inflaton energy per
    crossing. That is comparable to the perturbative 2πΓ/m per period (2.1×10⁻⁸, 8.6×10⁻⁹).
  - Integrated until the blocking lifts (y\|Φ\| = m_I/2 at a ≈ 51 and 44, i.e. 3.9 and 3.8 e-folds after S1), it is
    2×10⁻⁷ and 5×10⁻⁸, so K is unchanged.
  - It does supply an early Standard Model plasma through N → LH, consistent with inf_higgs row 7. This is relevant to
    the post-inflation Higgs, which was not computed.
- **Gravitational (oscillating-R) production over the 12–13 matter e-folds.**
  - For X (coupling (2/3)R_osc, R_osc = 9H²) the injected energy is R ≈ 1×10⁻⁴.
  - Minimally coupled Goldstones number 0.03 / 0.02, i.e. R ≈ 1.0 / 1.7 against R_vec = 1796 / 3226.
- **The S1 channel's R-odd induced Aφ_I\|H\|²** at the TB couplings (seesaw y_ν² = m_νM_N/v², M_N = m_I/4) is
  0.15–3.9 PeV (TB1) and 0.09–2.4 PeV (TB2). The inflationary Higgs bound \|Aφ_I\| ≲ 0.75 H_*² on the plateau gives
  6.0 and 13 PeV. So it passes, but only by ×1.5 at the log end for TB1.
- **Fragmentation.** For α = 0.3–1 T-models self-resonance is weak (Lozanov–Amin arXiv:1608.01213) [R]. HT6t with a
  realisable slow width would oscillate long enough for its Floquet rate (μ/H ≈ 6, inf_potential skeptic) to
  fragment it. This is irrelevant to the release, which forms 0.05–0.6 H_i⁻¹ after the end [A].
- **Thermal corrections.** N and Higgs thermal masses at T_max ~ 10¹³ GeV are ≲ m_I/2, so there is no thermal blocking
  of φ → NN [E]. The hidden sector has no thermal portal.

---------------------------------------------------------------------------------------------------------------------

## 7. CP-sign bookkeeping and M1r assumptions (checked)

- **Sign.** N_X < 0 on every S1 benchmark, the same as M1r's. HT6t uses α_I = 2.192 − π because φ_I − φ_min < 0
  there. TB uses φ_I > 0 with α_I = 2.192.
- **The Send flip** (N_X > 0) survives the width-switch test.
- **TB is reflection-symmetric.** On the mirror branch (φ_I < 0) the slice has μ_eff = 2μ − 10i (\|u\| = 17.47 [Q
  inf_reheat skeptic]), so m_L² < 0 and Y condenses. That is a Send-like preparation, whose sign may flip (cf. C4).
  The sign is then a branch choice, as in M1r, but the two branches are not mirror images: μ breaks the reflection
  softly [A].
- **The prepared states stay heavy** (lightest m² = 10.15–10.19 H_i² at S1; ≥ 8.7 H² at exit [Q]). The isocurvature
  suppression assumption holds.
- **N_* and n_s** use the Planck 2018 target as the brief requires. Retrospectively, ACT DR6 puts the T-models at
  −3.6σ [R].

---------------------------------------------------------------------------------------------------------------------

## 8. Not checked

- N_X and N_Y of the TB engine runs (Γ = 0 to τ = 2×10⁴) and the a⁻² freeze-out fits.
- f_eff, the exit-factor F(s, h_r, ε) solve.
- The Born φ_I → YY counts (C12).
- The case-L p map beyond the entropy-per-particle consistency.
- Finite-k Y with the new inflatons.
- The post-inflation Higgs at TB scales.
- A lattice for T1.

## Code (`skeptic/code/`, own; ≈ 2 min on one core)

| file | purpose |
|---|---|
| sk_cmb_fluid.py | T-model background, MS spectra, N_* bookkeeping from K, x₀, decay fluid (TB1/TB2); K for HT6t/T1 with Γ from S1 or ε_H = 1 → `sk_cmb_fluid.json` |
| sk_vmodes.py | independent longitudinal-vector modes for TB1 in the matter era → `sk_vmodes.json` |
| sk_gamma_switch.py | patched engine copy (width gated at ε_H = 1), release sensitivity → `sk_gamma_switch.json` |
| sk_higgs_TB.py | Higgs Hawking–Moss bubble count with H(N) at TB scales (smpot read-only) → `sk_higgs_TB.json` |
| sk_ledger.py | ledger/S/Q/yield/H_i rebuild, closed-form asymptote, leak logs, iteration fixed points → `sk_ledger.json` |
| sk_estimates.py | fermion preheating, gravitational production, induced A-term [E] → `sk_estimates.json` |
| sk_summary.py | → `sk_summary.json` |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_chain/skeptic/code; export OMP_NUM_THREADS=1
python3 -B sk_cmb_fluid.py T K; python3 -B sk_vmodes.py; python3 -B sk_gamma_switch.py HT6t T1
python3 -B sk_higgs_TB.py; python3 -B sk_ledger.py; python3 -B sk_estimates.py; python3 -B sk_summary.py
```

The M_t = 173.15 rows of `sk_higgs_TB.json` and the HT6t_Send row of `sk_gamma_switch.json` were added by short inline
calls to the same functions.

## Credits

Planck 2018 X arXiv:1807.06211; BICEP/Keck arXiv:2110.00483; ACT DR6 arXiv:2503.14452; Kallosh–Linde arXiv:1306.5220;
Graham–Mardon–Rajendran arXiv:1504.02102; Greene–Kofman hep-ph/9807339; Lozanov–Amin arXiv:1608.01213;
Mantziris–Markkanen–Rajantie arXiv:2011.03763 (criterion form via inf_higgs) [R].

Our own tracks, read-only: inf_chain (results files; engine source for the patched copy; driver), inf_higgs
(regression values), higgs_reheat (smpot.py), inf_reheat and inf_potential skeptics (quoted points).
