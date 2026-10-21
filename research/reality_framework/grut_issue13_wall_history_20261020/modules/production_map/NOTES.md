# Module production_map: producer → (ξ, Y_N, μ_N), restoration by a fixed-number hidden bath, spectator-condensate relic

21 October 2026. Claude module of the issue #13 wall-history packet. It addresses two gaps listed in
`verification/v_critic/CRITIC.md`: R1(b) / §4 item 3 (producer and the underpopulated bath) and §4 item 6
(spectator-condensate relic). The action is **unchanged** (`BRIEF.md`). Inflation, reheating, the inflaton, its
branching into the hidden sector and the end-of-inflation Hubble rate are all **outside-action preparation data
[P]**. Nothing here validates or excludes Candidate 1. Every result is a **conditional history constraint**.

Labels: **[D]** derived here; **[C]** computed by the code here; **[Q]** quoted (BRIEF, sibling module, verifier);
**[E]** order-of-magnitude estimate; **[P]** outside-action preparation datum; **[R]** recalled from the literature,
not fetched in this session.

Verified corrections used:
- Y_N,max = **4.373e-12**, the strict bound (v_relic §2: every late quantum has m ≥ m_X).
- Resummed T_c = 1.5571–1.5641 v (v_thermal).
- ξ_kin(φ) = 6.6e-5 (thermal_potential).
- H_perc = 1617/1792 GeV (v_inflation). It is used only as context.
- inflation_spectator's (φ, |X|) Fokker–Planck variances, used as a check.

## 0. Result in brief

**A. A fixed-number (μ_N ≤ 0) hidden bath restores Z2_φ. The restored phase ends earlier, not later.**

- I solved the Hartree gap equations at φ = 0 with Bose functions, a common μ_N and quasi-particle masses
  (X massless at tree level, m_D0 = 100 GeV, φ massless on the boundary) [D, C].
  - At μ_N = 0 this gives T_c = 1.5575 v. That lies inside the resummed range.
  - At μ_N < 0 it gives T_c(z) ≈ 1.5476 v √(ζ(2)/Li₂(z)), with z = e^{μ_N/T_h}.
- For a relativistic bath with conserved yield Y_N = n_h/s and energy fraction ε = ρ_h/ρ_SM, the restored phase
  ends at T_SM,end(Y_N, ε) = G(z) Y_N^{−1/3}. For μ_N/T ≲ −2 this becomes, to < 1%, **T_SM,end = 92.3 GeV · √ε/Y_N · (v/TeV)**.
  - At fixed Y_N the exit moves to higher T_SM as μ_N decreases: an underpopulated bath stops restoring earlier.
  - The lowest exit of any relic-compatible (Y_N ≤ 4.373e-12) kinetically thermal bath is the μ_N = 0 one:
    **T_SM,min = 2.078e6 GeV (2.078–2.087e6 with resummed T_c) × (v/TeV)^{4/3}**. This agrees with relic_yield
    (2.04–2.10e6) and with v_relic (2.066–2.088e6).
- Underpopulated, relic-compatible baths do restore, if they exist above T_SM,end. But they stay **kinetically
  thermal at the exit only for μ_N/T_h ≳ −5, ξ ≲ 5e-3 and 2.08e6 ≲ T_SM,end ≲ 2.4e7 GeV**. Outside that region the
  restoration, if any, comes from a free-streaming gas, and the exit is not a thermal transition.
- The zero mode responds (Π − λ_φv² > H²) only if T_SM,end < 6.0e9 GeV.
- [E] A cold +v start keeps memory of its sign through the exit unless T_SM,end ≲ 7e8 GeV, or unless the bath damps φ.

**A2. Direct heavy-inflaton branching** [P] (B_h into the quartet, equal sharing, E = m_I/2) gives:

- **Y_N = 1.41 B_h T_RH/m_I** (C = 0.939); ε = B_h/(1 − B_h);
- formal kinetic ratio **ξ = m_I/(5.6 T_RH)**, independent of B_h;
- **μ_N/T_h ≈ −23 to −40** [C].

This bath is a sparse, hot, free-streaming gas. It is **never kinetically thermal while relativistic**: at the exit
Γ_el/H ≤ 4e-7 for every m_I ≥ 10 T_RH. Its forward-scattering curvature still restores Z2_φ [D, C] in the window

- **B_h ∈ [8.0e-17 (m_I/T_RH)², 3.1e-12 (m_I/T_RH)] (v = 1 TeV)**, for H_e = m_I/3;
- the lower edge rises to 1.1–1.4e-16 (m_I/T_RH)² for low H_e;
- the window is non-empty iff **m_I/T_RH ≲ (2–4)e4**.

Cases:

- m_I = 1e13 GeV needs T_RH ≥ 1e9 GeV, giving 8e-9 ≤ B_h ≤ 3.1e-8.
- m_I = 1e10 GeV with T_RH = 1e9 GeV gives 8e-15 ≤ B_h ≤ 3.1e-11.
- For a cold +v start, [E] memory erasure by redshift raises the lower edge to ≈ 6e-15 (m_I/T_RH)², so the window
  needs m_I/T_RH ≲ 500.
- The restored phase ends at T_SM ≈ 1e8–6e9 GeV.

This replaces relic_yield's thermal-form example "B_h ≲ 1.2e-14".

**B. Spectator relic, H_inf ≫ v, instant reheating.** In the radiation era χ = aφ obeys an autonomous quartic
equation, so the conformal energy is exactly conserved. The hidden energy fraction is then
**ε = Σ⟨φ_i²⟩/(6M_P²) + ⟨U⟩/(3H²M_P²) = 1.727 (H_inf/M_P)²** [D, C]. That is 74× the critic's (c/3)(H/M_P)².

- The term the critic missed is the frozen-field energy a'²φ₀²/2. It is 135× the potential term.
- The hidden number is bracketed by:
  - the μ_N = 0 number at that energy (floor);
  - the oscillation adiabatic invariant (no number-changing). This is 15× the floor (overpopulated).
- Result: **Y_N = (0.18 – 2.8) (H_inf/M_P)^{3/2} = (1.5e-9 – 2.4e-8)(H_inf/1e13 GeV)^{3/2}**.
- Y_N exceeds 4.373e-12 for **H_inf ≳ 2.0e11 GeV (floor) or 3.3e10 GeV (adiabatic invariant)**.
  - With running ξ and couplings at μ = H the thresholds are 1.78e11 and 2.92e10 GeV.
  - The critic's "1e-9 (H/1e13)^{3/2}, 3e11 GeV" lands near the floor because two errors partly cancel.
  - Its statement that "number freezes at amplitude ~ v" is not right. The comoving number is the adiabatic
    invariant fixed at the onset of oscillation. Fragmentation can lower it toward the floor.
- For non-instant reheating, Y_N ∝ T_RH: Y_N ≈ k H_e T_RH/M_P² with k = 1.0 (floor) or 17 (invariant). This
  requires **H_e T_RH ≲ 2.6e25 GeV² (floor), 1.5e24 GeV² (invariant)**.
- [E, R] The relic inherits an uncorrelated isocurvature perturbation from the spectator field values (power ≈ 2e-3,
  amplitude ≈ 0.045 on CMB scales). If the recalled Planck limit is applied retrospectively, the relic must be
  ≲ 2e-4 of the DM, which moves the thresholds down to about 1e8–7e8 GeV.

## 1. Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/production_map/code
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py        # about 10 min; deterministic (seed 20261021)
```

| Script | Runtime | Output |
|---|---|---|
| `a1_restoration.py` | 3 s | `results/a1_restoration.json` |
| `b1_spectator.py` | ≈ 8 min | `results/b1_spectator.json` |
| `b2_spectator_bath.py` | < 1 s | `results/b2_spectator_bath.json` |
| `b3_isocurvature.py` | 2 s | `results/b3_isocurvature.json` |
| `a2_producer.py` | 2 s | `results/a2_producer.json` (reads b1) |
| `run_all.py` | — | `results/headline.json` (assertions, sha256 of each script) |

All 12 assertions in `results/headline.json` pass. Among them: T_c(1) inside the resummed range, ξ_kin(z = 1) within 5% of 6.6e-5, the MB closed form, the single-field adiabatic-invariant calibration, energy conservation, the hidden_heating c5 cross-check, and Λ₁/√λ = 0.08892. numpy and scipy only; 2 threads.

## 2. Conventions

- Action and benchmark as in BRIEF [Q]: v = 1000 GeV, g_X = 0.1, λ_φ = 0.01, λ_X = 0.001, h_D = 0.01, λ_D = 0.001,
  m_D0 = 100 GeV, M_P = 2.435323e18 GeV.
- Real basis (φ, q₁, q₂, D), X = (q₁ + iq₂)/√2.
- SM: g* = g*s = 106.75 (constant, stated approximation). s_SM = (2π²g*/45)T³ and ρ_SM = (π²g*/30)T⁴.
- T_RH is defined by Γ_I = H(T_RH) = √(π²g*/90) T_RH²/M_P.
- Hidden dof g_h = 4.
- Relic bound: the strict Y_N ≤ 4.373e-12 [Q v_relic] (Ω h² ≤ 0.12 is a retrospective datum).
- ξ_formal = T_h/T_SM is the kinetic temperature that a bath with the given (n_h, ρ_h) would have.

---

# Part A

## 3. Part A1: restoration criterion for a kinetically thermal bath at fixed number

### 3.1 Why the number is fixed [Q, D]

In the action the hidden quartet couples to the SM only through gravity (hidden_heating: Γ₂₂/H ~ 1e-26 for the
gravitational freeze-in gas). A hidden bath able to restore Z2_φ must therefore be produced by outside-action
physics [P]. Once produced, its number changes only by 2→4 in the symmetric phase (equilibrium needs ξ ≳ 0.068,
band 0.0084–0.27, for μ_N = 0 [Q thermal_potential]) or broken-phase 3→2 (ξ > 3.48e-3 [Q relic_yield]). Both rates
scale with the fugacity (∝ z for 2→4, ∝ z² for 3→2 forward), so an underpopulated bath is even further from number
equilibrium. In the relic-compatible window ξ ≤ 7.6e-4 the hidden number is a production datum. The conserved
quantities of a relativistic hidden bath are then

- the yield **Y_N = n_h/s_SM** (comoving number / comoving SM entropy), and
- the energy fraction **ε = ρ_h/ρ_SM** (both ∝ a⁻⁴ while the bath is relativistic and g* = 106.75 is constant).

Elastic 2→2 scattering, including the conversions φφ ↔ XX̄ ↔ DD, conserves both. A kinetically thermal bath
with these two numbers is a Bose gas with kinetic temperature T_h and fugacity z = exp(μ_N/T_h) ≤ 1 fixed by
(Y_N, ε). Underpopulated means ε > ε_th(Y_N), i.e. more energy per quantum than the μ_N = 0 gas.

### 3.2 The criterion [D]

The Hartree (one-loop, self-consistent) curvature of V_eff at φ = 0 is

  ∂²V_eff/∂φ²|₀ = −λ_φ v² + Π_φ,  Π_φ = 3λ_φ⟨δφ²⟩ + g_X²(⟨q₁²⟩ + ⟨q₂²⟩) + h_D²⟨D²⟩,
  ⟨χ²⟩ = ∫ d³k f_χ(k)/((2π)³ E_χ(k)),  f = 1/(e^{(E−μ_N)/T_h} − 1).

**Restoration ⇔ Π_φ(T_h, μ_N) > λ_φ v² = 1e4 GeV² (v/TeV)²** (same Π definition as hidden_heating §1 and
inflation_spectator §4.4). The thermal check (3λ_φ + 2g_X² + h_D²)/12 = 0.004175 is reproduced.

What I do with the masses at φ = 0 (all [D]):

- **X** is massless at tree level at φ = 0. Its quasi-particle mass comes from the gap equation
  m_q² = g_X² I_φ + 4λ_X I_q (each real component; I_χ = ⟨χ²⟩).
- **D** keeps m_D0 = 100 GeV: m_D² = m_D0² + h_D² I_φ + 3λ_D I_D.
- **φ** is tachyonic at tree level (−λ_φ v²). Its quasi-particle mass² at φ = 0 is m_φ² = −λ_φ v² + Π_φ, which is
  exactly zero on the restoration boundary. I evaluate Π_φ with m_φ = 0. Since Π_φ decreases monotonically with
  m_φ, the sign of the self-consistent m_φ² is fixed by Π_φ(m_φ = 0) − λ_φ v², so this locates the boundary
  exactly at Hartree level. The tachyonic side (broken phase) is never expanded around φ = 0.
- A common μ_N for φ, q₁, q₂ (kinetic contact through g_X²). D is not in energy contact with φ for ξ ≲ 0.015
  [Q thermal_potential]; I give it the same fugacity (true for equal-sharing producers). Its weight in Π_φ is
  h_D²/(3λ_φ + 2g_X² + h_D²) = 0.2%, so this choice is immaterial.
- μ_N ≤ 0 = m_φ(boundary) ≤ m_q, m_D keeps every species below Bose condensation.
- Bose functions with quasi-particle masses are integrated by fixed Gauss–Legendre quadrature (checked to 1e-15
  against scipy `quad`). Vacuum Coleman–Weinberg terms are not included: thermal_potential's resummed schemes put
  their effect inside the quoted 1.5571–1.5641 v range.

**Results [C]** (`a1_restoration.json → Tc_of_z`):

| z = e^{μ_N/T_h} | 1 | 0.9 | 0.5 | 0.1 | 0.01 | 1e-6 |
|---|---|---|---|---|---|---|
| T_c(z)/v (Hartree) | **1.5575** | 1.7437 | 2.6021 | 6.1966 | 19.825 | 1984.9 |
| m_q/T_h at T_c | 0.0340 | 0.0303 | 0.0203 | 0.0085 | 0.0027 | 3e-5 |
| m_D/T_h at T_c | 0.0661 | 0.0591 | 0.0396 | 0.0166 | 0.0052 | 5e-5 |
| ⟨E⟩/T_h | 2.705 | 2.757 | 2.890 | 2.981 | 2.998 | 3.000 |

- T_c(1) = 1.5575 v lies inside the resummed range 1.5571–1.5641 v [Q] and equals charge_phase's T_c0 = 1.5575 v.
- As z → 0, T_c(z) → 1.547646 v √(ζ(2)/Li₂(z)) to 4e-11: the thermal masses vanish and the leading form is
  exact. thermal_potential's leading scaling (×1.681 at z = 0.5, ×4.00 at z = 0.1) becomes ×1.671 and ×3.978 with
  thermal masses.

### 3.3 In terms of Y_N and ε (or ξ) [D, C]

With (n_h, ρ_h) = (Y_N s_SM, ε ρ_SM) at the exit and Π_φ(T_h, z) = λ_φ v², scale invariance of the massless gas
gives two one-variable functions:

  ε/Y_N^{4/3} = F(z),  T_SM,end = G(z) Y_N^{−1/3},

with F(1) = 16.53 (the minimum: ε_th(Y_N) = 16.53 Y_N^{4/3}), G(1) = 339.9 GeV, and F increasing, G increasing as z
decreases (G ∝ z^{−1/6} in the MB limit). For z ≲ 0.1 (MB) the closed form is

  **T_SM,end = √(λ_φ v² ε/(K Y_N²)) = 92.33 GeV · √ε / Y_N · (v/TeV),  K = (3/2) c_eff s²/ρ_SM = 1.1730,**

with c_eff = (3λ_φ + 2g_X² + h_D²)/4 = 0.012525. The exact (Bose, Hartree) exit agrees with it to < 1e-3 for
z < 1e-6 and lies 10% below it at z = 1. Equivalently, in terms of the kinetic ratio ξ = T_h/T_SM: z follows from
n_h/T_h³ = 46.83 Y_N/ξ³, and T_SM,end = T_c(z)/ξ (table `T_SM,end(Y_N, xi)`).

Exit table [C] (`exit grid`):

| Y_N | ε/ε_th | μ_N/T_h | ξ at exit | T_SM,end [GeV] | Γ_φ,el/H at exit | A = Π/H [GeV] |
|---|---|---|---|---|---|---|
| 4.373e-12 | 1 | 0 | 7.49e-4 | **2.078e6** | 128 | 1.65e9 |
| 4.373e-12 | 2 | −1.65 | 1.37e-3 | 3.23e6 | 20.5 | 6.9e8 |
| 4.373e-12 | 10 | −6.41 | 6.8e-3 | 7.26e6 | 0.37 | 1.35e8 |
| 4.373e-12 | 1e4 | −27.1 | 6.8 | 2.30e8 | 1.2e-8 | 1.35e5 |
| 4.373e-12 | 1e6 | −41.0 | 676 | 2.30e9 | 1e-13 | 1.35e3 |
| 1e-12 | 1 | 0 | 4.58e-4 | 3.40e6 | 48 | 6.2e8 |
| 1e-13 | 1 | 0 | 2.13e-4 | 7.32e6 | 10 | 1.3e8 |
| 1e-14 | 1 | 0 | 9.9e-5 | 1.58e7 | 2.2 | 2.9e7 |

**Answer to "does an underpopulated bath restore?"** Yes. A kinetically thermal bath with any Y_N ≤ 4.373e-12 and
any μ_N ≤ 0 is in the restored phase at every T_SM > T_SM,end(Y_N, ε), provided it already exists then (production
at T_SM,prod > T_SM,end). Nothing in the criterion forbids it. What μ_N < 0 changes:

1. **The restored phase ends earlier**, at a higher T_SM: at fixed Y_N, T_SM,end ∝ z^{−1/6} (MB), or ∝ √ε at
   fixed Y_N. Relative to the μ_N = 0 bath of the same Y_N, the exit is at higher T_SM, never lower. The critic's
   phrase "restores later or not at all" should read "restores only at higher visible temperatures".
2. **Minimum visible temperature of a relic-compatible restored phase [C]:** the minimum of T_SM,end over
   Y_N ≤ 4.373e-12 and ε ≥ ε_th(Y_N) is at Y_N = Y_max, μ_N = 0:
   **T_SM,min = 2.078e6 GeV (Hartree), 2.078–2.087e6 GeV with T_c rescaled to the resummed range, ∝ (v/TeV)^{4/3}.**
   With relic_yield's thermal-composition Y = 4.30e-12: 2.090e6 (2.089–2.099e6) GeV. This agrees with relic_yield
   (2.04–2.10e6 GeV (v/TeV)^{4/3}, μ_N = 0) and v_relic route B with the strict Y (2.066–2.088e6 GeV). The μ_N = 0
   value is a strict lower bound for every kinetically thermal bath with μ_N ≤ 0. (v-scaling: Y_max ∝ 1/m_X ∝ 1/v
   and T_SM,end ∝ v Y^{−1/3} at fixed couplings.)
3. **Kinetic equilibrium at the exit is lost quickly** because the elastic rate scales with the density (∝ z).
   Using thermal_potential's rate formula (Γ_φ = n_dof Σ|M|²S_f/(16π·2⟨E⟩²), which reproduces its ξ_kin = 6.6e-5 at
   z = 1 to 0.5% [C]), a kinetically thermal exit of a relic-compatible bath needs (`kinetically thermal exit
   region`):

   | Y_N | ε/ε_th max | μ_N/T_h min | ξ_end max | T_SM,end range [GeV] |
   |---|---|---|---|---|
   | 4.373e-12 | 6.68 | −5.21 | 4.5e-3 | 2.08e6 – 5.93e6 |
   | 1e-12 | 4.51 | −4.03 | 1.9e-3 | 3.40e6 – 7.96e6 |
   | 1e-13 | 2.44 | −2.22 | 4.7e-4 | 7.32e6 – 1.26e7 |
   | 1e-14 | 1.33 | −0.58 | 1.2e-4 | 1.58e7 – 1.96e7 |

   Below Y_N ≈ 3e-15 even the μ_N = 0 bath is not kinetic at its exit (ξ < 6.6e-5). **A kinetically thermal exit
   of a relic-compatible bath therefore happens at 2.08e6 ≲ T_SM ≲ 2.4e7 GeV with μ_N/T_h ≳ −5 and ξ ≲ 5e-3.**
   Outside this region the restoration, if any, is by a free-streaming gas (same Π formula with the actual f, see
   §4), and the exit is not a thermal transition.
4. **Zero-mode response (Hubble damping) [D].** In RD with conserved (Y_N, ε) the bath is scale invariant, so
   Π_φ = A H with A constant (A = λ_φ v²/H(T_SM,end)). The zero mode can respond, Π_φ − λ_φ v² > H² at some time, iff
   **A > 2√λ_φ v = 200 GeV ⇔ T_SM,end < 5.97e9 GeV (v/TeV)^{1/2}**, and the bath must exist at H = A/2. Every
   kinetically thermal exit above satisfies this by many orders (A ≥ 1e7 GeV).
5. **Memory of a cold +v start [E].** If φ sat at +v before restoration (low-scale inflation) and the bath does
   not damp the zero mode, the coherent amplitude left at the exit is ≈ 1.7 v·(√λ_φ v/A) (WKB in m_eff² = AH − λ_φv²).
   It drops below the vacuum soft-mode level ≈ 0.11 √λ_φ v only for A ≳ 15 v, i.e. T_SM,end ≲ 6.9e8 GeV. Above that,
   the exit can inherit a single sign over the whole patch (adiabatic ζ only shifts the local clock;
   inflation_spectator §4.2) and walls need not form. The p(+) = p(−) = 1/2 statement for such a history needs
   either thermal damping Γ_φ > H or A ≳ 15 v. Not computed beyond this estimate.

### 3.4 Overpopulated bath (ε < ε_th) [D, E]

Kinetic equilibrium is then a gas at μ_N = m_min plus a condensate of the excess number. The gas alone gives
Π_φ ≥ Π_th(T_E), with T_E the μ = 0 temperature of the same energy, so the restored phase lasts at least until
T_SM = T_c(1)/ξ_E, ξ_E = (ε g*/g_h)^{1/4}. The condensate adds a positive contribution, so 2.078e6 GeV is **not** a
lower bound for an overpopulated bath. A placeholder [E]: a cold, uncharged X condensate oscillating along one
direction at φ = 0 restores with g_X²⟨q₁²⟩ = λ_φ v², i.e. a quantum density 0.040 v³, which with Y_N ≤ 4.373e-12 would
allow restoration down to T_SM ≈ 5.8e5 GeV. Such a condensate has no charge protection and fragments (X self- and
φ-induced resonance); this number is not a computed bound. (The charged analogue is relic_yield's 5.4e5 GeV
placeholder.)

## 4. Part A2: a concrete outside-action producer, heavy-inflaton branching [P]

### 4.1 Model and equations [D]

**Inputs [P]:**
- m_I and T_RH. The decay rate is Γ_I = H(T_RH).
- The branching B_h into hidden pairs φφ, XX̄ and DD, shared equally among the four real dof. XX̄ carries no net
  charge.
- Each quantum is injected with E = m_I/2.
- The end-of-inflation Hubble rate H_e. Two variants: H_e = m_I/3 (as in hidden_heating), and a low-scale
  H_e = max(1e3 GeV, 30Γ_I).
- No hidden back-reaction (B_h ≪ 1). An MD inflaton era with perturbative decay.

Comoving bookkeeping in N = ln a (code docstring of `a2_producer.py`) tracks ρ_I, ρ_r, n_h, ρ_h and
J = Σ_dof ∫d³k f/E. The source of J is 2 quanta × (2/m_I) per hidden decay, and J ∝ a⁻² while free streaming.

The hidden quanta free-stream (see 4.3), so the curvature at φ = 0 is **Π_φ = c_eff J**, c_eff = 0.012525 (same
Π formula as §3.2 with the actual f).

During MD, Π_φ = A H with **A = 24 c_eff B_h Γ_I M_P²/m_I²** [D]. The ODE reproduces this to 3e-5. In RD after
reheating A stays constant. So the Hubble-damped criterion is A > 2√λ_φ v for every T_RH (§3.3 item 4), with the
optimum at H ≈ 100 GeV. The B_h needed for that is **independent of H_e for H_e ≫ 100 GeV** (fit factor k = 1.00).
For H_e = 1 TeV, only part of the asymptotic MD population exists at H = 100 GeV (fraction 1 − (100/1000)^{1/3}),
so k = 1.8. For H_e = 30Γ_I, k = 1.4–1.5.

### 4.2 Map to (Y_N, ε, ξ, μ_N) [C]

| Quantity | Result |
|---|---|
| Y_N | **1.5 C B_h T_RH/m_I, C = 0.939** (0.920 when H_e is only 30Γ_I) |
| ε | **B_h/(1 − B_h)** exactly (both sectors sourced by Γ_Iρ_I, both radiation-like) |
| ⟨E⟩/T_SM after reheating | 0.532 m_I/T_RH |
| ξ_formal = ⟨E⟩/(3T_SM) | **m_I/(6C T_RH) = 0.178 m_I/T_RH**, independent of B_h |
| z = e^{μ_N/T_h} | π²n_h/(4T_h³) = 115.6 Y_N/ξ³ ≈ 2.9e4 B_h (T_RH/m_I)⁴, so μ_N/T_h ≈ −23 to −40 in the windows |
| Π_free/Π_thermalized at fixed (n, ρ) | 1.51: the continuous-decay spectrum has ⟨E⟩⟨1/E⟩ = 2.27 > 3/2 |

The bath is **strongly underpopulated**: ε/ε_th(Y_N) is huge whenever m_I ≫ T_RH. It is also **hotter** than the
SM (ξ_formal > 1) whenever m_I > 5.6 T_RH.

### 4.3 It never thermalizes while relativistic [C, E rates]

With thermal_potential's contact-rate formula and ⟨s⟩ = 2⟨E⟩²:

- At the exit, Γ_φ,el/H ≤ 4e-7 for every row with m_I ≥ 10 T_RH.
- Γ/H ∝ a in RD, so formal kinetic equilibrium would be reached only when the hidden kinetic temperature has fallen
  to ≤ 70 GeV, i.e. after the quanta are non-relativistic. This agrees with hidden_heating's Γ₂₂/H = 1e-9 to 1e-25
  for the branching gas [Q].
- **"ξ after kinetic thermalization" is therefore formal.** The restoration is by the free-streaming gas, whose Π
  exceeds the thermalized one by 1.51 at fixed (n, ρ).
- Number-changing rates are suppressed further (∝ z), so Y_N is frozen at production.

### 4.4 The B_h window [C]

| Edge | Formula (v = 1 TeV) | Meaning |
|---|---|---|
| lower, B_min | **8.0e-17 (m_I/T_RH)²** for H_e = m_I/3 (k = 1.00); **1.1–1.4e-16 (m_I/T_RH)²** for low H_e (k = 1.4–1.8) | zero mode responds: max_t[Π − λ_φv² − H²] > 0 |
| upper, B_max | **3.10e-12 (m_I/T_RH)** (3.17e-12 at C = 0.920) | Y_N ≤ 4.373e-12 |
| lower, [E] memory | ≈ **6e-15 (m_I/T_RH)²** (H_e = m_I/3), 6.7–8.6e-15 (m_I/T_RH)² (low H_e) | A > 15 v; needed only for a cold +v start without damping |

Comparisons:
- The no-Hubble (formal) threshold max Π > λ_φv² is much lower and depends on H_e: 1.1e-21 at
  (m_I, T_RH, H_e) = (1e10, 1e9, 3.3e9) GeV. It is not a restoration, because the zero mode is frozen (v_heating
  caveat).
- Check [C]: with Γ_I fixed at H(1e8 GeV), H_e = 1e9 GeV and m_I = 3e9 GeV this module gives 7.21996e-14
  (Hubble-damped) and 3.40697e-20 (no Hubble), against hidden_heating c5's 7.21996e-14 and 3.40691e-20 [Q].
  Agreement is to 1e-6.

**Window exists iff m_I/T_RH < 3.9e4 (H_e = m_I/3), 2.2–2.8e4 (low H_e), ∝ (v/TeV)^{−2}.** With the [E] memory
edge it needs m_I/T_RH ≲ 500.

Windows for the requested examples, H_e = m_I/3 [C] (`producer grid`). "n.h." means m_I ≤ T_RH, which is outside
the heavy-inflaton assumption; the formulas are listed but not claimed.

| m_I \ T_RH [GeV] | 1e7 | 1e8 | 1e9 | 1e10 | 1e11 | 1e12 |
|---|---|---|---|---|---|---|
| 1e10 | 8.0e-11 – 3.1e-9 | 8.0e-13 – 3.1e-10 | 8.0e-15 – 3.1e-11 | n.h. | n.h. | n.h. |
| 1e11 | 8.0e-9 – 3.1e-8 | 8.0e-11 – 3.1e-9 | 8.0e-13 – 3.1e-10 | 8.0e-15 – 3.1e-11 | n.h. | n.h. |
| 1e12 | **none** | 8.0e-9 – 3.1e-8 | 8.0e-11 – 3.1e-9 | 8.0e-13 – 3.1e-10 | 8.0e-15 – 3.1e-11 | n.h. |
| 1e13 | **none** | **none** | 8.0e-9 – 3.1e-8 | 8.0e-11 – 3.1e-9 | 8.0e-13 – 3.1e-10 | 8.0e-15 – 3.1e-11 † |

† With H_e = m_I/3 = 3.3e12 GeV and T_RH = 1e12 GeV, the spectator relic of Part B is 0.14 (floor) to 2.3 (adiabatic
invariant) times Y_max. That cell can be relic-excluded by the inflationary spectators alone. With H_e = m_I/3 the
cell (1e13, 1e11) carries 0.014–0.24 Y_max; every other cell carries ≤ 0.024 Y_max. In the low-H_e variant every
cell carries ≤ 3e-5 Y_max.

**State at the window edges [C]:**

| Quantity | B_h = B_min | B_h = B_max |
|---|---|---|
| exit T_SM,end | 1.9e8 GeV (T_RH = 1e7) to 6.0e9 GeV (T_RH ≥ 1e11) | 3e6–1.5e9 GeV |
| Γ_el/H at the exit | 2e-14 – 2e-9 | ≤ 4e-7 (heavy rows) |
| μ_N/T_h at the exit | — | −23 to −40 |
| ξ_formal at the exit | — | 1.8 – 2300 |

**Validity flags.**
- With H_e = m_I/3 the maximum SM temperature during reheating exceeds m_I/2 in most cells (T_max = 1.3e10–2.4e13
  GeV). The low-H_e variant has T_max = 3.2e8–1.4e12 GeV, and is valid in every cell with m_I ≥ 10 T_RH.
- Thermal blocking or plasma corrections to the SM channels (which would change the effective B_h) are not modelled.
  The hidden-sector map itself is unaffected.

**Relevance to walls.**
- For H_e ≳ 1.4–1.8 TeV the inflationary spectator already leaves both signs of φ (inflation_spectator regime D
  [Q]). The producer's restoration then re-forms the walls at its exit; it does not create the possibility.
- For low-scale inflation (cold +v start, regimes A–C) the producer is what restores the symmetry. p(+) = p(−) at
  its exit then additionally needs the memory criterion [E] or thermal damping. Neither is supplied by this
  free-streaming gas at B_min.

---

# Part B

## 5. Spectator-condensate relic for H_inf ≫ v (instant reheating) [P for the inflationary inputs]

### 5.1 End-of-inflation state [C, Q]

The Starobinsky–Yokoyama equilibrium P ∝ exp(−8π²U/(3H⁴)) [R] is computed for the full quartic U on a
140×220×220 grid in (φ, |X|, D):

| Variant | ⟨φ²⟩/H² | ⟨q₁²+q₂²⟩/H² | ⟨D²⟩/H² | Σ⟨χ²⟩/H² | ⟨U₄⟩/H⁴ |
|---|---|---|---|---|---|
| μ0 couplings, ξ = 0 (exact quartic, H-independent) | 0.681 | 5.458 | 4.148 | **10.29** | 0.0380 |
| couplings and Hubble masses 12ξ_a H² at μ = H [Q inflation_spectator], H = 1e10 GeV | 0.868 | 7.194 | 4.561 | 12.62 | 0.0567 |
| same, H = 1e13 GeV | 0.956 | 8.023 | 4.730 | 13.71 | 0.0666 |

Checks:
- inflation_spectator's 2D FP at 1e13 GeV gives 0.97 and 8.0 [Q].
- The single-field quartic gives 1.3176/√(λ/0.01) analytically, and D alone gives 4.167.
- hidden_heating's c = 0.042–0.094 [Q] is the same order as ⟨U⟩/H⁴ = 0.038–0.067.
- Reaching equilibrium needs ≳ 100–350 e-folds at H_inf (1/Λ₁ for λ = 0.01 and 0.001). A 60-e-fold-from-zero
  variant is given in 5.3.

### 5.2 Energy fraction: exact in RD [D, C]

After instant reheating a = a_e(1 + a_e H_e τ) and a'' = 0. With χ_i = aφ_i,

  χ_i'' − ∇²χ_i + ∂V₄(χ)/∂χ_i = 0 (the ξR term vanishes because R = 0; the (1 − 6ξ)a''/a term vanishes),

as long as |φ| ≫ v, m_D0. V₄ is a homogeneous quartic, so the conformal system is autonomous and
E_c = χ'²/2 + V₄(χ) is **exactly conserved**. The initial data are χ = a_eφ₀ and χ' = a_e²(H_eφ₀ + φ̇₀). They give the
late-time ratio

  **ε = ρ_h/ρ_SM = Σ_i⟨φ_{0,i}²⟩/(6M_P²) + ⟨U⟩/(3H_e²M_P²) + O(slow roll)**.

| | coefficient c_ε in ε = c_ε (H/M_P)² |
|---|---|
| analytic, minimal couplings | 1.7272 |
| sampled initial data incl. slow-roll velocity [C] | **1.7268** |
| critic's (c/3), c = 0.07 [Q CRITIC] | 0.0233 |

The frozen-field term is 135× the potential term. Physically, the frozen field keeps its potential energy while the
radiation redshifts, until H falls to the oscillation frequency. The conformal variables make this exact without
specifying an onset time.

### 5.3 Number [C]

The homogeneous conformal 3-field system (φ, signed |X| ray, D) is integrated by velocity Verlet: 1500 samples from
P_eq, dt = 0.02, τ ≤ 3000, averages over τ ∈ [1000, 3000]. The energy drift is 8e-7. Since χ' ∥ χ initially, the X
charge is exactly zero and X stays on a ray.

The number estimator is J = Σ_i⟨χ_i'²⟩/ω_i, with ω_i from zero crossings. For one quartic field it equals the exact
adiabatic invariant (4/3)E/ω, ω = 0.8472√λ amplitude, to ≤ 0.2% (calibration). The curvature-frequency estimator
would give 0.72× instead.

**Results:**
- J = 20.9 ± 0.3 H_e³ (comoving, a_e = 1).
- Shares: φ 11%, X 42%, D 47%.
- Mean quantum energy E_c/J = 0.25 H_e (critic: 0.08 H_e).

**Bracket for the final hidden number.** The coherent oscillation fragments: X preheating has q = 1, Floquet
μ_max = 0.147 [Q inflation_spectator], and φ and D self-resonate. A turbulent cascade can then lower the number
before the occupations reach O(1) and number-changing freezes (thermal 2→4 needs ξ ≳ 0.068).

- **Floor**: the μ_N = 0 number of a bath with energy ε. ξ_E = 2.605 (H/M_P)^{1/2}, Y_floor = 0.184 (H/M_P)^{3/2}.
- **Upper**: J conserved. **Y_J = 2.83 (H/M_P)^{3/2}**. J/n_eq(ρ) = **15.4**, i.e. overpopulated.

The cascade itself was not simulated.

| | Y_N at 1e13 GeV | threshold H_inf for Y_N = 4.373e-12 |
|---|---|---|
| floor (μ_N = 0 at the same energy) | 1.53e-9 | **2.01e11 GeV** |
| adiabatic invariant J | 2.35e-8 | **3.26e10 GeV** |
| critic's formula (c/3)(3/4)(T_RH/0.08H), re-derived [Q] | 9.8e-10 | 2.70e11 GeV (2.67e11 with Y = 4.30e-12) |
| running ξ and couplings at μ = H | — | 1.78e11 (floor), 2.92e10 (J) GeV |
| 60 e-folds from zero (Gaussian, N/(4π²) per real field) | — | 2.54e11 (floor), 4.51e10 (J) GeV |

Equivalent statement for the floor: the summed spectator variance must satisfy **√(Σ⟨φ_i²⟩) ≲ 6.5e11 GeV** at the
end of inflation, independent of the quartic couplings.

**Check of the critic's derivation.**

1. Its energy fraction (c/3)(H/M_P)² is the end-of-inflation ratio. It omits the frozen-field energy a'²φ₀²/2,
   which dominates by 135×, so the true ε is 74× larger.
2. Its quantum energy 0.08 H is 3.1× too small.
3. "Quartic oscillation redshifts like radiation until amplitude ~ v; then number freezes" is not correct as stated.
   The comoving number of the coherent oscillation is the adiabatic invariant, fixed at oscillation onset (the
   conformal system is autonomous). The later quartic-to-quadratic passage at amplitude ~ v is adiabatic and leaves
   it unchanged. Only fragmentation and rescattering change it, toward the floor.

Errors 1 and 2 partly cancel. The critic's 1e-9 (H/1e13)^{3/2} lies 1.6× below the floor and 24× below the
invariant. The threshold "H_inf ≳ 3e11 GeV" should read **3.3e10–2.0e11 GeV**.

### 5.4 Non-instant reheating [C]

The same 300-sample system is run in an MD + perturbative-decay background (a''/a = a²ρ_I/(6M_P²)) for
γ = Γ_I/H_e = 10 … 1e-4. The ratio R = Y(γ)/Y(instant) at fixed H_e is:

| γ | 10 | 1 | 0.1 | 0.01 | 1e-3 | 1e-4 |
|---|---|---|---|---|---|---|
| T_RH/T_RH,inst | 3.16 | 1 | 0.316 | 0.1 | 0.032 | 0.01 |
| R (J) | 1.00 | 0.994 | 0.805 | 0.332 | 0.107 | 0.0338 |
| R (floor) | 1.00 | 0.99 | 0.765 | 0.300 | 0.097 | 0.0307 |

The fit for γ ≤ 1e-2 is R = A γ^p with p = 0.496 (derived: 1/2, i.e. Y ∝ T_RH), A = 3.27 (J) and 2.94 (floor). So for
T_RH ≲ 0.1 T_RH,inst:

  **Y_N ≈ k H_e T_RH/M_P², k = 1.0 (floor), 17 (J); Y_N ≤ 4.373e-12 ⇔ H_e T_RH ≲ 2.6e25 GeV² (floor), 1.5e24 GeV² (J)**.

H_e,max(T_RH) [C] (`spectator bound H_e,max(T_RH)`, grid up to 1e16 GeV):

| T_RH [GeV] | ≤ 1e8 | 1e9 | 1e10 | 1e11 | 1e12 |
|---|---|---|---|---|---|
| H_e,max (floor) | > 1e16 | > 1e16 | 2.2e15 | 2.3e14 | 2.3e13 |
| H_e,max (J) | > 1e16 | 1.3e15 | 1.3e14 | 1.4e13 | 1.4e12 |

For T_RH ≲ 1e9 GeV the spectator abundance alone does not constrain any H_e of interest.

### 5.5 Is the spectator bath a restoring bath? [D, C]

Its energy fraction makes it a hidden bath with ξ_E = 2.605 (H_inf/M_P)^{1/2} and J/n_eq = 15. The gas part alone
restores (§3.4) until **T_SM = T_c/ξ_E ≈ 9.3e6 GeV (H_inf/1e10 GeV)^{−1/2}**; the condensate excess prolongs it.

| H_inf | Effect |
|---|---|
| ≳ 2.4e4 GeV | zero mode responds |
| ≳ 1.6e9 GeV | a μ = 0-thermalized version is kinetically thermal at its exit |

Two cases:
- For H_inf ≳ 2.4e4 GeV the spectator bath itself restores Z2_φ after inflation (gas bound), and walls re-form at
  its exit with p = 1/2, because P_eq is P-symmetric.
- For 1.6–1.8 TeV ≲ H_inf ≲ 2.4e4 GeV inflation alone already leaves p_- = 1/2 (inflation_spectator regime D [Q]).

Either way the network-forming branch needs no further hidden producer. **The relic bound caps these histories at
H_inf ≲ 3e10–2e11 GeV (instant reheating).**

### 5.6 Isocurvature of the spectator relic [E, R, retrospective]

The relic density in each patch is set by the local φ₀. It therefore carries an isocurvature perturbation
uncorrelated with ζ.

**Method.** Single-field Fokker–Planck spectral expansion per species, with f = |φ|^{3/2} and
Δ²_f(k) = Σ c_n² 2Λ_n (k/a_eH)^{2Λ_n} [R: Starobinsky–Yokoyama; Markkanen et al.]. The code reproduces
Λ₁/√λ = 0.08892. At N_* = 50–60 it gives P_S ≈ (1.8–2.2)e-3.

**Bound.** I take the Planck 2018 limit on uncorrelated CDM isocurvature as β_iso ≲ 0.04 with P_ζ = 2.1e-9 [R, not
fetched]. Then the relic can be at most a fraction f ≲ 2e-4 of the dark matter, i.e. Y_N ≲ 9e-16. That moves the
thresholds to **H_inf ≲ 7e8 GeV (floor), 1.1e8 GeV (J)**.

This is an order-of-magnitude estimate: the φ–X coupling is ignored in the spectra and the Planck number is
recalled. It does not replace the abundance threshold. It flags that the abundance bound is probably not the
strongest one.

---

## 6. Conditional history constraints (no exclusion, no validation)

**C-A1 (any kinetically thermal hidden bath with μ_N ≤ 0, unchanged action):**

1. Z2_φ is restored while T_SM > T_SM,end(Y_N, ε). In the MB regime T_SM,end = 92.3 GeV √ε/Y_N (v/TeV). Near
   μ_N = 0 use the table.
2. With Y_N ≤ 4.373e-12, every such restored phase ends at T_SM ≥ 2.08e6 GeV (v/TeV)^{4/3}. It ends at that value
   only for μ_N = 0 and Y_N = Y_max.
3. A thermal (kinetically equilibrated) exit additionally needs μ_N/T_h ≳ −5 and T_SM,end ≲ 2.4e7 GeV. The zero-mode
   response needs T_SM,end < 6e9 GeV.
4. T_c = 1.557–1.564 v holds only at μ_N = 0. For μ_N < 0 the kinetic T_c is 1.5476 v √(ζ(2)/Li₂(z)) × (1 + ≤ 0.7%).

**C-A2 (heavy-inflaton branching [P]):**

1. Z2-restoring and relic-compatible only for 8e-17 (m_I/T_RH)² ≲ B_h ≤ 3.1e-12 (m_I/T_RH) (v = 1 TeV), which needs
   m_I/T_RH ≲ (2–4)e4.
2. The resulting bath is non-thermal, with μ_N/T ≈ −23 to −40 and ξ_formal ≈ 0.18 m_I/T_RH. Its restored phase ends
   at T_SM ≈ 1e8–6e9 GeV.
3. For a cold +v start, p(+) = p(−) at the exit is not established unless B_h ≳ 6e-15 (m_I/T_RH)² [E] or the
   φ condensate is damped.
4. The end-of-inflation H_e must also satisfy the C-B bounds.

**C-B (spectator relic, minimal coupling, equilibrium variances [P]):**

1. With instant reheating Y_N ≤ 4.373e-12 requires H_inf ≲ 3.3e10 GeV if the oscillation number survives
   fragmentation, or ≲ 2.0e11 GeV if it relaxes to the μ_N = 0 number.
2. With MD reheating it requires H_e T_RH ≲ 1.5e24–2.6e25 GeV².
3. [E, R] Isocurvature, if applied retrospectively, gives H_inf ≲ 1e8–7e8 GeV.
4. Above these, such a history is inconsistent with the known DM density (retrospective) unless outside-action
   entropy injection dilutes it by Δ ≥ Y_N/Y_max.

## 7. Free inputs

- **[P]** Hidden production in Part A: Y_N, ε (equivalently ξ and μ_N), the production time, the D population,
  and the pre-restoration φ state (cold +v or random).
- **[P]** The heavy-inflaton producer: m_I, T_RH (Γ_I), B_h, equal sharing, the injection energy m_I/2, H_e, and
  perturbative decay.
- **[P]** Inflation in Part B:
  - H_inf (= H_e) and instant vs MD reheating;
  - the number of e-folds (equilibrium vs 60 from zero);
  - N_* for isocurvature;
  - the minimal coupling, with the running ξ variant quoted.
- **Retrospective data:** ω_DM = 0.12 through Y_max [Q]; the Planck isocurvature limit [R].
- **Scale v:** T_SM,min ∝ v^{4/3}; T_SM,end ∝ v at fixed (Y_N, ε); B_min ∝ v; B_max ∝ 1/v; H_inf thresholds
  ∝ Y_max^{2/3} ∝ v^{−2/3}.

## 8. What failed or was corrected

1. **Convention mismatch with hidden_heating.** My first comparison with c5 differed by 2.44. hidden_heating fixes
   Γ_I = H(1e8 GeV) and reports a ρ_r-based T_RH = 6.4e7 GeV for the same run; I had used 6.4e7 as T_RH. With the
   same Γ_I the two codes agree to 1e-6. The remaining difference is the T_RH convention only.
2. **My own first Part B estimate** multiplied the critic's ratio by (H_e/H_osc)² with a guessed H_osc ≈ 0.1 H_e.
   It is superseded by the exact conformal energy (§5.2), which needs no onset time.
3. **The first running-ξ rows** put the inflationary ξR potential into the post-inflation conformal energy. In RD,
   R = 0, so only U₄ enters. Corrected; the effect is < 1%.
4. **The number estimator is not unique for the coupled system.** The zero-crossing and curvature-frequency
   estimators differ by 25% (20.9 vs 15.7). I use the zero-crossing one, which is exact for one field. The
   difference is a systematic on the J edge of the bracket.
5. **No-Hubble restoration criterion.** It is not a physical restoration when H > √Π (frozen zero mode). It is kept
   only as a formal comparison number.

## 9. Not done, and limits

**Part A**
- One-loop Hartree only, without Coleman–Weinberg terms. No Ginzburg-region treatment near T_c; thermal_potential's
  resummed range covers the scheme spread at μ_N = 0.
- Rates are thermal_potential's order-of-magnitude contact formula, without Bose enhancement, collinear logs or
  thermal masses in the rates.
- The overpopulated (condensate) case is bounded only from one side. The 5.8e5 GeV uncharged-condensate number is a
  placeholder [E].
- The memory criterion is a WKB estimate [E]. Damping of the φ zero mode by a sparse gas is not computed.
- Producer: perturbative decay only (no preheating, no thermal blocking), equal sharing, injected energy exactly
  m_I/2, masses neglected in J (E ≫ m), g* constant.

**Part B**
- No lattice simulation of fragmentation or the turbulent cascade, so the final number is only bracketed (factor 15).
- Homogeneous-patch dynamics; gradient energy of the frozen superhorizon field is neglected.
- Equilibrium variances assume long inflation.
- Gravitational particle production of the quartet at the end of inflation and inflaton-condensate scattering
  (hidden_heating c4: n ~ 8e-3 H_e³ with E ~ m_I) are ≲ 1e-3 of J and were not added.
- The isocurvature estimate is [E, R].

**Neither part**
- No wall persistence, network or observable. Those belong to other modules.

## 10. Literature (no text copied; nothing fetched in this session)

| Reference | Read level | Used for |
|---|---|---|
| Starobinsky & Yokoyama, PRD 50 (1994) 6357, astro-ph/9407016 | recalled [R] | equilibrium P_eq, Λ₁ = 0.08892√λ (reproduced numerically) |
| Markkanen, Rajantie, Stopyra, Tenkanen, JCAP 08 (2019) 001, arXiv:1904.11917 | recalled [R] | spectral expansion of correlators |
| Greene, Kofman, Linde, Starobinsky, PRD 56 (1997) 6175, hep-ph/9705347 | via inflation_spectator [Q] | Lamé resonance, q = 1 band |
| Planck 2018 X (inflation), arXiv:1807.06211 | recalled [R] | uncorrelated CDM isocurvature β_iso ≲ 0.04 (retrospective, approximate) |
| Micha & Tkachev, PRD 70 (2004) 043538, hep-ph/0403101 | recalled [R] | turbulent cascade lowers the particle number (qualitative only) |
| Sibling modules and verifiers (BRIEF, CRITIC, relic_yield, v_relic, thermal_potential, hidden_heating, inflation_spectator, v_heating, v_inflation) | read | quoted inputs [Q] |

## 11. Files

| File | Contents |
|---|---|
| `code/common.py` | constants, Bose integrals, Hartree gap equations at φ = 0, T_c(z), elastic rate |
| `code/a1_restoration.py` | Part A1 → `results/a1_restoration.json` |
| `code/b1_spectator.py` | Part B moments, conformal dynamics, MD reheating → `results/b1_spectator.json` |
| `code/b2_spectator_bath.py` | spectator bath as a restoring bath → `results/b2_spectator_bath.json` |
| `code/b3_isocurvature.py` | isocurvature estimate → `results/b3_isocurvature.json` |
| `code/a2_producer.py` | heavy-inflaton producer → `results/a2_producer.json` |
| `code/run_all.py` | runs everything, `results/headline.json` |
