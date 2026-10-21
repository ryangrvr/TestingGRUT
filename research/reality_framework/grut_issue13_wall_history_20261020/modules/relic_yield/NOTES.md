# Module relic_yield: conserved hidden charge, stable relics of a hot hidden sector, and FLS prior art

21 October 2026. Claude module of the issue #13 wall-history packet. The action is unchanged
(`BRIEF.md`; `grut_candidate_1_20261018/GRUT_CANDIDATE_1.md` sections 2-3). This module adds no bias, portal or
inflaton coupling. Inflation, reheating, hidden-sector production and charge generation are **outside-action
preparation data**, and they are labelled where they enter. Nothing here validates or excludes Candidate 1. The
results are conditional history constraints. Known cosmological data are used retrospectively, and that use is
labelled.

Labels: **[D]** derived in closed form here; **[C]** computed by the code here; **[Q]** quoted (literature, the
repository, Codex commits or sibling modules); **[E]** order-of-magnitude estimate; **[P]** outside-action
preparation input; **[R]** recalled, not fetched in this session.

## 0. Result in one paragraph

The comoving net X charge is exactly conserved by the unchanged action. This follows from Noether's theorem, there
is no anomaly, and an FRW lattice check holds it to 4e-16 [D, C]. Candidate 1 therefore cannot generate Y_Q. Y_Q is
an outside-action initial datum supplied at a free time T_Q. It cannot predate the last ≈14–20 inflationary e-folds
for r = 1 [C, P].

In the +v vacuum the action has three stable species: X (by U(1)_X), D (by Z2_D) and h. h is stable because
h → XX̄ and h → DD are closed and it has no portal. Its gravitational lifetime is estimated at > 1e38 s, against an
age of 4.4e17 s [D, E]. A hidden bath that was ever relativistic therefore leaves a symmetric relic, and that relic
is not removed:

- the forbidden channel XX̄ → hh only converts quanta into the stable, secluded h;
- the broken-phase 3→2 (cannibal) processes reach equilibrium only for ξ > 3.5e-3 [C].

With a μ = 0 thermal reference number, Ω_sym h² = 0.272 (ξ/10⁻³)³ near the bound. The retrospective bound
Ω_sym h² < 0.12 gives **ξ_max = 7.58e-4** (8.36e-4 if D was never populated) [C]. A Z2-restoring thermal history then
needs the hidden bath above T_c at **T_SM ≥ 2.04–2.10e6 GeV (v/TeV)^{4/3}** [C]. That is possible: a kinetic
window 6.6e-5 ≲ ξ ≤ 7.6e-4 exists. The bath can never be in number equilibrium there, so its number chemical
potential is a production input. For r ≳ 1e-3 the hidden charge in such a bath is Bose-condensed at T_c and keeps
φ = 0 after the thermal restoration ends. The exit is then the charge_phase spinodal route, with p(+) = p(−) [C].

None of the eight FLS/soliton-cosmology papers checked analyses an exact-Z2 wall network. Each uses a potential
without Z2 or an explicit bias. FGGK 1988 and Macpherson–Campbell 1995 say explicitly that the exact-Z2 limit gives
the domain-wall problem [Q].

## 1. Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/relic_yield/code
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py      # 200 s total, deterministic seeds
```

`run_all.py` runs five scripts in order and writes `results/headline.json`, which holds the headline numbers,
24 assertions and the sha256 of every script.

| Script | Runtime | Output |
|---|---|---|
| `r1_charge.py` | 2 s | `charge_conservation.json` |
| `r2_spectrum.py` | < 1 s | `spectrum_stability.json` |
| `r3_rates.py` | 114 s | `rates.json` |
| `r4_boltzmann.py` | 84 s | `boltzmann.json` |
| `r5_charge_history.py` | < 1 s | `charge_history.json` |

**All 24 assertions pass.** Only numpy and scipy are used.

## 2. Conventions

- **Action.** Benchmark and action as in `BRIEF.md` [Q]: v = 1000 GeV, g_X = 0.1, λ_φ = 0.01, λ_X = 0.001,
  h_D = 0.01, λ_D = 0.001, m_D0 = 100 GeV and M_P = 2.435323e18 GeV. The |X|²D² portal α is zero at μ0. Couplings
  are fixed at their μ0 values, with no running inside the rates.
- **Fields and metric.** Real basis (h, q1, q2, D) with φ = v + h and X = (q1 + i q2)/√2. The action uses
  (−,+,+,+). Amplitudes use mostly-minus internally, and every output is metric independent.
- **Cosmological inputs [Q].**
  - s0 = 2891.2 cm⁻³ and ρ_c/h² = 1.05375e-5 GeV cm⁻³, so Ω h² = m·Y·2.7437e8 GeV⁻¹.
  - ω_b = 0.02233 gives **Y_B = 8.674e-11**.
  - g*s = g* = 106.75, held constant (stated approximation).
  - ω_DM = 0.120 ± 0.001: Planck 2018, abstract read. It is a **retrospective known datum**, and it applies only if
    the relic behaves as CDM.
- **The hidden temperature ratio ξ.** ξ ≡ T_h/T_SM, defined when the hidden sector is still relativistic: at
  T_h = 4 m_X with μ = 0, the "thermal reference". Between T_c ≈ 1.56 TeV and 400 GeV ξ is constant to ~1%,
  because g_h = 4 and the SM has g*s = 106.75 at T_SM ≥ 400/ξ GeV ≫ 200 GeV, with no number-changing, as shown
  below.

## 3. Item 1: exact U(1)_X, conservation, anomaly, chronology

**Noether [D].** U depends on X only through |X|². Under X → e^{iα}X the action is invariant. With the current
j^μ = i(X ∇^μX* − X* ∇^μX) and the field equation □X = ∂U/∂X* = U′(|X|²)X, the divergence is

∇_μ j^μ = i(X □X* − X* □X) = i(|X|²U′ − |X|²U′) = 0 on shell.

In FRW, d(a³ n̄_Q)/dt = −(volume average of a spatial divergence) = 0. So the comoving net charge N_Q = a³ n̄_Q is
constant, and Y_Q = n_Q/s changes only through entropy production.

- **Algebraic check [C].** q1 ∂U/∂q2 − q2 ∂U/∂q1 = 0 to 1.3e-16 relative, at 2000 random field points. The check
  includes the U(1)_X-symmetric counterterms m_X0²|X|² and α|X|²D², which EFT rule 2 admits.
- **Lattice FRW check [C].**
  - Setup: periodic 16³ lattice; all four fields with the full potential; random inhomogeneous initial data with
    net charge; radiation era a ∝ t^{1/2}, with a growing by 3.16 over the run.
  - The kick-drift-kick leapfrog conserves the comoving Q = Σ(q2 π1 − q1 π2) dx³ to **4.4e-16**. Each kick and
    drift conserves it exactly for a U(1)-invariant force and lattice Laplacian.
  - The physical density falls as a⁻³ to 1e-13.
  - The fields stay finite (max |field| is recorded).

**No anomaly [D].** U(1)_X is global and acts only on the gauge-singlet scalar X. No fermion carries X charge.

- There is therefore no ABJ anomaly, and no mixed U(1)_X–gravitational anomaly, since both need chiral fermions.
- No gauge field couples to j_X, so there is no instanton or sphaleron violation.
- EFT rule 2 (GRUT_CANDIDATE_1.md §3) admits only symmetric counterterms, so running cannot generate violation.
- Gravity-only production is C-symmetric, so it gives zero net charge (hidden_heating module [Q]).

**Verdict:** Candidate 1's dynamics cannot generate a hidden yield. Global-charge violation by quantum gravity
(black holes, wormholes, Planck-suppressed charged operators) is **not in the action**. Adding it would be an action
extension; postulate P4 claims no quantum-gravity protection.

**Inflationary dilution [C, P].** Charge present N e-folds before the end of inflation is diluted by e^{−3N}.

Requiring the charge to be subdominant at that time (the bound ρ_Q,i ≥ m*|n_Q,i| is the classical floor [Q]), with
instantaneous reheating as the most permissive case, gives

N < N_max = (1/3) ln[3H²M_P² / (m* r Y_B s(T_RH))].

| H_inf [GeV] [P] | 1e2 | 1e3 | 1e6 | 1e9 | 1e13 |
|---|---|---|---|---|---|
| N_max, r = 1 | 13.9 | 14.3 | 15.4 | 16.6 | 18.1 |
| N_max, r = 1e-3 | 16.2 | 16.6 | 17.7 | 18.9 | 20.4 |

Observable scales leave the horizon ≈ 50–60 e-folds before the end. A pre-inflationary charge is therefore excluded
as the source of Y_Q.

**Chronology:** Y_Q must be supplied by outside-action physics during the last ≲ 20 e-folds, at reheating or later.
The following are all **free inputs**, none fixed by the action:

- the supply temperature T_Q;
- the spatial form of the supplied charge (homogeneous condensate, plasma, or localized);
- any later entropy dilution.

## 4. Item 2a: spectrum and stability in the +v vacuum

**Masses [D].** m_X = 100, m_h = 141.421 and m_D = 100.499 GeV.

**Broken-phase vertices [D].** These come from expanding U at φ = v + h:

- λ_φ v h³ (Feynman vertex 6λ_φ v = 60 GeV) and λ_φ/4 h⁴;
- 2g_X² v h|X|² (20 GeV) and g_X² h²|X|²;
- h_D² v h D² (vertex 0.2 GeV) and ½h_D² h²D²;
- λ_X|X|⁴ and λ_D/4 D⁴.

The symmetric phase has **no cubic vertex**, so there particle number changes only by 2→4 (thermal_potential [Q]).

| Species | Verdict | Reason |
|---|---|---|
| X, X̄ | exactly stable | lightest U(1)_X-charged state |
| D | exactly stable | lightest Z2_D-odd state |
| h | no non-gravitational decay | h → XX̄ needs 200 GeV and h → DD needs 201.0 GeV, but m_h = 141.4 GeV; portals are zero at μ0 and stay zero under non-gravitational running [Q hidden_heating] |

**Gravitational decay of h** (estimates):

| Channel | Formula | τ_h |
|---|---|---|
| Naive dimensional (upper rate) | Γ = m_h⁵/M_P⁴ = 1.6e-63 GeV | **4.1e38 s** |
| h → graviton pair, loop [E] | L = (c/v) h R_abcd R^abcd with c = Σ_i n_i v ∂ln m_i²/∂h /(32π²·180) = 1.23e-4 (X: 4, D: 0.020, h self: 3); Γ = c² m_h⁷/(4π v² M_P⁴) | 1.7e49 s (band 1.7e48–1.7e50) |
| h → SM via the trace, nonminimal ξ_φ [E] | −ξ v h R → (ξ v/M_P²) h T: a dilaton-like coupling κ = ξ v v_H/M_P²; Γ ≈ κ² Γ_SM-Higgs(141 GeV) × (1–10) | 4e45–4e50 s for \|ξ_eff\| = 1e-3–1e-5 |

Notes on the table:

- Tree-level h → gg vanishes for minimal coupling. The potential gives no vertex linear in h (U′(v) = 0), and h R
  vanishes on shell.
- The loop row uses Γ_SM-Higgs(141 GeV) ≈ 8.6e-3 GeV [R]. The low-energy theorem is used with
  m_h²/4m_X² = 0.5, hence the ×10 band.
- The ξ row: ξ_φ(μ0) = 0 by matching and ξ_φ(m_h) = −3.7e-5 from running [Q hidden_heating β].

**Verdict:** h is a cosmologically stable relic of the action. All in-action estimates exceed the age by more than
1e20.

**Outside-action comparison.** None of these operators is in the action; they are shown for scale only:

- a Z2_φ-odd (φ/4M_P) F F would give τ ≈ 2.8e8 s, and it would also bias the walls;
- φ²FF/(4M_P²) gives 4.1e38 s;
- a graviton-induced Higgs portal δλ ≈ 3.4e-34 [Q] gives 2.0e41 s.

## 5. Item 2b: processes allowed at tree level and their rates

**Method [D, C].** Tree amplitudes are built by Berends–Giele recursion over the exact cubic and quartic tensors
∂³U and ∂⁴U at the vacuum. Each diagram contributes M = −Πλ Π 1/(p² − m²). The recursion is checked as follows
(all in `rates.json: checks`):

- the analytic tensors agree with finite differences to 1e-12;
- hh → XX̄ agrees with the hand-written s/t/u/contact amplitude to 1e-15;
- charge-violating (XX → hh) and Z2_D-violating amplitudes are exactly 0;
- Lorentz invariance holds to 6e-16 and identical-particle symmetry to 3e-16.

**Enumeration [C].** Charge- and Z2_D-conserving channels open at threshold:

- **Conversions** (2→2, number conserving): hh → XX̄, hh → DD and DD → XX̄. The reverse directions XX̄ → hh and
  XX̄ → DD are **forbidden channels**: closed at rest, open only in the Boltzmann tail. XX̄ → DD is closed by
  only 1.0 GeV.
- **Number-changing 3→2:** 22 channels (11 up to charge conjugation). Examples: XXX̄ → Xh, XX̄h → hh,
  Xhh → Xh, hhh → hh, DDD → hD and XDD → Xh.
- 4→2 processes exist but are not included, since their rate is ∝ n³.

**Two exact threshold zeros [D, C]:**

1. **Pure-h 3→2 at threshold.** |M|² = 7e-39 for hhh → hh, against 3.2e-9 for XXX̄ → Xh. This is the known
   threshold nullification of broken λφ⁴ [R: Brown; Argyres–Kleiss–Papadopoulos; Voloshin, 1992–93; not
   re-read]. The thermal average is nonzero.
2. **hh → XX̄ at threshold.** M_thr = −4g_X²(1 − g_X²/λ_φ) = 0 at the benchmark ratio a = g_X²/λ_φ = 1 [D]. This
   is a benchmark coincidence, not a symmetry, and it makes the conversion p-wave-like:
   ⟨σv⟩(x = 20) = 1.4e-12 GeV⁻² = 1.6e-29 cm³ s⁻¹ [C].

**Thermal averages [C].**

- **3→2:** ⟨σv²⟩ = ⟨(8E1E2E3)⁻¹ ∫dΦ₂|M|²/S_f⟩, by Monte Carlo over relativistic Maxwell–Boltzmann momenta at 12
  temperatures (x = m_X/T = 0.25–300). The result tends to the threshold value at x = 3000 (ratio 0.9957 ± 0.0007).
  The dominant process depends strongly on energy: its thermal average is only 0.046 of the threshold value at
  x = 3, and 0.33 at x = 10.
- **2→2:** Gondolo–Gelmini integrals over σ(√s) from 40-point angular quadrature. They agree with a direct MC of
  ⟨σ v_Møl⟩ to 0.5%.

| Process | x = 3 | x = 10 | x = 30 | threshold |
|---|---|---|---|---|
| ⟨σv²⟩ XXX̄ → Xh [GeV⁻⁵] | 4.4e-19 | 3.1e-18 | 6.3e-18 | 9.5e-18 |
| ⟨σv²⟩ XX̄h → hh | 5.8e-19 | 1.6e-18 | 1.9e-18 | 1.8e-18 |
| ⟨σv²⟩ hhh → hh | 9.6e-19 | 6.8e-19 | 1.6e-19 | 0 |
| ⟨σv⟩ hh → XX̄ [GeV⁻²] | (x = 1: 1.5e-11) | (x = 5: 8.9e-12) | (x = 20: 1.4e-12) | 0 |
| ⟨σv⟩ DD → XX̄ | (x = 1: 3.1e-16) | (x = 5: 3.9e-15) | (x = 20: 6.3e-15) | |
| ⟨σv⟩ hh → DD | (x = 1: 3.6e-15) | (x = 5: 1.6e-14) | (x = 20: 2.4e-14) | |

## 6. Item 2c: Boltzmann equations and Ω_sym(ξ)

**System [D].** The species are X (X + X̄, g = 2, common number chemical potential), h (g = 1) and D (g = 1).
They share one kinetic temperature T_h, use Bose–Einstein n, ρ, P, and their state is (ln T_h, ψ_X, ψ_h, ψ_D)
with ψ_s = (μ_s − m_s)/T_h. With N = ln a:

- dn_s/dN = −3n_s + C_s/H
- dρ_h/dN = −3(ρ_h + P_h). The sector is closed; gravitational exchange with the SM is negligible [Q hidden_heating].
- C_s = Σ_p R_p^eq(T_h) [e^{Σ_A η} − e^{Σ_B η}] (ν_s^B − ν_s^A), with η = μ/T and R^eq the μ = 0 rate per volume.
- H² = (ρ_SM + ρ_h)/(3M_P²), with T_SM ∝ a⁻¹.

The series for n, ρ, P and their analytic derivatives are checked by finite differences to 1e-9. A collisionless run
conserves the comoving hidden entropy and every particle number to 4e-10.

**Initial state [P].** T_h = 4m_X, μ = 0 (thermal reference) and T_SM = T_h/ξ. The end point is x = 2000.

| ξ | Ω_sym h² | Ω / (frozen initial number) | max Γ_3→2/H | max Γ_conv/H |
|---|---|---|---|---|
| 1e-5 | 3.02e-7 | 0.999 | 8.3e-6 | 1.4e-3 |
| 1e-4 | 2.94e-4 | 0.972 | 8.3e-4 | 0.14 |
| 3.16e-4 | 8.96e-3 | 0.937 | 8.3e-3 | 1.3 |
| **7.58e-4** | **0.120** | 0.91 | 0.048 | 7.5 |
| 1e-3 | 0.272 | 0.900 | 0.083 | 13 |
| 3.16e-3 | 7.66 | 0.80 | 0.83 | 128 |
| 1e-2 | 197 | 0.65 | 8.3 | 1.3e3 |
| 0.1 | 1.31e5 | 0.43 | 828 | 1.3e5 |
| 1 | 9.53e7 | 0.32 | 8.1e4 | 1.2e7 |

**What the table shows [C]:**

- Γ_3→2/H ∝ ξ² (to 2% over the whole scan), so **3→2 reaches equilibrium only for ξ > ξ_32 = 3.48e-3**.
- Below ξ_32 the hidden particle number is frozen, and the relic equals the initial number, redistributed by
  conversions.
  - hh → XX̄ is active for ξ ≳ 3e-4. It moves the heavier h into X and lowers Ω by up to ~10%.
  - DD ↔ XX̄ is slow.
- Above ξ_32, cannibalism lowers Ω by at most a factor ≈ 3 at ξ = 1, with 3→2 freezing out at x ≈ 1.4–6.6. In this
  regime T_h/T_SM rises after the hidden sector becomes nonrelativistic (to ≈ 2.9 at x ≈ 26 for ξ = 1): this is
  cannibal heating.
- At ξ_max the composition is Ω_X = 0.083, Ω_h = 0.0066 and Ω_D = 0.030. The total yield is ΣY = 4.30e-12 and the
  mean mass is 101.7 GeV.

**Hot-relic closed form [D].** Y_s = 45ζ(3)g_s ξ³/(2π⁴ g*s), so Ω h² = 3.154e8 ξ³ with every species frozen. It is
2.858e8 ξ³ if h converts fully to X, and 2.141e8 ξ³ if D is absent. The corresponding ξ_max values are 7.25e-4,
7.49e-4 and 8.25e-4, which bracket the ODE values.

- **ξ_max (ODE, Ω_sym = 0.12) = 7.581e-4** with all species populated.
- **ξ_max = 8.364e-4** if D was never populated. D–φ energy exchange needs ξ ≳ 0.015 [Q thermal_potential], so at
  these ξ the D population is a production input.
- **Y_N,max = 4.30e-12 ≈ 0.050 Y_B** for any production history.
- Scaling: ξ_max ∝ (v/TeV)^{−1/3}.

**Forbidden and cannibal prior art.**

- **Forbidden dark matter** (D'Agnolo–Ruderman 1505.07107 [Q abstract]) works because the heavier annihilation
  products are in equilibrium with an external bath. Here h and D are stable and secluded. XX̄ → hh (forbidden)
  only converts quanta and conserves the hidden number. It removes nothing.
- **Cannibal dynamics** (Carlson–Machacek–Hall 1992; Pappadopulo–Ruderman–Trevisani 1602.04219; Farina et al.
  1607.03108 [Q abstracts]) is present (XXX̄ → Xh and others), but only for ξ > 3.5e-3, far above the relic bound.

**Floor inequality** (any symmetric relic only strengthens the conditional H0-floor exclusion) [D].

- Assumption: at late times the hidden population is dilute and nonrelativistic. The bound charge carries at least
  m* per unit, and every other hidden quantum is a free particle of mass ≥ m_X.
- Rigorous weaker form, without that assumption: ω_h ≥ max[(m*|r|/m_b)ω_b, ω_sym].

ω_h ≥ (m*|r|/m_b) ω_b + ω_sym  ⇒  H0 ≥ 100 √(ω_b + ω_h,min + ω_sym) km/s/Mpc.

The maser D_A envelope (Codex 2ee7f82 [Q]) depends only on the total pressureless density, which strictly lowers
D_A. Its necessary condition becomes

m_X|r| + (m_X/m*)(m_b/ω_b) ω_sym < 63.93 GeV (primary) / 74.64 GeV (stress).

- The coefficient is 79.28 GeV per unit ω_sym, with the envelope's conservative ω_b = 0.02125 and m_b = 0.94738.
- For ω_sym = 0.12: |r| < 0.544 (primary) and 0.651 (stress).
- The r = 1 floor rises from 113.25 to 118.43 km/s/Mpc.
- Separately, ω_asym ≥ 1.338 r. The retrospective DM datum would require |r| ≤ 0.0897 if ω_h ≤ 0.12.

## 7. Item 3: is a Z2-restoring thermal history compatible with the relic bound?

**Minimum visible temperature [C].** With ξ conserved between T_c and the nonrelativistic epoch, as argued in §2,
restoration at T_h = T_c happens at T_SM = T_c/ξ.

| T_c/v | ξ_max all species | ξ_max no D | hot-relic all |
|---|---|---|---|
| 1.5476 (leading) | **2.04e6 GeV** | 1.85e6 | 2.14e6 |
| 1.5571–1.5641 (resummed [Q]) | **2.05–2.06e6 GeV** | 1.86–1.87e6 | 2.15–2.16e6 |

**Robustness to the number chemical potential [D, C].** At fixed hidden number density,

Π_φ ∝ n^{2/3} Li₂(z)/Li₃(z)^{2/3}

(massless bath, fugacity z ≤ 1). This is monotonically increasing in z (checked on 9 points), so μ_N = 0 restores
with the fewest quanta. Restoration therefore needs n ≥ n_th(T_c) = 4ζ(3)T_c³/π² = 1.81–1.86 v³ for any μ_N ≤ 0.
With the relic bound Y_N ≤ 4.30e-12, this gives **T_SM ≥ 2.08–2.10e6 GeV (v/TeV)^{4/3}** [C]. The bound does not
depend on how the bath was produced.

A coherent X population could restore with far fewer quanta. Using the cold charged-condensate density
n_c ≈ 0.0316 v³ as a **PLACEHOLDER**, the bound for that route is T_SM ≳ 5.4e5 GeV [E].

**Window [C, Q].**

| Condition | Value of ξ | Status |
|---|---|---|
| Kinetic equilibrium of φ at T_c [Q thermal_potential] | ξ ≳ 6.6e-5 (v/TeV)^{1/2} | satisfied: a window up to ξ_max = 7.6e-4 **exists** |
| 2→4 chemical equilibrium [Q] | ξ ≳ 0.068 (band 0.0084–0.27) | **incompatible** with the bound |
| Broken-phase 3→2 equilibrium [C] | ξ > 3.5e-3 | **incompatible** with the bound |

**Answer:** a restoring thermal history is compatible with the relic bound only with small ξ ≤ 7.6e-4 and a large
visible temperature, T_SM ≥ 2.0e6 GeV at restoration. Such a bath is kinetically thermal but **never chemically
equilibrated**. Its number chemical potential, and hence the zero-μ T_c premise, are production inputs [P]. The
needed outside-action hidden reheating is, for example, a branching B_h = (g_h/g*)ξ⁴ ≲ 1.2e-14 at T_max ≥ 2e6 GeV
(see hidden_heating §5 for the formula).

The relic bound can be relaxed by outside-action entropy injection into the SM after restoration, with dilution
factor Δ: ξ_max → ξ_max Δ^{1/3} and T_SM,min → T_SM,min Δ^{−1/3}.

## 8. Item 4: hidden charge in a thermal hidden plasma

**Capacity [C, D].** In the restored phase the X thermal mass is M = √Π_X T = 0.03416 T. A complex scalar gas
holds at most n_Q,max = n_Q(μ = M) = 0.011386 T³. The leading form is M T²/3; the exact BE integral agrees to 3e-5.
Excess charge Bose-condenses. The charge needs n_Q/T_h³ = r Y_B (2π²g*s/45)/ξ³ = 4.06e-9 r/ξ³, so

**ξ_BEC = 7.09e-3 r^{1/3}** [C].

| r | ξ | μ/T at T_c | Condensed fraction at T_c |
|---|---|---|---|
| 1 | 7.58e-4 | = M/T (0.0342) | 0.9988 |
| 0.0897 | 7.58e-4 | M/T | 0.986 |
| 0.01 | 7.58e-4 | M/T | 0.878 |
| 1e-3 | 7.58e-4 | 0.0282 | 0 |
| 1 | 0.01 | 0.0124 | 0 |
| 1 | 1 | 1.24e-8 | 0 |

**Asymmetric versus symmetric component.**

- The symmetric pairs are not annihilated away. XX̄ → hh conserves the hidden number, and 3→2 is inactive for
  ξ < 3.5e-3. So the relic is asymmetric **plus** symmetric, unlike standard asymmetric dark matter.
- At ξ_max the symmetric yield is 4.3e-12. The charge yield is 8.7e-11 r, so the effective asymmetry is
  η_eff = r Y_B/(r Y_B + Y_sym).

**Charge outlives thermal restoration [C].** At the end of a relic-compatible restoration (ξ = ξ_max, T_SM = T_c/ξ),
n_Q/n_c = 1.1e3 r, using the cold-sector n_c as a **PLACEHOLDER**. For r ≳ 1e-3 the charge alone therefore keeps
φ = 0 after T_h drops below T_c, until T_SM = T_n(r) = 198 r^{−1/3} TeV.

- Equivalently, the charge holds whenever ξ < T_c/T_n(r) = 7.8e-3 r^{1/3}.
- The exit is then the charge_phase route [Q]: a spinodal at n ≈ n_c, with p(+) = p(−) by the exact Z2 and a
  sign-correlation length of 1e2–1e3 ℓ.

**Solitogenesis [E].** At the spinodal exit the charge per Kibble–Zurek correlation volume is n_c ξ_KZ³ ≈ 1.9e10
(r = 1) or 8.5e8 (r = 0.01), using charge_phase ξ_KZ = 850 and 300 ℓ [Q]. Both are far above the thin-wall
binding threshold Q ≈ 1.3e4 [Q]. Charge-restored regions can therefore leave bound Q-balls. With an exact Z2, the
vacuum between them takes random signs, and both signs percolate (p = 1/2 > 0.31). This is the FGGK "exact
degeneracy" case.

**Prior art: how each work handled the Z2 walls.** (a) assumed a bias; (b) used a potential without a degenerate Z2;
(c) analysed the walls.

| Work | Read level | Class | What was done |
|---|---|---|---|
| Frieman, Gelmini, Gleiser, Kolb, PRL 60, 2101 (1988), "Solitogenesis: primordial origin of nontopological solitons"; SLAC-PUB-4538 / FERMILAB-PUB-88/13-A | **full** (Fermilab preprint scan, OCR text) | **(b) + (a)** | The FLS potential has an h\|φ\|²(σ − σ0)² coupling and a σ³ term (a = λ2/λ1 = 0.15); the reflection symmetry is explicitly broken and the vacua are non-degenerate. They state that in the limit h = λ2 = 0 the symmetry is exact and stable walls form and dominate unless σ0 ≲ 10 MeV. Solitogenesis needs the bias, whose vacuum pressure shrinks or pinches false-vacuum regions. The percolating (0.31 < p+ < 0.69) case is treated only with the bias. |
| Griest, Kolb, Massarotti, PRD 40, 3529 (1989); FERMILAB-PUB-89-110-A | introduction + model section (OCR of scan, equations partly illegible) + abstract | **(b)** | Same FGGK potential with λ2/λ1 = 0.15. The soliton wall (surface) energy is ignored. Statistical charge fluctuations replace the cosmic asymmetry. No exact-Z2 wall analysis. |
| Frieman, Olinto, Gleiser, Alcock, PRD 40, 3241 (1989); FERMILAB-PUB-89-128-A | §§I–III + figure captions (OCR) + abstract | **(b)** | Generic NTS evolution (evaporation, accretion, freeze-out). Formation is discussed with schematic potentials that have a unique true vacuum, plus "pinch-off in the asymmetric vacuum model". Exact-Z2 walls are not analysed. |
| Macpherson, Campbell, PLB 347, 205 (1995), hep-ph/9408387, "Biased discrete symmetry breaking and Fermi balls" | abstract + introduction + formation section | **(a) + (c)** | An "approximate" discrete symmetry: bias term A(φ) with vacuum splitting Λ. They note that an exact Z2 gives wall domination (Zel'dovich–Kobzarev–Okun). They analyse the biased percolating network, following Gelmini–Gleiser–Kolb, and its fragmentation into Fermi balls. |
| Krylov, Levin, Rubakov, PRD 87, 083528 (2013), arXiv:1301.0354 | abstract + §2 including footnote 1 | **(a) assumed, not analysed** | Tree potential λ(φ² − v²)² + h²φ²\|χ\|², the same structure as Candidate 1. Footnote 1: the φ → −φ symmetry makes the model non-viable through the domain-wall problem; they assume explicit breaking or a charged φ, "irrelevant for what follows". |
| Pontón, Bai, Jain, JHEP 09 (2019) 011, arXiv:1906.10739 | abstract + model section and section list (searched) | **(b)** | The confining field is the SM Higgs doublet: a connected, gauged vacuum manifold with no discrete-vacuum walls. Strong first-order EW transition via a singlet. |
| Bai, Long, Lu, PRD 99, 055047 (2019), arXiv:1810.04360, "Dark quark nuggets" | abstract + introduction (searched) | **(b)** | Dark-QCD first-order transition with no real Z2 scalar. Walls are mentioned only for other authors' axion quark nuggets. |
| Hong, Jung, Xie, PRD 102, 075028 (2020), arXiv:2008.04430, "Fermi-ball dark matter from a first-order phase transition" | abstract + §II + toy-model section | **(b)** | Toy model U = ½(μ² + cT²)φ² + (μ3/3)φ³ + (λ/4)φ⁴ with a tree-level φ³, so no Z2. A Z2 appears only as an R-parity-like stabiliser of the fermion. |

**Summary [Q].** No work analyses a persistent wall network of an **exact** Z2 in an FLS-type formation history.
Six of the eight remove the degeneracy in the potential, either by construction or by a cubic term. Macpherson–Campbell
and Krylov–Levin–Rubakov keep the Z2 structure but break it explicitly. FGGK and Macpherson–Campbell state that the
exact limit gives the domain-wall problem. BRIEF rule 5 forbids every such bias for Candidate 1, so no prior-art
solitogenesis mechanism disposes of its walls. The prior art supplies no persistence calculation either.

## 9. Item 5: charge supplied at T_Q into a cold hidden sector

This uses the homogeneous critical density n_c ≈ 0.0316 v³ as a **PLACEHOLDER**. charge_phase confirms it as the
T_h = 0 lower spinodal; the forced-restoration threshold is n_s1 = 0.0428 v³ [Q].

No charge-driven restoration occurs if and only if

n_Q(T_Q) = r Y_B (2π²g*s/45) T_Q³ < n_c ⇔ **T_Q < T_n(r) = 198.2 r^{−1/3} (v/TeV) TeV** [D].

| r | T_Q,max (n_c) | T_Q,max (n_s1) |
|---|---|---|
| 1 (benchmark) | **198.2 TeV** | 219.2 TeV |
| 0.6393 (maser primary: m_X\|r\| < 63.93 GeV) | **230.1 TeV** | 254.4 TeV |
| 0.7464 (maser stress) | 218.5 TeV | 241.6 TeV |
| 0.5442 (primary, with ω_sym = 0.12) | 242.8 TeV | 268.5 TeV |
| 0.0897 (retrospective DM, ω_h ≤ 0.12) | 442.8 TeV | 489.7 TeV |
| 0.01 | 920.0 TeV | 1017 TeV |

Because T_n decreases with r, **T_Q < 230 TeV (v/TeV)^{4/3}** guarantees no charge-driven restoration for every
maser-allowed r. Here r_max = 63.93 GeV/(g_X v).

These conditions apply to charge supplied homogeneously into a cold hidden sector (T_h ≪ T_c) whose φ sign is
already selected. Codex's production preparation corresponds to T_vis = 175.6 TeV < 198 TeV, which is charge_phase's
H2 case.

## 10. Conditional history constraints (no prediction, no exclusion, no validation)

- **R1 (charge).**
  - The hidden net charge is an outside-action datum: supplied at a free T_Q, after or within the last ≲ 14–20 e-folds.
  - Supplied homogeneously into a cold sector at T_Q < 198 r^{−1/3} (v/TeV) TeV, it cannot restore Z2 by itself.
  - Above that, it can, and the exit is the random-sign spinodal (charge_phase H1).
- **R2 (symmetric relic).**
  - Any hidden bath that was ever relativistic leaves X, h and D relics, which neither forbidden annihilation nor
    cannibalism remove below ξ = 3.5e-3.
  - Ω_sym h² < 0.12 (retrospective) requires Y_N ≤ 4.3e-12, i.e. ξ ≤ 7.6e-4 (8.4e-4 without D) for a thermal
    number.
  - Any symmetric relic only tightens the H0-floor condition: m_X|r| + 79.28 GeV ω_sym < 63.93 GeV.
- **R3 (restoring history).**
  - A thermal Z2 restoration compatible with R2 needs T_SM ≥ 2.0–2.1e6 GeV (v/TeV)^{4/3} at restoration,
    6.6e-5 ≲ ξ ≤ 7.6e-4, and a number chemical potential fixed by production.
  - If r ≳ 1e-3, the charge then dominates the exit (§8).
  - Formation of walls in such a history is the charge_phase / thermal_potential statement, p(+) = p(−) with both
    signs percolating. **Persistence is not addressed here.**

## 11. Free inputs on which these results depend

- **[P]** Hidden production:
  - ξ, defined at the relativistic epoch;
  - the produced number, i.e. μ_N, set here to 0 as a thermal reference;
  - whether D was populated;
  - T_h,max.
- **[P]** Y_Q = r Y_B, T_Q, and the spatial form of the supplied charge.
- **[P]** H_inf and the reheating history (for N_max); later entropy dilution Δ.
- The scale v: ξ_max ∝ v^{−1/3}, T_SM,min ∝ v^{4/3}, T_n ∝ v r^{−1/3}. Also the benchmark couplings, in particular
  the coincidence a = g_X²/λ_φ = 1, which zeroes the hh → XX̄ threshold amplitude.
- **Retrospective data:** ω_DM = 0.120 (Planck 2018) and the Codex maser product bounds. The ω_b inputs are BBN
  calibrations.
- **Placeholders:** n_c = 0.0316 v³ (charge_phase); ξ_kin, ξ_chem and T_c from thermal_potential.

## 12. What failed or was corrected

1. **First Boltzmann formulation failed its own conservation test.** It used η = μ/T variables, and a collisionless
   run drifted by 5% in entropy and 6–9% in numbers over 8 e-folds, because η ≈ x ≈ 1e5 lost precision. It also
   overflowed: the log of the MB density was clipped once e^{−x} underflowed at x > 745. It was replaced by
   ψ = (μ − m)/T variables and log-space densities. Conservation is now 4e-10, and every result here comes from the
   corrected code.
2. **Threshold-only ⟨σv²⟩ was wrong by up to 20×.** It overestimated the dominant 3→2 rate (by 20× at x = 3 and
   3× at x = 10) and missed the nonzero thermal pure-h rate. It was replaced by Monte Carlo thermal averages.
3. **My first hand estimate of σv(hh ↔ XX̄)** (~1e-9 GeV⁻², WIMP-like) was wrong by about three orders of
   magnitude. At a = 1 the threshold amplitude vanishes exactly (derived in §5).
4. **Finite-difference tensor check with a 0.01 GeV step failed** (3e-4 absolute in the h⁴ entry, from roundoff).
   With a 1 GeV step, exact for a quartic polynomial, it agrees to 1e-12.
5. **Literature access.** The guessed Fermilab URL for the GKM preprint returned 404; it was obtained through
   INSPIRE's file link. The SLAC copy of FGGK failed TLS verification and was not needed. GKM and FOGA are OCR
   scans with partly illegible equations, which is why their read levels are partial.

## 13. Not done

- 4→2 processes and Bose-enhanced collision kernels.
- Thermal masses and φ_min(T) inside the rates. Their effect at T_h ≤ 4m_X is ≈ 3%.
- D kinetic decoupling. A common T_h is assumed; this does not affect Ω at ξ_max, where the D number is frozen.
- g*(T) variation below 200 GeV. It enters only after freeze-out for ξ ≲ 1e-2.
- Coupling running in the rates.
- The relativistic 2↔3 rates between T_c and 4m_X are not computed. I argued they are smaller: Γ/H ∝ ξ²/T, and the
  cubic vertex is ∝ φ_min(T) → 0 at T_c.
- No production mechanism, hence no prediction of ξ, μ_N or r.
- No charged-condensate dynamics at T_h ≈ T_c, beyond the capacity and threshold estimates.
- **No wall formation or persistence calculation, network evolution or observable.** Those belong to other modules.
- No quantum-gravity charge violation.

## 14. Literature and sources (no paper text copied; downloads kept in the session scratchpad)

| Reference | Read level | Used for |
|---|---|---|
| Frieman, Gelmini, Gleiser, Kolb, PRL 60, 2101 (1988), doi:10.1103/PhysRevLett.60.2101 (preprint FERMILAB-PUB-88/13-A) | full (OCR text of preprint) | prior-art table: potential, bias, wall statements |
| Griest, Kolb, Massarotti, PRD 40, 3529 (1989), FERMILAB-PUB-89-110-A | abstract + introduction/model (OCR, partial) | prior-art table |
| Frieman, Olinto, Gleiser, Alcock, PRD 40, 3241 (1989), FERMILAB-PUB-89-128-A | abstract + §§I–III, captions (OCR) | prior-art table |
| Macpherson, Campbell, PLB 347, 205 (1995), arXiv:hep-ph/9408387 | abstract + introduction + formation section | prior-art table |
| Krylov, Levin, Rubakov, PRD 87, 083528 (2013), arXiv:1301.0354 | abstract + §2 (footnote 1) | prior-art table |
| Pontón, Bai, Jain, JHEP 09 (2019) 011, arXiv:1906.10739 | abstract + model section (searched) | prior-art table |
| Bai, Long, Lu, PRD 99, 055047 (2019), arXiv:1810.04360 | abstract + introduction (searched) | prior-art table |
| Hong, Jung, Xie, PRD 102, 075028 (2020), arXiv:2008.04430 | abstract + §II + toy-model section | prior-art table |
| D'Agnolo, Ruderman, PRL 115, 061301 (2015), arXiv:1505.07107 ("Forbidden Dark Matter" on arXiv) | abstract | forbidden-channel logic (§6) |
| Carlson, Machacek, Hall, ApJ 398, 43 (1992), doi:10.1086/171833 | abstract (INSPIRE) | cannibal logic |
| Pappadopulo, Ruderman, Trevisani, PRD 94, 035005 (2016), arXiv:1602.04219 ("Cannibal Dark Matter" on arXiv) | abstract | cannibal logic |
| Farina, Pappadopulo, Ruderman, Trevisani, JHEP 12 (2016) 039, arXiv:1607.03108 | abstract | cannibal phases |
| Planck 2018 VI, arXiv:1807.06209 | abstract | ω_c = 0.120 ± 0.001 (retrospective datum) |
| Gondolo, Gelmini, NPB 360, 145 (1991) | not read; formula recalled [R] | thermal average; verified by direct MC to 0.5% |
| Birrell & Davies (heat-kernel a_2, Riem²/180) | recalled [R] | h → gg loop estimate |
| Threshold nullification in broken λφ⁴ (Brown; Argyres–Kleiss–Papadopoulos; Voloshin, 1992–93) | recalled [R], not re-read | interpretation of the hhh → hh zero |
| LHC Higgs XS WG SM width near 141 GeV ≈ 8.6 MeV | recalled [R] | ξ-trace decay estimate (O(1) uncertainty irrelevant) |
| Gelmini, Gleiser, Kolb, PRD 39, 1558 (1989); Zel'dovich, Kobzarev, Okun 1974 | cited via Macpherson–Campbell / FGGK only, not read | context |
| Repository: BRIEF.md; GRUT_CANDIDATE_1.md §§1–7; sibling NOTES (charge_phase, thermal_potential, hidden_heating); Codex 135d535 README, DERIVATION, FEASIBILITY, CLAIM_EVIDENCE (git show); Codex 2ee7f82 REPORT.md (git show, sections) | read | inputs, thresholds, floor and maser bounds |

## 15. Files

| File | Contents |
|---|---|
| `code/common.py` | constants, potential and vertex tensors, BE thermodynamics (η and ψ forms), MB densities |
| `code/amplitudes.py` | Berends–Giele tree amplitudes, kinematics, 2→2 σ, MC thermal ⟨σv²⟩, process enumeration |
| `code/r1_charge.py` | item 1 → `results/charge_conservation.json` |
| `code/r2_spectrum.py` | spectrum, stability, h lifetime → `results/spectrum_stability.json` |
| `code/r3_rates.py` | all processes, checks, thermal rates → `results/rates.json` |
| `code/r4_boltzmann.py` | Boltzmann ODE, Ω(ξ), ξ_max, restoration table → `results/boltzmann.json` |
| `code/r5_charge_history.py` | items 3–5 support, floor inequality, prior art → `results/charge_history.json` |
| `code/run_all.py` | everything, plus `results/headline.json` (24 assertions, code hashes) |
