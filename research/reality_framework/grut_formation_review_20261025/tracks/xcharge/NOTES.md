# Track `xcharge`: independent check of the explicit U(1)_X-breaking test in 8986d34

Claude review, 25 October 2026. Read-only on Codex's branches. No commits.

The packet was read with `git archive 8986d34` into the scratchpad. This covers
`research/reality_framework/grut_cstar_audit_20261024/` plus the parent inputs it imports at the same commit:
`grut_core_conversion_20261023/code/scattering.py` and the `grut_portal_debris_window_20261023` profile files.
I read CLAUDE_HANDOFF.md, REPORT.md, MATHEMATICS.md, GLOBAL_CHARGE.md, CSTAR_HYPOTHESIS.md, LITERATURE_AUDIT.md,
REPRODUCE.md, `code/charge_leakage.py`, `code/verify.py`, `code/report.py` and CHARGE_LEAKAGE_RESULTS.json.

**Independence.** Every number marked [C] comes from my own code in `code/`. I solved the profile myself and did
not use Codex's CSV, and I derived the linearization and normalization myself. I solved the forced problem with
two numerical methods that differ from Codex's `solve_bvp` collocation, and I added a full nonlinear time-domain
simulation that Codex does not have. I ran no Codex code, so there are no reproduction numbers. Codex's numbers are
quoted [Q] from its JSON.

Labels: [D] derived, [C] computed here, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

Scope. This covers one model with supplied fields, one frozen Q = 200000 ground object and a hypothetical operator
with an unknown coefficient. Nothing here is a quantum-gravity prediction, a lifetime certificate or evidence that
the fields exist.

## 0. Verdict summary

| # | Claim (8986d34) | Codex [Q] | Mine | Verdict |
|---|---|---|---|---|
| 1 | Canonical normalization X = vψ/√2, Q = −2∫Im(X*Ẋ), source in the anti-X channel, 3 channels at ((d+1),(d−1),d)Ω, all open for d = 5, 6 | as stated | Re-derived [D]. All K² > 0 [C] | CONFIRMED |
| 2 | Rate prefactor Ṅ_i = 4π m_X/(g²ħ)·\|ε_C\|²·K\|A\|², ε_C = d·2^{1−d/2}c v^{d−2}/(m_X²M^{d−4}); √2 real-scalar normalization; 4π angular flux | as stated | My convention gives 8π m_X/(g²ħ)·\|ε_m\|²·K\|A\|² with ε_m = ε_C/√2, which is identical [D]. Golden-rule free limit agrees [D] | CONFIRMED |
| 3 | Core ledger −Q̇_core = (d+1)Ṅ₊ + (d−1)Ṅ₋ + dṄ_s; total violation d·ΣṄ | as stated | Re-derived [D]. Equal to the time-domain drift [C] | CONFIRMED |
| 4 | First-order torque is a bounded oscillation with zero cycle average; secular loss is O(\|c\|²) | argued | Full nonlinear simulation, ±η runs: the odd part is a pure oscillation with amplitude 2ηI_d/Ω, matching to 4e-5 to 3e-3. The even drift scales as η² (×4.00, ×4.01) [C] | CONFIRMED, for the single-frequency ground state (caveat in §2) |
| 5 | Initial fractional core-charge loss, X⁵, unit c, M = 2.435e18 GeV (coherent single-quantum channel) | 1.5360e-12 /s | 1.53600e-12 /s (two solvers agree to 7e-12). Time-domain drift matches to 0.03–0.2% [C] | CONFIRMED (that channel) |
| 6 | Same, X⁶ | 3.9196e-50 /s | 3.91960e-50 /s; time domain 0.2–0.5% [C] | CONFIRMED (that channel) |
| 7 | \|c₅\| bound for < 1% in 13.8 Gyr | 1.23e-4 (CSTAR_HYPOTHESIS.md, user summary); 0.000122 (REPORT.md); JSON 1.22268e-4 | 1.22268e-4 from the coherent channel alone [C] | CORRECTED (rounding): 1.22e-4, not 1.23e-4 |
| 8 | Unit-coefficient scale for the 1% target | 1.99e22 / 8.8e10 GeV | 1.9915e22 / 8.80e10 GeV [C] | CONFIRMED (coherent channel) |
| 9 | Peak-to-peak first-order excursion fraction | 2.673e-16 / 1.583e-34 | 2.6730e-16 / 1.5834e-34 [C] | CONFIRMED |
| 10 | The leading-response forced radiation is the leading O(\|c\|²) charge drain ("leading-response approximation") | implicit | **No.** Two-quantum anti-X pair emission, from the C(d,2)λX₀^{d−2}δX² term, is kinematically open. It is not suppressed by the surface form factor and exceeds the coherent rate by about 45× (d = 5) and 2.7e3× (d = 6) [E, free-Born + LDA] | CORRECTED |
| 11 | X⁵ unit c at M_red fails the 1% target; X⁶ passes | fail / pass | Still fail / pass with pairs included. Total ≈ 7.1e-11 /s (d = 5), 1.04e-46 /s (d = 6) [E] | CONFIRMED (outcome) |
| 12 | \|c₅\| ≲ 1.2e-4 suffices | 1.22e-4 | ≈ 1.8e-5 with pairs (range 1.76–1.84e-5 over estimators; O(1) medium uncertainty → roughly 1e-5 to 3e-5) [E] | CORRECTED |
| 13 | Rate ratio d = 5 : d = 6 ∼ (M/f)² | not stated | Coherent: 3.92e37 = 19.4 (M/\|X₀(0)\|)². Pairs: 6.7e35 = 0.33 (M/\|X₀(0)\|)², with \|X₀(0)\| = 1.712 GeV [C/E] | scaling CONFIRMED, prefactor channel-dependent |
| 14 | Survival through formation imposes a stronger bound (task question) | "hot core untested" | No, within the supplied history. The homogeneous 20 MeV condensate has Γ/H ≈ 9e-38 (d = 5). The ground family gives pair rates of 1.5–9.1e-11 /s for Q = 5.6e4–5.3e5. Beating today's 13.8-Gyr exposure needs a ≳ 4e11–4e17 rate boost [E/A] | NEW (not stronger) |
| 15 | Z₆ discrete gauge removes X⁵; X⁶ at M_red passes | candidate | X⁶ total Γt = 4.5e-29, a margin of ~2e26 below 0.01. The protection logic matches Krauss–Wilczek / Banks–Seiberg [A] | CONFIRMED as candidate (UV not built) |
| 16 | C* side identities: S_dS = 8π²M²/H_Λ² ≈ 3.31e122; closed dS a = cosh(Ht)/H satisfies the constraint with p = −ρ; a ∝ t² has w = −2/3 and event horizon d_E = t | as stated | 3.306e122 [C]; identities re-derived [D] | CONFIRMED (side check only) |
| 17 | Remaining C* clauses (Hilbert-space counting, observer dictionary, q-equilibrium) | — | outside this track | NOT CHECKED |

## 1. Model and inputs [Q]

L = |∂X|² + ½(∂φ)² − g²φ²|X|² − λ_X|X|⁴ − (λ_φ/4)(φ² − v²)². Here g = 0.1, v = 0.924089637375362 GeV and
m_X = gv = 92.41 MeV. A = λ_X/g² = 0.1 and B = λ_φ/g² = 3, so m_s = √6 m_X. The object is the Q = 200000 ground
state. M = 2.435e18 GeV and t_ref = 13.8 Gyr (Julian years, 4.3549e17 s). The model comes from the
portal-debris MATHEMATICS.md and `scattering.py`; Q, M and t_ref come from `charge_leakage.py`.

**My profile [C]** (`code/profile.py`) uses collocation at fixed Ω with Robin outer conditions, then a Brent root
in Ω for g²Q = 4πΩ∫ρ²F². It gives Ω = 0.8611951308219, F₀ = 2.6199231342, S₀ = 0.0060123211 and
E/(Qm_X) = 0.9233161354. These agree with Codex's refined table to about 2e-13.

## 2. Derivation [D]: how X^d drains a rotating core

Write X₀ = (f/√2)e^{−iωt} with f = vF, and δV = λX^d + λ*X*^d with λ = c/M^{d−4}.

*Ward identity:* ∂_μ j^μ = −2d·Im(λX^d). On X₀, X₀^d ∝ e^{−idωt}, so ∫ j̇⁰ = 2d|λ|I sin(dωt) with
I = ∫|X₀|^d d³x. The first-order charge change is −(2|λ|I/ω)cos(dωt). This is bounded, with peak-to-peak
4|λ|I/ω and zero cycle average.

*Second order:* the response δX ∝ λ* beats against λX₀^{d−1} to give a DC part, and the surviving DC part is the
outgoing radiative flux. Linearize with δX = e^{−iωt}[a e^{−iνt} + b* e^{iνt}], δφ = c e^{−iνt} + c.c.,
ν = dω and fields in units of v. The source enters only the b equation: −ε_m F^{d−1}, with
ε_m = d·2^{(1−d)/2}λv^{d−2}/m_X². The lab frequencies are (d+1)ω (X), (d−1)ω (anti-X) and dω (s); ν = dω is
the rotating-frame frequency. For U = ρ(a, b, c), the potential matrix is real and symmetric, with off-diagonal
entries AF² and √2SF. That makes Σ Im(U*U′) the conserved number flux with equal weights.

*Physical flux:* a channel with amplitude A e^{ikr}/r carries 8πk|A|² quanta per unit time, for the complex field
and for the real field normalized as φ = c e^{−iνt} + c.c. (direct T^{0r} calculation). Converting to
dimensionless units gives Ṅ_i = (8π m_X/(g²ħ))|ε_m|²K_i|Â_i|². Codex's (4π, ε_C = √2ε_m) is identical. In the free
limit this equals Fermi's golden rule Γ = k|J̃(k)|²/(2π) for J = dλX₀^{d−1}.

*Ledger:* adiabatically dE_core = ω dQ_core, so −Q̇_core = Σ(Ω_i/Ω)Ṅ_i = (d+1)Ṅ₊ + (d−1)Ṅ₋ + dṄ_s. Including
the emitted charge, the total violation is d·ΣṄ.

*Caveat on "averages to zero":* this holds for a single-frequency (circular) condensate. A mixed-frequency
condensate can have a DC first-order torque, as in Affleck–Dine: with components at ±ω, X^d has a DC piece
a^{d/2}b^{d/2} for even d. That is relevant only to non-maximally-charged formation states.

## 3. Coherent single-quantum rate [C] (`code/forced.py`)

Two solvers were used: (A) particular-plus-homogeneous DOP853 shooting with outgoing matching; (B) the Wronskian
formula A = W⁻¹∫Φᵀ S dρ, which needs only regular solutions. They agree to 7e-12 (d = 5) and 4e-11 (d = 6). The
source-work identity holds to 7e-12. Results do not depend on R = 30, 40 or 55.

| d | K (X, anti-X, s) | Unit flux K\|Â\|² (X, anti-X, s) | Core-loss share | Γ_Q /s, mine | Codex |
|---|---|---|---|---|---|
| 5 | 5.069, 3.296, 3.541 | 1.398e-5, 0.080701, 0.012799 | 0.02%, 83.4%, 16.5% | 1.53600e-12 | 1.53600e-12 |
| 6 | 5.945, 4.188, 4.550 | 1.109e-6, 0.017483, 0.0012978 | 0.01%, 91.8%, 8.2% | 3.91960e-50 | 3.91959e-50 |

With V = 0, the free control fluxes are 0.143139 and 0.027155, the same as Codex's analytic Green integral. Medium
distortion lowers the coherent rate by about 0.68×.

## 4. Full nonlinear time-domain check [C] (`code/timedomain.py`)

The method evolves the complete classical equations with δV switched on smoothly. Settings: radial 4th-order
finite differences, h = 0.04, RK4 dt = 0.01, polynomial sponge on ρ > 50, core charge inside ρ < 25, T = 420. Runs
at +η, −η and 0 cancel the background breathing. The odd part gives the O(η) term and the even part the O(η²) term.
The fit window is τ ≥ 160.

| d | η | Even drift d(g²Q)/dτ | Predicted −8πΣ(W/Ω)K\|Â\|²ε² | Odd oscillation amplitude / 2ηI_d/Ω |
|---|---|---|---|---|
| 5 | 2.5e-4 | −9.4888e-7 | −9.4956e-7 | 1.76078 / 1.76084 |
| 5 | 5e-4 | −3.7969e-6 (h = 0.03: −3.7971e-6) | −3.7982e-6 | 3.52115 / 3.52168 |
| 5 | 1e-3 | −1.5222e-5 | −1.5193e-5 | 7.0388 / 7.0434 |
| 6 | 5e-4 | −6.745e-7 | −6.730e-7 | 5.4939 / 5.4972 |
| 6 | 1e-3 | −2.704e-6 | −2.692e-6 | 10.965 / 10.994 |
| 6 | 2e-3 | −2.49e-4 (23× prediction) | −1.08e-5 | — |

At larger drive the classical system goes parametrically unstable. For d = 6 at η = 2e-3 the drift is 23× the
prediction. For d = 5 at η = 4e-3, an earlier single-sign run gave a drift 1.5e4× the prediction. For d = 5 at
η = 1.6e-2 and d = 6 at η ≥ 8e-3 the runs blow up (`results/timedomain_large_eta_*.log`). This is the classical,
stimulated counterpart of the pair channel in §5. At physical η ≈ 1e-17 the instability is absent, but spontaneous
pair emission remains.

## 5. Final states and the missing pair channel [D/E] (`code/pairs.py`)

Removing m condensate quanta yields final states with charge m − d and energy mω (ω = 0.8612 m_X). For the direct
operator, with anti-X only, n anti-X are open for n ≤ dω/(ω + m_X): n = 1, 2 for both d = 5 and d = 6, and n = 3 is
closed. Through background vertices, more channels open, and they enter at higher order:

- X at (d+1)ω and s at dω are open.
- For d = 5: anti-X + s at 4ω = 3.4448 is closed against a threshold of 3.4495, and ss is closed.
- For d = 6: anti-X + s and ss are open.

The pair term C(d,2)λX₀^{d−2}δX² creates anti-X pairs at total energy (d−2)ω (2.58 and 3.44 m_X). Pairs can be
emitted back to back, so they draw on G̃(q) at small q. The single quantum must carry K = 3.3–4.2 from a smooth
edge, which costs about 3e-8 in |form factor|².

Free-Born rate, with vacuum dispersion and the finite source:
Γ = 2C²|λ|²v^{2(d−2)}2^{−(d−2)}m_X^{−3}J, with core loss (d−2) per pair.

| d | Pair Γ_Q /s, free-Born | LDA vacuum | LDA with local m_eff² = S² + 2AF² | Pair / coherent |
|---|---|---|---|---|
| 5 | 6.93e-11 | 7.23e-11 | 6.65e-11 | 45 |
| 6 | 1.040e-46 | 1.049e-46 | 1.013e-46 | 2.7e3 |

The Bogoliubov/medium corrections are O(1) and not computed. The single-quantum analog was 0.68×, and the
splitting-vertex interference is about AF²/2 ≈ 0.34 in amplitude. The conclusion survives anything short of a
>40× suppression.

Hertzberg's oscillon analysis (PRD 82, 045022 (2010), arXiv:1003.3459, abstract read) is the precedent:
classical emission is exponentially suppressed, while quantum emission is power-law and often dominant.

## 6. Bounds and scaling

| d | Total Γ_Q /s [E] | Γt | \|c\| for 1% | Unit-c scale | \|c\| for the Poulin-type target fΓt ≈ 0.087 |
|---|---|---|---|---|---|
| 5 | 7.09e-11 | 3.09e7 | 1.80e-5 | 1.35e23 GeV | 5.3e-5 |
| 6 | 1.04e-46 | 4.5e-29 | 1.5e13 | 6.3e11 GeV | — |

For the last column, Poulin–Serpico–Lesgourgues (JCAP 08 (2016) 036, abstract read) give fΓ < 6.3e-3 /Gyr for DM
decaying into dark radiation. Treating the emitted anti-X as dark radiation is approximate, because they become
nonrelativistic later.

Γ ∝ |c|²/M^{2(d−4)} holds exactly at this order. Amplitude scaling: the coherent rate goes as
|X₀|^{2(d−1)}·|FF(K)|² and the pair rate as |X₀|^{2(d−2)}·volume. Per charge, the pair rate saturates for large
thin-wall Q.

The ground-family scan (`results/scan_ground_family.json`, Ω = 0.825–0.95, Q = 5.6e4–5.3e5, F₀ = 2.15–2.62) shows:

- The coherent d = 5 rate swings by about 400× across the family, from 3.5e-15 to 1.5e-12. These are diffraction
  nodes of the edge form factor, and Q = 2e5 sits near a maximum.
- The pair rate is smooth, from 1.5e-11 to 9.1e-11, and dominates everywhere.

## 7. Formation era [E/A]

- **Supplied homogeneous condensate (T = 20 MeV, n = 1.82e-13 GeV³, |X| = 9.9e-7 GeV).** Coherent emission is
  forbidden, because the source has k = 0. The pair rate per charge is 2.4e-35 /s for d = 5, so Γ/H = 8.8e-38;
  for d = 6, Γ/H ≈ 5e-86. The dilute stage is harmless.
- **Earlier stages.** Γ/H ∝ T⁴ for d = 5. A 1% loss would need a homogeneous high-amplitude phase near
  T ≳ 1e7 GeV, which is outside the model and not part of GRUT's supplied history. Any Affleck–Dine-type origin
  of the X asymmetry would make these operators the charge source.
- **Formed cores.** Denser family members do not radiate much faster per charge; the per-charge rate tracks F₀⁴,
  and F₀ ≤ 2.62. Hot or excited cores could raise the coherent channel by up to the ~400× family spread, and
  stimulated emission is negligible at hidden-sector temperatures ≪ m_X. Exceeding today's 13.8-Gyr exposure
  needs a boost of t₀/t_form ≈ 4e17 (1 s) to 4e11 (1e6 s).

**Formation imposes no stronger bound in the supplied history. Today's object sets it.**

## 8. What quantum gravity suggests for |c_d| [E]

- **Exact global symmetries.** Banks–Seiberg (PRD 83, 084019 (2011)) and Harlow–Ooguri (PRL 122, 191601 (2019);
  arXiv:1810.05338) argue against exact global symmetries in quantum gravity. They supply no coefficient
  (abstracts read).
- **Naive estimate.** O(1) Planck-suppressed operators, the axion-quality argument (Kamionkowski–March-Russell,
  PLB 282, 137 (1992), abstract read), give c₅ ≈ 1. That fails by 5.5e4 in c, or 3e9 in rate.
- **Wormholes.** Breaking goes as e^{−S}. Kallosh–Linde–Linde–Susskind (PRD 52, 912 (1995), abstract read) give
  S ≈ 10 in Einstein gravity, with string theory adding e^{−8π²/g²}. The required S is 9.0 from the coherent
  channel alone and 10.9 with pairs, so plain Einstein-gravity wormholes are marginal and the outcome hangs on O(1)
  action details. Giddings–Strominger (1988) is [R]. The Hebecker–Mikhail–Soler review (Front. 2018,
  arXiv:1807.00824) and Alonso–Urbano (arXiv:1706.07415) were abstract-read.
- **Swampland lower bound.** Daus–Hebecker–Leonhardt–March-Russell (NPB 960 (2020) 115167, abstract read) put
  breaking at no less than e^{−M_P²/Λ²} for Higgsed-gauge-plus-axion U(1)s. That is a lower bound on violation,
  for a different symmetry class.
- **Discrete gauge symmetry.** Krauss–Wilczek (PRL 62, 1221 (1989); Crossref metadata and abstract summary)
  show that discrete gauge charges survive black holes and wormholes. With Z₆ the X⁵ operator is absent exactly
  and X⁶ passes with any c ≲ 1e13. Codex's Σ†X⁶ route gives c₆ ≈ v_Σ/M. Banks–Seiberg treat Z_p gauge theory
  explicitly.
- **Q-ball prior art.** Kawasaki–Nakatsuka (JCAP 04 (2020) 017, arXiv:1912.06993, abstract read) treat Q-ball
  decay through a U(1)-breaking A-term. Q-ball formation and fragmentation follow Coleman,
  Kusenko–Shaposhnikov, Enqvist–McDonald and Kasuya–Kawasaki [R].

## 9. Limits

- The pair rate is a free-Born/LDA estimate. It needs a distorted-wave (Bogoliubov) calculation summed over
  partial waves, including anti-X X, XX, Xs and the d = 6 anti-X s and ss channels.
- The background is frozen, with no evolved depletion and no population average over Q.
- Lower-dimensional operators are assumed absent, as Codex states.
- The classical time domain checks only the coherent channel and is spherical (ℓ = 0).
- The quantum-gravity estimates are order-of-magnitude.

## Reproduce

From `research/reality_framework/grut_formation_review_20261025/tracks/xcharge/code/`:

```
python3 profile.py 200000
python3 forced.py ../results/forced_Q200000.json          # ~30 s
python3 pairs.py ../results/pairs_Q200000.json
python3 scan.py ../results/scan_ground_family.json          # ~5 min
python3 timedomain.py 5 0.00025,0.0005,0.001 ../results/timedomain_d5.json      # ~8 min
python3 timedomain.py 6 0.0005,0.001,0.002 ../results/timedomain_d6.json
python3 timedomain.py 5 0.0005 ../results/timedomain_d5_fine.json 0.03 0.0075 420
python3 summary.py
```

The large-η logs came from an earlier single-sign version of the `timedomain.py` driver. They are kept as
qualitative evidence of the instability only.
