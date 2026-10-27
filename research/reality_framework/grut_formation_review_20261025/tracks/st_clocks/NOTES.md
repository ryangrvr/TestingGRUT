# Track st_clocks: independent decay clocks and kinematics for Codex b73de64 (source timing), issue #13

27 October 2026 · Claude · this is AI review work in the same session, not outside replication.

Packet reviewed: `grut_source_timing_20261026` at b73de64, read from the git archive in the scratchpad. I read
MATHEMATICS.md, REPORT.md, PROTOCOL.md, CLAUDE_HANDOFF.md, `code/source_clocks.py` and `code/timing_common.py`.
**No Codex code was run or imported.** Codex's SOURCE_CLOCKS.json and RESPONSE_RESULTS.json were read only as data:
they supply the comparison targets and the sampled crossing times. Every number marked "own" comes from this track's
code:

- `code/background800.py` extends **our own** reheating engine (`tracks/rp_reheating`, rpreheat.c, the engine that
  reproduced Codex's e160356 runs) to τ = H_i t = 800 at rtol 1e-10, with a tight control at 1e-12.
- `code/clocks.py` re-derives the widths from the couplings and computes the clocks on our own a(τ). The source
  dictionary `kin_source/results/source_interface.json` is our own file, and Codex's `ks_common.SRC` loads that same file.

Scope: decay clocks and kinematics only. There is no nonlinear evolution, no population, no f_sym, and no claim of
validation or exclusion. The model is a supplied one with supplied fields. Labels: [D] derived, [C] computed,
[Q] quoted, [E] estimate, [A] argued.

## 1. Widths and their physics [D, C]

The spot-point couplings come from rp.py and hw_gates spotmodel.py [Q own]: H_i = 1.3334962e11 GeV, M_Y = 30,
h = 40, μ = 7.159456, κ = 3.4586847e-3 and λ_Σ = 1. The potential terms are 2Re(κY*X³) and 2Re(μΣ*Y²), the charges
are X:Y:Σ = 1:3:6, and Σ = h/√2.

- **Y_L.** Write Y = (Y₁+iY₂)/√2. The μ term splits the masses to m² = M_Y² ± √2μh.
  - The light mode is the real direction with arg Y² = arg Σ + π (for μ > 0). Its mass is m_L = 22.2486 H_i.
  - Each real mode couples as (κ/√2)(X³ + X*³), up to a phase. Per channel (XXX and X̄X̄X̄) this gives
    |M|² = (3!κ/√2)² = 18κ², with an identical-particle factor 1/3!.
  - Massless three-body phase space: own Dalitz quadrature / (m²/256π³) = 1 − 1e-16.
  - Total: Γ_L = 3κ²m_L/(256π³), split half to each channel. **Γ_L = 13413.6819 GeV, γ_L = 1.005903e-7.** Codex and
    the dictionary agree to 2e-16. The heavy mode has the same formula with m_H; it cascades first, see §5.
  - Conditions: X is effectively massless here. It is in the φ = 0 phase, the condensate mass is ≪ m_L late, and the
    vacuum m_X is irrelevant. |Σ| stays at its vacuum value (hw_gates: to 1e-4).
- **Vector.** The Σ kinetic term with charge 6 and Σ = h/√2 gives m_A = 6eh. A → XX̄ is one channel (X is a single
  charge-1 complex scalar, distinct particles), with Σ|ε·(p₁−p₂)|²/3 = 4p*²e²/3.
  - Width: Γ_A = e²m_A β³/(48π). **Own: 2.1223251e-10 GeV (e = 1e-7) and 2.1223251e-7 GeV (e = 1e-6), equal to Codex
    and the dictionary to 2e-16.** β³ − 1 = −5e-15.
  - Late X/Y condensate contributions to m_A are negligible: e²|X|² against (6eh)² is about 0.1/a².
  - A → YȲ and A → σ are closed (m_A ≪ m_Y, m_σ). Branching to XX̄ is 1 [Q gp_leftovers].
- **Scaling with e.** Γ_A = e³h/(8π) at fixed h, so Γ ∝ e³. The own ratio of medians is 1000.0000000000. This makes
  sense only because h is held fixed while e changes.

## 2. Medians and the radiation clock [C]

| Quantity | Codex | Own | Verdict |
|---|---:|---:|---|
| τ½, Y_L | 6.890794e6 | 6.890794e6 | CONFIRMED |
| τ½, A at e = 1e-6 / 1e-7 | 4.355172e17 / 4.355172e20 | same to 5e-15 | CONFIRMED |
| C = a²H/H_i, (3, 0.5) | 1.9858052758 | **1.9858052754** (rel −2.1e-10) | CONFIRMED |
| C, (5, 1) | 1.4758677865 | **1.4758677875** (rel +7.0e-10) | CONFIRMED |
| Drift of C over τ 64–800 (max−min)/C | 4.0e-11, 1.3e-10 (Codex "range") | 7.3e-12, 7.7e-12 | CONFIRMED (tighter) |
| max \|a² − 2C(τ−τ*)\|/a², τ 64–800 | — | 8.2e-12, 1.1e-11 | a² = 2C(τ−τ*) is exact to roundoff |
| τ* (matched at 800) | 0.193674103, −0.169061875 | 0.193673989, −0.169061502 (fit 64–800: same to 1e-9) | CONFIRMED |
| a(800) | 56.36061784, 48.59925395 | 56.36061784, 48.59925395 (8 significant figures) | CONFIRMED |

C settles to within 1e-6 of its late value by τ ≈ 22 on (3, 0.5) and τ ≈ 11 on (5, 1); within 1e-3 by τ ≈ 9.7 and
4.6. The tight run changes C by 6e-15. At τ = 20, C on (3, 0.5) is 1.985796, 4.5e-6 low (the inflaton remnant).

Physical validity of extending the radiation tail [A]: a²H stays constant while g* stays constant, and kin_source
gives T_vis > 9 TeV through every listed decay. The hidden matter components (Y_L, A) enter the Friedmann equation
suppressed by η² = (H_i/M_r)² = 3e-15 and stay negligible until long after their decays. So a matched C is the right
clock out to the vector epochs. C = 1 is not.

## 3. Scale factors at the medians [C]

| Parent | Codex (3,0.5) / (5,1) | Own (3,0.5) / (5,1) | Verdict |
|---|---|---|---|
| Y_L | 5231.40 / 4509.97 | 5231.40 / 4509.97 (rel −1.1e-10 / +3.5e-10) | CONFIRMED |
| A, e = 1e-6 | 1.31518e9 / 1.13381e9 | same | CONFIRMED |
| A, e = 1e-7 | 4.15897e10 / 3.58543e10 | same | CONFIRMED |

The 1%, 10%, 90% and 99% number quantiles and the 10%, 50% and 90% energy quantiles also agree. The worst case is
3.5e-10 relative in a and 1e-13 relative in τ.

## 4. Fractions at the response samples [C]

These use the own a(τ) (cubic spline on Δτ = 0.05, checked with Simpson at grid-aligned τ) and the own tail
(mpmath upper incomplete gamma, checked against direct mpmath quadrature from τ = 800 to ∞).

| Sample | number fraction Codex / own | energy fraction F_E Codex / own |
|---|---|---|
| (3,0.5), 0.03 rad unit, τ = 134 (a = 23.0527 both) | 1.34790e-5 / 1.34790e-5 | 3.716111e-8 / **3.716112e-8** (rel 2.3e-7) |
| (5,1), 0.03 rad unit, τ = 147.5 (a = 20.8777 both) | 1.48370e-5 / 1.48370e-5 | 4.306670e-8 / **4.306671e-8** (rel 1.9e-7) |
| 1e-6 rad unit, τ = 772 / 782 | 7.7653e-5 / 7.8659e-5, same | 5.14581e-7 / 5.24964e-7, same (rel 1.6e-8) |
| Vector fractions (all 36 samples) | same | same (rel ≤ 4.7e-7) |

- Simpson agrees with the spline to 1e-9 relative. The 1.6–4.7e-7 relative offset from Codex is interpolation
  error, from Codex's PCHIP on its own log-spaced grid (steps up to 3.2 in τ). Codex's quoted "orders 4 and 8 agree to
  3.4e-16" tests only the quadrature of that interpolant. The difference is immaterial at the quoted digits.
- Closed tail: ∫_t^∞ γ√(2C(s−τ*))e^{−γs}ds = √(2C/γ)e^{−γτ*}Γ(3/2, γ(t−τ*)), derived by the substitution
  u = γ(s−τ*). Closed form against quadrature agrees to double precision. **CONFIRMED.**
- Late-decay limit: total mean a = 5568.653 (3,0.5) and 4800.711 (5,1). √(πC/2γ) differs by 1.9e-8 and 1.7e-8
  (vectors: < 1e-15). **CONFIRMED.** (Old C = 1: 3951.68.)

## 5. Other prescribed sources before τ ≈ 150–800, and the ordering statement [C, A, E]

Enumeration of the source dictionary (`source_interface.json` = Codex `ks_common.SRC`), with own widths on own
backgrounds:

| Source | Own width γ = Γ/H_i | Timing on the quiet histories | What it deposits |
|---|---|---|---|
| Y_H → Y_L + A_L (Goldstone equivalence, coupling √2μ G Y₁Y₂) [E width] | 0.035042 (4.673e9 GeV, = dictionary) | Γ = H at τ = 14.5 / 14.1 (a = 7.5 / 6.5); τ½ = 19.8; **99.1% / 99.4% done by τ = 134 / 147.5** | one Y_L plus one longitudinal vector of E = 11.2 H_i, sterile until its own decay at τ ~ 1e17–1e20; **no X-sector heat** |
| σ → A_L A_L / σ → Y_L Y_L | 1.125 / 0.0111 | prompt (τ ~ 1) | N_σ = 3e-9 (dictionary), 3e-10 late here: negligible |
| Y_L → XXX, X̄X̄X̄ | 1.0059e-7 | 1% at τ = 9.99e4 (a = 630 / 543); median 6.89e6 | **first X-sector heat in the dictionary** |
| A → XX̄ (primordial and cascade vectors) | 1.59e-18 / 1.59e-21 | medians 4.36e17 / 4.36e20 | the dominant energy (dictionary R = 4.05 / 40.6) |
| Inflaton | — | decays only into the visible radiation fluid in both engines (own rpreheat.c; Codex reheat.c/reheating.py read) | no hidden branch is prescribed; the hidden fields feel it only through gravity and ξR\|X\|², already in the computed history |

- **Verdict on "linear control is lost well before substantial prescribed decay heating": CONFIRMED, qualified to
  X-sector heating.**
  - The earliest prescribed X-sector heat source is Y_L → 3X. Its 1% number point is τ ≈ 1.0e5, 130–750 times later
    than the unit samples.
  - The one prescribed decay that does act before τ ≈ 150 is the Y_H → Y_L + A_L cascade, which is essentially
    complete by then. It feeds sterile vectors and Y_L and puts no heat into X.
  - The visible reheating bath is also complete then. Its contact with X (through φ and the Higgs portal,
    κ_portal ≲ 1e-9 [Q hw_gates]) is not prescribed and not computed.
  - "Half the mediators decay around H_i t = 6.89e6" is correct for Y_L's rest clock.
- **Not in the dictionary [E]: a condensate-assisted Y_L → XX / X̄X̄ (one leg on the X condensate).**
  - Rate: Γ₂ = 9κ²|X|²/(8πm_L), with |X| taken from our run. Γ₂ = 1.8e-6 at τ ≈ 134, which is 18 times the
    three-body rate (Γ₂/H ≈ 5e-4).
  - Integrated from τ = 1, it gives a Y_L decayed fraction of 1.1e-3 by τ = 134 and 1.5e-3 by τ = 800. That is about
    80 times the three-body 1.35e-5.
  - sp_yledger had already flagged this channel and neglected it.
  - It is still ≪ 1, so the ordering does not change. It means the quoted 1.35–1.48e-5 is the dictionary channel
    only, not a bound on all Y_L decays.
- **Population note [C] (diagnostic; Codex correctly does not transfer normalisations).**
  - On these histories the homogeneous late counts are N_YL = 0.00682 / 0.0178 and N_YH = 9.7e-6 / 4.3e-6. The
    dictionary has 0.334 and 0.0173.
  - With E_X = 1.052e5 / 6.46e4, that population would give R_Y = 0.0080 / 0.029, against the dictionary's 0.806.
  - Injection up to the unit sample would be 3.0e-10 / 1.3e-9 of E_X (4.2e-8 with the dictionary population).
  - Quantum Y production is absent from a homogeneous run, so this is not a source normalisation.

## 6. Wavelength cube [D, C]

Units [D]. N_X = a³n_X/H_i³ is the classical Noether charge, with ħ = 1, counted in the comoving volume whose
physical size was H_i⁻³ at the a = 1 start. The rp engine's NX and Codex's q₁p₂ − q₂p₁ (with q = √2az) are both
a³·2Im(X*Ẋ). A comoving k (in H_i units) has physical wavelength 2πa/k. So Q_cube = N_X(2π/k)³ exactly, and it is
independent of a while N_X is conserved.

| Quantity | Codex | Own |
|---|---:|---:|
| N_X, 64–80 normalisation | 8156.643 / 5236.091 | 8156.735 / 5236.084 (τ-weighted window; rel 1.1e-5) |
| N_X at τ = 800 | 8116.07284 / 5209.71567 | same to 3e-12 (drift −0.50%) |
| Q_cube at k = 2.13 | 2.0937e5 / 1.3440e5 | **2.0937e5 / 1.3440e5** (with N_X at 800: 2.083e5 / 1.337e5) |
| Q_cube at k = 5.8 | 1.0370e4 / 6.657e3 | same |

**CONFIRMED** as a geometric number. The choice of cell moves it by an order of magnitude: a half-wavelength cube
gives 2.6e4 / 1.7e4, a sphere of diameter λ gives 1.1e5 / 7.0e4. A single cosine mode makes sheets, not cubes. Codex
already labels it geometric and not a fragment charge; that label is necessary.

## 7. Verdict table (`results/verdicts.json`)

| # | Claim | Codex | Own | Verdict |
|---|---|---|---|---|
| 1 | Γ_YL formula and value | 1.34137e4 GeV | 1.34137e4 GeV (own derivation and phase space) | CONFIRMED |
| 1 | Γ_A formula and values; Y_L identification | 2.1223e-10 / 2.1223e-7 GeV | same | CONFIRMED |
| 1 | 1000× median ratio for 10× e | 1000 | 1000.0000000000 (Γ ∝ e³ at fixed h) | CONFIRMED |
| 2 | Medians | 6.89e6, 4.36e17, 4.36e20 | same (≤ 5e-15) | CONFIRMED |
| 3 | C | 1.986, 1.476 | 1.98580528, 1.47586779 (own engine; drift 7e-12) | CONFIRMED |
| 4 | a at medians | 5.23e3/4.51e3, 1.32e9/1.13e9, 4.16e10/3.59e10 | same (≤ 3.5e-10) | CONFIRMED |
| 5 | Number fractions at unit samples | 1.35–1.48e-5 | 1.3479–1.4837e-5 | CONFIRMED |
| 5 | Energy fractions at unit samples | 3.72–4.31e-8 | 3.7161–4.3067e-8 (rel ≤ 2.3e-7) | CONFIRMED |
| 5 | Closed tail; mean a ≈ √(πC/2γ) | formulas | quadrature exact; 1.9e-8 | CONFIRMED |
| 6 | Q_cube, N_X units | 1.34–2.09e5 | 1.3440–2.0937e5; N_X = a³n/H_i³ count | CONFIRMED (geometric only) |
| 7 | "well before substantial prescribed decay heating" | statement | earliest X heat Y_L → 3X (1% at τ = 1.0e5); Y_H cascade 99% done by τ 134 but no X heat; [E] condensate-assisted Y_L 1e-3 by τ 134 | CONFIRMED, qualified |

No CORRECTED or REFUTED items.

Not checked here (other tracks):
- the response crossings themselves (τ = 134 / 147.5 / 772 / 782);
- the direct-vector adaptation and the phase-current sign;
- the quartic Floquet spectrum.

## 8. Limits

- All widths are tree-level vacuum widths. No medium effects (Bose enhancement in the X condensate, plasma masses)
  and no thermalisation are included. The Y_H cascade width relies on Goldstone equivalence [E].
- The rest clock assumes parents present at τ = 0. Relativistic vectors make early vector numbers smaller (Codex
  says this too).
- The condensate-assisted channel is an order-of-magnitude estimate. Its |X| is the homogeneous amplitude, with no
  phase space or frequency corrections.
- The background is the homogeneous quiet history only. Inhomogeneous release, which is the point of Codex's
  conclusion, would change both the Y counts and N_X.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
cd research/reality_framework/grut_formation_review_20261025/tracks/st_clocks/code
RP_SCRATCH=$SP/review/st_clocks/rp python3 -B background800.py   # ~6 s; writes results/background800.json
python3 -B clocks.py                                             # ~10 s; writes results/clocks.json
python3 -B verdicts.py                                           # writes results/verdicts.json
```
