## A. Classical charge floor (independent derivation)

Derived from the Candidate 1 action without reading Codex's DERIVATION.md first.

With X = (q₁ + iq₂)/√2 and the (−+++) canonical action, the proper energy density for any
observer is ρ = ½φ̇² + |Ẋ|² + ½Ḋ² + (non-negative gradient terms) + U, and the U(1) charge density
is n_Q = q₂q̇₁ − q₁q̇₂ (sign convention immaterial). Cauchy–Schwarz gives |n_Q| ≤ 2|X||Ẋ|, so
|Ẋ|² ≥ n_Q²/(4|X|²). The D terms are non-negative at the benchmark (m_D0 = g_X v, λ_D > 0), so
ρ ≥ n_Q²/(4s²) + V(s) with s = |X| and V(s) = min_φ U(φ, s, D = 0). By AM–GM,
ρ ≥ |n_Q| · min_s √(V(s)/s²) ≡ m* |n_Q|. Minimizing U/s² over (φ, s): on the φ = 0 branch the
minimum is √(λ_φλ_X) v², at s² = v²√(λ_φ/(4λ_X)); on the φ = v branch the infimum is g_X²v², approached as
s → 0; interior stationary points never go lower. Hence m* = v · min(g_X, (λ_φλ_X)^{1/4}).

| Check (`code/review.py`, R1) | Result |
|---|---|
| m*/v, closed form vs direct 2D minimization | 0.056234132519 both |
| μ = m*/m_X; m* | 0.562341325190; 56.23413252 GeV |
| pointwise ρ ≥ m*\|n_Q\| on 800,000 random field/velocity/gradient samples | 0 violations |
| saturation by the homogeneous rotating condensate at φ = 0, ω = m* | residual 0 |

**Disposition: PASS** for the classical, minimally coupled, canonical floor. It agrees with Codex
to every printed digit.

## B. Flat-FRW floor and the six-maser envelope (independent implementation)

Pesce et al. 2020 Table 1 was transcribed here from the arXiv v2 PDF
(`data/pesce2020_table1_claude_transcription.json`), not from Codex's JSON. Velocities are
CMB-frame, optical convention, as the table note states. ω_b = 0.02233 ± 0.00036 (Mossa et al.
2020), m_b = 0.938 GeV with the 1% allowance, corners as in Codex's frozen protocol.

| Quantity | This review | Codex |
|---|---:|---:|
| ω_h,min (nominal) | 1.3387080801 | 1.3387080801 |
| H_floor nominal / conservative (km/s/Mpc) | 116.66353672 / 113.25181752 | same |

| Galaxy | D_low corner | D_A,max primary | D_A,max stress | Primary | Stress |
|---|---:|---:|---:|---|---|
| UGC 3789 | 39.50 | 35.20 | 41.58 | mismatch | ok |
| NGC 6264 | 81.10 | 90.97 | 97.00 | ok | ok |
| NGC 6323 | 40.40 | 72.08 | 78.23 | ok | ok |
| NGC 5765b | 96.90 | 77.84 | 83.95 | mismatch | mismatch |
| CGCG 074-064 | 66.00 | 67.03 | 73.21 | ok | ok |
| NGC 4258 | 7.25 | 12.55 | 19.07 | ok | ok |

Necessary product bound m_X|Y_Q0/Y_B0| < 63.9344 GeV (primary), 74.6429 GeV (stress): identical
to Codex. Convention checks: reading the velocities as radio convention changes the NGC 5765b
stress envelope from 83.95 to 85.92 Mpc (still a mismatch); using the atomic mass unit 0.9315 GeV
for m_b changes H_floor from 113.25 to 113.64 (tightens).

**Disposition: PASS** for the arithmetic and the decision as stated. This is a second-agent,
retrospective review of public data; the conditions below are scope, not defects.
