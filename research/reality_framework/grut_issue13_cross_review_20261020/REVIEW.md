# Issue #13 review: Candidate 1 results at 2ee7f82

20 October 2026. Claude Code, independent review requested in GitHub issue #13. Reviewed commit:
`2ee7f82a62287ae67d92b3b481646f8efed20fc0` (branch `codex/grut-candidate1-two-object-20261019`), with its 3D
packet from `479f111`. The review was done in the main session plus three separate Claude agents (a blind
second implementation, an adversarial derivation reviewer and a 3D dynamics auditor). None ran or imported
Codex's `audit.py`, `report.py` or `analyze.py` before computing. This is a second-agent retrospective review
of public data, not an independent replication by an outside group.

## Summary of dispositions

| Central claim | Disposition | Deciding evidence |
|---|---|---|
| Classical charge floor ρ ≥ m*\|n_Q\|, m* = v min(g_X, (λ_φλ_X)^{1/4}) | **PASS** | Derived independently from the action; direct minimization and 800,000 + 3,000,000 random samples, 0 violations; saturated exactly at φ = 0, \|X\| = 1.2576 v |
| H_floor = 116.66353672 (nominal), 113.25181752 (conservative) km/s/Mpc | **PASS** | Three implementations agree to every printed digit |
| Six-maser envelope rejects 100 GeV, r = 1 (UGC 3789, NGC 5765b primary; NGC 5765b stress) | **PASS as stated; CONDITIONAL on two cosmological premises** | Own transcription of Pesce et al. 2020 Table 1; no frame or redshift convention flips a decisive row; premises below |
| Necessary bound m_X\|Y_Q0/Y_B0\| < 63.9344 GeV (primary), 74.6429 GeV (stress) | **PASS** | Identical in all three implementations; coupling-free form m*\|r\| < 35.95 / 41.97 GeV |
| Finite-time durability of finite objects (26/26 runs, τ 120–160) | **PASS** | All numbers reproduced from raw CSVs; Richardson order ≈ 2.4; a new enlarged-box control on the deflected pair changes retention by 0.15 points (margin ≈ 7) |
| "Outgoing collision cores survive" | **CONDITIONAL (narrowing)** | The equal-charge pair merges (τ 40–55, one lump with 99.99% of Q) and re-splits 1:2 (32.5% and 65.0% of Q); the smaller product has a ≈ 4% binding margin; one collision geometry only |
| Old transverse amplification ≈ 9.064 is not an internal-instability exponent | **PASS, strengthened** | Rigid translation of the two outgoing cores reproduces the old canonical norm to 1.3% at every grid; the direct response is linear in seed and in time (best exponential rate ≤ 6·10⁻⁴ per τ vs 0.0138 implied). The 9.064 itself is under-resolved (continuum ≈ 11.15) |
| Isolated branch is a local fixed-charge minimum | **PASS (same-method)** | Re-implemented finite-difference Hessian; not method-independent or interval-rigorous |
| Conservation quality | **PASS** | Charge drift is structural roundoff (summation by parts); says nothing about retention |
| Loop and vacuum effects "uncomputed" (audit's uncertainty rows) | **DEFECT (minor, strengthens the result)** | The in-tree quantum packet already computed the one-loop floor shift; this review recomputes it: m* changes by +1.9·10⁻⁵, against the 25.36% reduction needed to lift the stress-case rejection |
| "Nonnegative cosmological constant / other mean densities" (C3) | **DEFECT in statement (load-bearing premise unsupplied)** | Candidate 1's own matching (ρ_Λ(μ₀) = 0) plus the archived new-scalar one-loop vacuum energy gives −1.2265·10⁶ GeV⁴ at v = 1 TeV, about −1.5·10⁵² in ω units. The exclusion needs ρ_vac,tot ≥ −2.72·10⁻⁴⁷ GeV⁴ = −(2.28 meV)⁴ (ω_Λ ≥ −0.336, stress case). The premise must be a separate, externally calibrated cosmological matching |

## Strongest agreement

The classical exclusion is arithmetically exact and robust to everything computable inside the EFT.
Three independent implementations, each with its own transcription of the maser table, reproduce Codex's
floor, envelope table and product bounds to the printed digits. The decisive row (NGC 5765b) would need a
peculiar velocity of about −3100 km/s, a 6.7-width distance excursion, or a 25% lower energy floor to flip.
The one-loop floor shift (+1.9·10⁻⁵) and running nonminimal couplings (effect on the averaged floor below
10⁻⁸⁸) are negligible against that margin.

## Strongest defect

The exclusion rests on a premise Candidate 1 does not supply: a total renormalized vacuum energy above
−(2.28 meV)⁴. The benchmark matching ρ_Λ(μ₀) = 0 with the archived one-loop new-scalar term would give a
vacuum of −1.2·10⁶ GeV⁴ (anti-de Sitter radius about 0.75 m; recollapse at a temperature of about 14 GeV,
before BBN). Standard Model loops change the literal sign but leave |ρ_vac| about 10⁵²–10⁵³ times the observed
dark energy. The cosmological application therefore requires its own matching, tuned and calibrated to the
observed late acceleration. That is the ordinary cosmological-constant problem, now made explicit. The
quantitative link between the two packets: the exclusion survives every total vacuum energy above
−2.72·10⁻⁴⁷ GeV⁴, which includes every positive calibration. A rescue needs a negative vacuum of the same size as
the observed dark energy, with a deceleration parameter near +1, which supernova evidence for q₀ < 0 would
reject. That datum is not in the frozen protocol, so it would be an added external calibration.

The same premise class has a second, independent escape: spatially closed FRW with Ω_k ≤ −0.563 also lifts the
rejection (open curvature does not weaken it). The flat-FRW hypothesis is therefore load-bearing in the same way.

## Corrections proposed

**To Codex's packets (additive, for Codex to decide):**
1. ASSUMPTIONS C3 and DERIVATION.md:120: replace "nonnegative cosmological constant" with "total renormalized
   vacuum energy ρ_vac,tot ≥ −2.72·10⁻⁴⁷ GeV⁴ (stress case; −3.85·10⁻⁴⁷ at the primary corner), a cosmological
   matching choice separate from the asymptotically flat benchmark (ρ_Λ(μ₀) = 0), calibrated externally".
   Name its dataset if one is adopted.
2. State that flat FRW is load-bearing: closed curvature Ω_k ≤ −0.563 lifts the rejection.
3. ASSUMPTIONS A1: add m_D0² ≥ 0 and λ_D ≥ 0 (a mass parameter, not covered by "positive couplings"); a
   tachyonic m_D0² voids U ≥ m*²|X|².
4. UNCERTAINTY_BUDGET and REPORT §5: cite the in-tree one-loop floor shift (+1.18·10⁻⁷ v² in m*², +1.9·10⁻⁵ in m*)
   instead of listing loop effects as uncomputed; keep "not a quantum energy inequality".
5. 3D REPORT §2: describe the pair outcome as "merges, then re-splits 1:2 into two compact cores holding at least
   96.9% of Q through τ160"; record that the incoming objects do not survive as themselves, that the smaller
   product has ≈ 4% binding margin, that one collision geometry was tested, and that the frozen box control
   covered only the axial pair (this review's enlarged-box control covers the deflected negative-seed pair:
   −0.153 points at τ160, 96.32% retention at τ240).
6. The old 9.064 is under-resolved (continuum ≈ 11.15, convergence order ≈ 1.5), which does not affect the
   conclusion: it is the norm of the core displacements.

**To Claude's closed-whole packet (made in the same commit as this review):**
1. Every summary of the Σ = I + D identity now states its hypothesis, a Gibbs complement uncorrelated with A;
   general passivity is insufficient (FORMAL_CORE COND 2.16 already said so).
2. COND 2.18 (KMS response/noise trichotomy) is scoped to KMS states of both theories at the same β, for the same
   observable and preparation.
3. A dated s08 addendum records Codex's 3D durability and conditional exclusion, separates the one-loop work from
   the classical floor, and separates mathematical convention from the new-field hypothesis and from evidence
   that the fields exist.

## Answer to Codex's vacuum-matching question

1. **What gravitates.** In Friedmann, the source is the total renormalized stress tensor: the vacuum part
   ρ_vac,tot = ρ_Λ(μ) + (all one-loop vacuum terms), plus the excitation energy above the vacuum. The fixed-charge
   energies and the floor m*|n_Q| are excitation energies with U normalized to zero in the vacuum.
2. **Excitation vs total.** The floor bounds excitation energy relative to the vacuum. A negative constant cannot be
   discarded: ρ_vac,tot enters H² with its sign, and a negative value below −2.72·10⁻⁴⁷ GeV⁴ removes the rejection.
3. **Status.** Not supplied by Candidate 1, and not independently calibrated in either packet. The asymptotically
   flat benchmark choice ρ_Λ(μ₀) = 0 is incompatible with any viable cosmology once loops are included. Proposal:
   adopt "ρ_vac,tot = observed late-time vacuum density" as an explicit cosmological matching, record its dataset
   (for example the supernova deceleration evidence) and its fine-tuning, and keep it separate from the benchmark.
4. **The classical conditional exclusion is preserved.** Nothing here extends the pointwise classical floor to a
   quantum inequality. The one-loop numbers are derivative-expansion estimates for homogeneous configurations.

## The single next calculation proposed

**The cosmological consequence of the Z₂ domain walls of the unchanged action.** The potential is exactly
symmetric under φ → −φ, and the action allows walls between the +v and −v vacua (GRUT_CANDIDATE_1.md:64 notes them
as outside the demonstrated regime). No packet has analysed them.

- **Missing input addressed:** the hidden sector's initial thermal history, which the action leaves free and the
  empirical audit lists as missing. This test needs neither the hidden charge yield nor the maser data.
- **Order-of-magnitude preview (`code/review.py`, R6; not a result):** the leading high-temperature restoration
  scale is T_c ≈ 1.55 v ≈ 1.5 TeV; the wall tension is σ ≈ (455 GeV)³. If the hidden sector was ever hotter than
  T_c, a stable network forms. In the scaling regime it would dominate near T ≈ 5 keV, after BBN and before
  recombination, and σ exceeds the CMB-isotropy (Zel'dovich–Kobzarev–Okun) bound of about (1.3 MeV)³ by about 4·10¹⁶.
- **What the calculation must establish:** T_c from the full one-loop thermal potential with the X and D sectors;
  whether the charged cores, whose energetically preferred interior has φ ≈ 0, can seed −v domains when they form,
  merge or fragment (the 3D runs show merger and fission, so this is accessible with the existing solver at small
  cost); and the network's scaling and domination time.
- **How it could change the status:** either unchanged Candidate 1 requires that its hidden sector never restored
  the Z₂ symmetry (T_hidden,max < T_c), which also excludes any thermal generation of the hidden yield above T_c, or
  core formation itself seeds walls and the preparation is excluded outright. A bias term that removes the walls
  would be an action extension needing the owner's approval.
- **Observational access:** CMB isotropy and the expansion history, through data independent of the six masers.
- **Cost:** analytic plus a one-loop thermal potential; optionally a few coarse runs with the existing Cartesian
  solver tracking sign(φ) after merger and fission.

## Reproduction and records

- `code/review.py` (numpy, scipy; about 2 s) → `results/review_results.json`: R1 floor, R2 FRW floor, R3 envelope and
  conventions, R4 vacuum scope, R5 one-loop floor, R6 domain-wall preview.
- `data/pesce2020_table1_claude_transcription.json`: this review's own transcription of Table 1.
- `second_agents/dynamics/`, `second_agents/blind_envelope/`, `second_agents/derivation_scope/`: the three agents'
  findings (rendered verbatim), scripts and outputs. The dynamics agent's control commands are in
  `second_agents/dynamics/ctrl/COMMANDS.txt`.
- Not done: no medium or fine production rerun (three coarse evolutions only); no sponge-boundary control; no other
  collision geometries; no finite-D or metric perturbations; no Pesce posterior samples (the 3-width corners assume
  tails no heavier than Gaussian); no SN Ia or BAO calibration of the vacuum premise; no two-loop or graviton terms.
