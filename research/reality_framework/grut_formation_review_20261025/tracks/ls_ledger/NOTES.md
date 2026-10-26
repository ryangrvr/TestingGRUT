# Track ls_ledger: injection ledger of Codex b448f25 (low-scale candidates 10 and 17)

26 October 2026 · Claude, issue #13 · same-session AI work, so this is not outside replication.

Packet under review: b448f25 (`origin/codex/grut-low-scale-debris-20261025`,
`research/reality_framework/grut_low_scale_debris_20261025/`). It was read only through `git archive b448f25` into the
scratchpad. No Codex code was run. Codex's JSON was read as data, for comparison only. Nothing was committed.

**Scope.** This is one supplied tree-level working model. Nothing here is evidence that X, Y, Σ or the vector exist. It
is not a validation, an exclusion or a GRUT prediction. This track checks the injection side only. The case-L
thresholds (R_L = 103.0 and 93.4) are quoted [Q] and were not rechecked here.

Labels: [D] derived, [C] computed, [E] estimate, [Q] quoted.

## Method (all our own code)

- **Releases [C].** Our sp_yledger integrator (`tracks/sp_yledger/code/yrelease.c`) and our own multi-start minimum
  (`ymodel.minimum`).
  - It uses cosmic time, physical fields and Dormand–Prince 5(4). Codex uses conformal fields and DOP853.
  - c_H = 10, H_i = 434867798206.6179 GeV, aligned radial minimum, zero velocity, instant release.
  - Ensemble phases: δ = δ0 + N(0, 0.001), drawn with numpy `default_rng(2026102510 / 2026102517)`. This reproduces
    all 48 published member phases exactly.
  - Runs: members 0, 3, …, 45 (16 per candidate), the nominal δ0, and controls (a = 3000; rtol 1e-12 and 1e-10).
  - Run with `code/ls_release.py`.
- **Late projection [C].** `code/ls_project.py` is a separate Python projection of the a = 1000 state. It uses the
  complex-inverse forced Y response and Y modes in the frame co-rotating with arg Σ/2.
  - The Σ radial decomposition continues the undamped release for 20 Σ periods.
  - It then fits δr = r − h/√2 to a static offset, the free m_σ line and the Y-forced lines (2ω_H, 2ω_L, ω_H ± ω_L).
- **Ledger [D, C].** `code/ls_ledger.py`.
  - Widths come from our own formulas. Y_H → Y_L A uses Goldstone equivalence.
  - Decay average: ⟨a⟩ = Γ(3/2)√(2H_i/Γ) (sp_critic). It is checked by direct quadrature with the actual birth epoch
    and, for vectors, the two-stage σ → A → X X̄ chain with time dilation.
- **Decay-timing model [E].** `code/ls_damped.c` + `ls_damped.py`. This is yrelease.c plus one radial friction term,
  Γ_σ ṙ, on Σ (the vacuum width), with counters for the decayed σ number:
  - an upper count, energy removed / m_σ;
  - a lower count, energy removed / ω(E), where ω(E) is an anharmonic frequency estimate.

## Results versus Codex

Our value is first, Codex's second. Ours are 16-member means; "subset" means Codex's own values for the same 16 members.

| Quantity | Candidate 10 | Candidate 17 | Verdict |
|---|---|---|---|
| Ensemble net ε | 0.1574113 ± 7e-6 / 0.1574128 (subset 0.1574113); nominal 0.15741241 / 0.15741241 | 0.1470609 ± 5e-5 / 0.1470462 (subset 0.1470609); nominal 0.14707361 / 0.14707361 | CONFIRMED [C] |
| Raw yield / Y_Q | 5.69004 / 5.69010 | 5.30997 / 5.30942 | CONFIRMED [C] |
| N_YH, N_YL | 16.778, 22.185 / 16.780, 22.184 | 17.814, 24.441 / 17.811, 24.446 | CONFIRMED [C] |
| N_σ radial | 24.979 / 24.977 (subset 24.979) | 36.625 / 36.653 (subset 36.627) | CONFIRMED [C] |
| Same-member maximum deviation | ε 1.2e-7, N_X 8e-12, E_X 3e-11, N_Y 2e-5, N_σ 9e-10 | same order | CONFIRMED [C] |
| Decay-averaged R_Y | 113.07 / 113.07 | 122.80 / 122.81 | CONFIRMED [C] |
| R_Σ, two vectors per late σ | 86.34 / 86.34 | 126.79 / 126.88 (subset sampling) | CONFIRMED as computed [C] |
| Conservative R (Y + Σ) | 199.41 / 199.40 | 249.59 / 249.69 | CONFIRMED arithmetic; lower-bound status: see below |
| Central R | 228.42 / 228.42 | 280.43 / 280.53 | CONFIRMED [C] |
| Pairs per net charge (central) | 0.04647 / 0.04647 (9.43 × the 1% allowance) | 0.06143 / 0.06146 (12.46 ×) | CONFIRMED [C] |

**Controls [C].**
- Long run, nominal, between the 600–1000 and 2000–3000 windows: ε changes by 2e-8. N_YH drifts 0.05% and N_σ 2e-6.
- The rtol 1e-12 and 1e-10 runs agree with the main runs to ≤ 1e-4 in N_Y and 1e-8 in ε.
- The snapshot and the window means agree to ≤ 7e-5.

## Σ radial subtraction [C]: what is counted is the free oscillation

At a = 1000 the continuation fit splits the formula count N_σ = a³[ṙ² + λ(r² − h²/2)²]/m_σ as follows (candidate 10;
candidate 17 is alike):

| Part of N_σ | Value |
|---|---|
| Formula count | 24.97946 |
| Free m_σ line | 24.97950 |
| Static offset | 3.7e-5 quanta (δr = −3.4e-7; this is the anharmonic −(3/2r0)⟨δr²⟩ shift) |
| Y-forced lines | 1e-8 |
| Predicted centrifugal shift | 1.1e-10 |
| Predicted Y-forced static shift | 1.6e-8 |

The fit residual is 5e-4. The charge-carrying angular motion is correctly excluded (a⁴E_ang = 1.93, a negligible
share of E_X = 33077). **Codex's subtraction is CONFIRMED.**

## Disjointness and decay averaging [D, C]: CONFIRMED

**Widths.** Ours and Codex's agree to ≤ 1e-12. Codex's λ^{3/2} form of Y_H → Y_L A equals the Goldstone form
μ²p/(4π m_H²).

**Decay averages.**
- ⟨a⟩_L = 43299.7 and ⟨a⟩_A = 4.1434e11, both equal to Γ(3/2)√(2H_i/Γ).
- Starting the L clock at a = 1 or at a = 23 changes ⟨a⟩ by ≤ 2e-7.
- Vector deposit with time dilation:

  | Case | Deposit / m_A⟨a⟩_A |
  |---|---|
  | Born at a = 10 from σ | 1 + 7e-9 |
  | Born at a = 1000 from Y_H | 1 + 3.7e-5 |
  | Full σ-decay-averaged chain | 1 + 2e-8 |

**Ledger structure.**
- Each Y is counted once, at m_L. The Y_H mass excess goes only into its cascade vector.
- Vector birth momentum redshifts and is not added to the parent.
- Σ's angular charge is excluded.
- E_X is a separate Hamiltonian piece.
- No energy or quantum is counted twice.

## Main finding [C, E]: the late σ and Y_H counts are not the populations that decay

**Codex's assumption.** DECAY_RESULTS assumes "late projected numbers conserved to decay".

**Where it fails [C].** In our nominal histories N_σ and N_YH settle within 10% only for a ≳ 67–75, and within 1% only
for a ≳ 120–180. The physical decays happen far earlier:

| Population | Mean decay epoch ⟨a⟩ | Surviving to a = 67 |
|---|---:|---:|
| σ | 15.6 | 5e-7 |
| Y_H cascade | 23.2 | 1e-3 |

Before freeze-out the instantaneous Σ count swings between 17 and 735 (median 222–306 for a = 5–40).

**Where it holds.** The total N_Y is stable within 10% from a ≈ 17 and within 3% from a ≈ 40. Y_L decays at a ≈ 4e4, so
R_Y is well defined.

**Damped model (vacuum Γ_σ as radial friction; 4 members per candidate) [E].** The ratio is decayed σ number over late
undamped N_σ.

| Quantity | Candidate 10 | Candidate 17 |
|---|---|---|
| Decayed / late N_σ (ω(E) count to m_σ count) | 2.63–6.28 | 2.99–6.29 |
| Decays at a ≤ 17 | 87% | 78% |
| Backreaction: ε | 0.1574 → 0.1508 (−4.2%) | 0.1471 → 0.1412 (−4.0%) |
| Backreaction: N_Y | −6.2% | −9.1% |
| R_Σ | 227–543 | 379–797 |
| R_Y with backreaction | 106.0 | 111.7 |
| Conservative R | 333–649 | 491–909 |

**Sensitivity [E].** Γ_σ × 1/3 and × 3 still give decayed/late ratios ≥ 1.9 and a conservative R ≥ 276. At × 3,
N_Y falls by 28%.

**Reading.**
- The conservative label is plausible in this model: the true Σ injection is larger.
- It is not demonstrated by Codex's late-count method.
- The model is phenomenological in the anharmonic, nonadiabatic stage where most decays occur. Goldstone production
  from a large-amplitude Σ oscillation may be nonperturbative.

**Cascade vectors.** Y_H is likewise repopulated by coherent H–L mixing before freeze-out. The bracket [0, N_Y] for
cascade vectors is therefore not a strict upper bound. This affects only the central and upper ledgers.

## Margins against Codex's quoted raw-density thresholds (103.0, 93.4) [Q thresholds]

| Ledger | Candidate 10 | Candidate 17 |
|---|---|---|
| Conservative R / threshold | 1.94 | 2.67 |
| Y alone | 1.10 | 1.31 |
| Y alone, damped model | 1.03 | 1.20 |

- The Y-only failure of candidate 10 rests on a thin margin.
- The Y + Σ failure is robust, and our model makes it larger.
- The backreaction lowers ε by about 4%, which would also lower R_L. That is not computed here.

## Open items (not certified)

1. A kinetic or coupled treatment of σ → A_L A_L and Y_H → Y_L A_L during a ≈ 5–70. This includes nonperturbative
   Goldstone production and the backreaction on ε and N_Y. Our friction model is [E].
2. The case-L threshold and its dependence on the reduced ε. This was not checked in this track.
3. Primordial vectors: R ≈ 0.011 is quoted [Q], not rechecked (sp_critic confirmed the evolving-mass method).
4. Release-induced Goldstone production (Codex's open item as well).

## Reproduction (about 25 min on 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/ls_ledger
python3 -B $T/code/ls_release.py $T/results/ls_release.json 3 extras
python3 -B $T/code/ls_damped.py  $T/results/ls_damped.json
git archive b448f25 research/reality_framework/grut_low_scale_debris_20261025 | tar -x -C <scratch>
python3 -B $T/code/ls_ledger.py $T/results/ls_release.json $T/results/ls_damped.json \
        <scratch>/research/reality_framework/grut_low_scale_debris_20261025 $T/results/ls_ledger.json
```

The binaries are built in `$SCRATCHPAD/review/ls_ledger/build`. Nothing is written into the sp_yledger track.
