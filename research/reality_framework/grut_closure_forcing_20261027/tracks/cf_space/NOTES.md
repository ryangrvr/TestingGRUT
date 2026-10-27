# Track cf_space: closure in space and capacity

27 October 2026. Claude, parallel pass for issue #13 (brief: `../../BRIEF.md`). Same-session AI work; not an
outside replication. Git read-only; nothing committed by this track. Codex's packets were read through `git show` only.

Labels: [D] derived, [C] computed, [Q] quoted from a source read in this pass, [E] estimate, [A] argued,
[R] recalled. Every comparison with existing data is **retrospective**. Nothing here validates or excludes a model of
the real world, and no new theory is claimed. Each postulate below is an **added choice**, not a consequence of the
word "closed" (FORMAL_CORE §4 row "compact spatial topology … SCOPE"; Codex FOUNDATION_GATE, 2908e91).

Files: `code/cfs_common.py` (constants, sources), `code/p2_capacity.py`, `code/p1_space.py`,
`code/derivations.py` (sympy checks); `results/p2_capacity.json`, `results/p1_space.json`,
`results/derivations.json`, `results/summary.json`.

## 0. Answer

No postulate in this family gives GRUT content that is both distinctive and testable.

- **Spatial compactness (P1)** forces only standard global constraints. Each gauged charge must total zero, the
  translations or rotations of a symmetric background act as constraints, flux is quantized, and S³ carries a
  curvature sign.
  - GRUT's homogeneous runs already impose zero hidden charge, because of homogeneity. Compactness adds nothing
    there.
  - In the Higgs phase, the only extra content is a mod-6 condition on the universe's net X and Y number, and
    nobody can observe it.
  - The S³ curvature sign is the one observable consequence. It is standard closed-FLRW cosmology. Data are
    ambiguous on it (retrospective): DESI DR2+CMB puts the forbidden side at 2.1σ, while Planck+lensing+BAO is
    neutral.
- **Finite state capacity (P2)** comes in readings that contradict each other. Closure selects none of them.
  - **Banks–Fischler finite Hilbert space.** Λ > 0 is a retrodiction, and only if a future crunch is also excluded.
    Inflation is capped at about 26 e-folds above the horizon minimum, at any scale. This reading allows the heavy
    window, but it contradicts the q-theory final state Λ_∞ = 0.
  - **TCC.** It forces V^{1/4} < 6.9×10⁸ GeV and r < 2×10⁻³¹. This **excludes the GRUT heavy-mediator window** and
    every GRUT inflationary X-charge branch on record. The exclusion holds even for the most permissive history
    after inflation (H_end < 9.4×10⁹ GeV, against H_i = 1.33×10¹¹) and for a two-stage rescue (dilution ≥ 10⁵⁴).
- **The postulate with the most force acts against GRUT's own branch.** TCC is a published swampland conjecture,
  and the only test that could refute it, a detection of primordial B-modes, cannot tell GRUT apart from anything
  else.

## 1. Inputs (supplied unless marked)

- Heavy spot point [S], from `REVIEW_HEAVYWINDOW.md` and `tracks/hw_gates/NOTES.md`:
  - H_i = 1.3335×10¹¹ GeV is the Hubble rate at release, at the end of inflation;
  - M_Y = 30 H_i, h = 40 H_i, κ = 3.4587×10⁻³, e = 10⁻⁷;
  - c6 = 1.2504×10⁶, with a 1% survival ceiling of 1.486×10¹³.
- Reheating [S], from `REVIEW_REHEATING.md`: within about one Hubble time, T_RH ≈ 2–5×10¹⁴ GeV.
- Scaling family [Q, own track], from `tracks/ls_scaling`: Y ∝ ε c_H^{3/2} H^{3/2}/λ_X and c6 ∝ H⁻² under
  (M_Y, h, μ, H) → s(…).
- Other GRUT inflationary points [S]:
  - scan-point branch, H ≈ 3–6×10¹¹ GeV (`REVIEW_SCANPOINT.md`);
  - gauge point, H = 4.7×10¹³ GeV (`REVIEW_GAUGEPOINT.md`).
- Candidate 1 [S]: m_X = 92.41 MeV and v = 0.9241 GeV. Its one-loop scalar vacuum term is −1.2265×10⁻⁶ v⁴ and is
  uncancelled (vacuum SYNTHESIS §6).
- Codex architecture [S], 04c53f5 THEORY.md and ACTION.md:
  - hidden U(1) with charges X, Y, Σ = 1, 3, 6;
  - T³ is a supplied branch choice;
  - "No … de Sitter information-capacity rule is inserted into this action."
- Cosmology:
  - H₀ = 67.36 km/s/Mpc and Ω_Λ = 0.6847 [R, Planck 2018];
  - A_s = 2.099×10⁻⁹ [R];
  - T₀ = 2.7255 K [R];
  - M_P = 2.435×10¹⁸ GeV, reduced [R].
  - With the DESI DR2+CMB H₀ = 68.17 [Q, arXiv:2503.14738 eq. 21], N_min shifts by −0.01.

## 2. P1: spatial compactness

### 2.1 The postulates

- **P1-S3.** The Cauchy slices are S³ (k = +1).
- **P1-T3.** The Cauchy slices are a flat T³ of side L. This is Codex's branch.

Both are SCOPE items in FORMAL_CORE §4. A spatially infinite FLRW universe also has no exterior. WORKING_THEORY s09
already records that compactness does not imply k = +1, and that closed de Sitter gives Ω_K = −1/(e^{2N} − 1). Neither
point is repeated here.

### 2.2 What compactness forces

1. **Every gauged charge totals zero [D].** Integrate ∇·E = ρ over a compact slice with no boundary: Q = 0.
   - Check: on a periodic 16³ lattice with a random E, Σ div E = −1.1×10⁻¹⁴, machine zero. The same field in an
     open box gives 0.057, the boundary flux [C].
   - For GRUT's hidden U(1): N_X + 3N_Y + 6N_Σ = 0 on every slice (D7 in `results/derivations.json`).
   - For electromagnetism: zero net charge. Data bound a uniform charge asymmetry below 10⁻²⁶ e per baryon [Q,
     Caprini–Ferreira hep-ph/0310066 v2, search level]. This is consistent, and it is also forced by homogeneity
     (item 2).
2. **The homogeneous GRUT cosmology already imposes it, with no help from compactness [D].**
   - For spatially constant fields, ∇·E ≡ 0 (D6). Gauss's law then gives ρ = 0 at every point, on R³ as well as T³.
   - Our runs enforce N_X + 3N_Y + 6N_Σ = 0 to rounding: 6×10⁻¹² (hw_gates), ≤ 5×10⁻¹⁰ (hw_verify) and 7×10⁻¹⁰
     (REVIEW_TIMING).
   - Codex's sentence "the total gauge charge, including Sigma, remains zero on a compact slice" (THEORY.md §5) holds
     without compactness for every homogeneous history.
   - **Nothing new is forced in the homogeneous sector.**
3. **In the Higgs phase the residue is a mod-6 condition [D].** The Σ condensate absorbs charge only in units of 6. A
   compact slice therefore forces (N_X + 3N_Y) ≡ 0 (mod 6) for the whole universe.
   - This allows 6 of the 36 residue classes of (N_X, N_Y) mod 6 [C].
   - On R³ a net Z6 charge would be allowed, detectable only by Aharonov–Bohm phases around Z6 strings.
   - It applies to a total of order the number of dark-matter quanta, so it is **unobservable**.
4. **Flux and holonomy sectors [D].**
   - A charge-1 field is present, so ∫_{T²} F ∈ (2π/e)ℤ on every 2-cycle.
   - Each Z6 string carries flux 2π/(6e), so the net number of wrapping strings through each 2-cycle is ≡ 0 (mod 6).
   - Flat Z6 bundles on T³ give 6³ = 216 holonomy sectors.
   - At the heavy spot point strings do not form, because the U(1) stays broken through inflation (hw_gates gate 3).
     The trivial sector is therefore consistent and nothing is observable.
5. **Symmetries of a symmetric background become constraints [R].** These are the linearization-stability
   conditions of Moncrief (J. Math. Phys. 16, 493 (1975)) and Fischer–Marsden–Moncrief (1980), and Higuchi's
   quantum version for de Sitter (Class. Quantum Grav. 8 (1991)).
   - On T³ FLRW, the total momentum of all perturbations vanishes.
   - On S³ FLRW, the six SO(4) charges vanish.
   - These conditions are unobservable from inside.
6. **Quantum gravity on a closed slice, in a supplied framework [Q].**
   - With the gravitational path integral, a closed universe without an observer has a **one-dimensional** Hilbert
     space. Adding an observer gives an effective Hilbert space whose dimension is set by the observer's entropy.
     Source: Harlow–Usatyuk–Zhao, arXiv:2501.02359, read at search level.
   - The de Sitter static patch with an observer is Type II₁ and has a maximum-entropy state (empty dS)
     (Chandrasekaran–Longo–Penington–Witten, arXiv:2206.10780, abstract).
   - This is the closest thing in the literature to "compactness forces finite capacity". It needs the gravitational
     path integral as a supplied input, and its errors are suppressed by e^{−S_obs}. It is **untestable** and fits
     Codex's relational architecture.
7. **The S³ curvature sign [D].** Ω_K < 0 strictly, the volume is 2π²R_c³, and the spectrum is discrete. For T³:
   Ω_K = 0 exactly and k_min = 2π/L.
8. **One door opens, nothing is supplied [Q, FORMAL_CORE COND 5.3].** The Fewster–Verch natural-state no-go uses
   noncompact Cauchy surfaces. A category restricted to compact slices escapes the theorem as proved. That only
   *permits* a covariant state rule; it does not supply one.

### 2.3 Discriminators and data (retrospective)

Distance to last scattering: χ_rec = 3.116 c/H₀ = 13.87 Gpc [C].

**Curvature, Ω_K.**

| Dataset (source) | Ω_K | P(Ω_K < 0), the S³ side | One-sided 95% R_c lower limit if closed |
|---|---|---|---|
| Planck 2018 TT,TE,EE+lowE (1807.06209 eq. 46b) [Q] | −0.044 +0.018/−0.015 | 0.996 | 17 Gpc |
| + lensing (eq. 47a) [Q] | −0.0106 ± 0.0065 | 0.949 | 31 Gpc |
| + lensing + BAO (eq. 47b) [Q] | 0.0007 ± 0.0019 | 0.356 | 90 Gpc |
| DESI DR2 + CMB, ΛCDM+Ω_K (2503.14738 Table 5) [Q] | 0.0023 ± 0.0011 | **0.018 (2.09σ)** | none (95% lower limit > 0) |
| DESI DR2 + CMB + SNe, w0wa+Ω_K (Table 5) [Q] | 0.0006 to 0.0011, ± 0.0013 | 0.20–0.32 | 114–138 Gpc |

The S³ sign "test" depends on the dataset and the model, and gives answers in opposite directions. A DR1 reanalysis
reports that the PR4 HiLLiPoP/LoLLiPoP likelihoods remove the DESI deviation [Q, search level, arXiv:2511.20757].
**The current data contain no fixed discriminator.**
- The cosmic-variance floor is about 1.5×10⁻⁵ at 1σ (0804.1771); the practical floor is about 10⁻⁴ (0901.3354,
  1604.01410). Both are [Q, search level].
- Inflation with N ≥ N_min + 5 drives |Ω_K| below that floor.

**Topology.**
- Matched circles find no pairs. The shortest loop *through us* is longer than 0.985 d_LSS = 6.14 c/H₀ [Q, COMPACT
  Part I, arXiv:2211.02603, JCAP 01 (2023) 030].
- Loops that avoid us can be shorter by a factor of 2 to at least 6 [Q, same source].
- Planck 2015 XVIII (1502.01593): cubic torus R_i > 0.97 χ_rec at 99% CL, i.e. L > 6.04 c/H₀ [Q].
- COMPACT (2210.11426, PRL 132, 171501): topology might still be detectable for scales somewhat larger than the
  horizon [Q, abstract].

### 2.4 Foundation gate for P1

| Item | P1-S3 | P1-T3 |
|---|---|---|
| 1. Postulate | Cauchy slices S³ | Cauchy slices flat T³ (side L supplied) |
| 2. Compatibility | GR + gauge theory on a compact slice: Q = 0, SO(4) charges zero; consistent | same, total momentum zero; consistent |
| 3. Forbidden baseline | infinite flat or open FLRW (Ω_K ≥ 0); any state with net gauge charge or net Z6 charge | R³, S³, H³ (any Ω_K ≠ 0); same charge statements |
| 4. Discriminator and comparator | sign of Ω_K. Comparator: DESI DR2+CMB 0.0023 ± 0.0011 (2.1σ against, retrospective) and Planck+lens+BAO 0.0007 ± 0.0019 (neutral) | Ω_K = 0 exactly (shared with flat ΛCDM); matched circles and correlation-matrix searches, needing L < ~6.1 c/H₀ |
| 5. Failure condition | Ω_K > 0 at ≥ 3σ in a combination robust to the CMB likelihood, with the central value above the 10⁻⁴ floor | a robust detection of Ω_K ≠ 0 above the floor, or a detected non-T³ topology |
| Supplied | R_c (or L), all couplings and the state | L and the holonomy sector |
| Credit | standard closed FLRW; Gauss law; Moncrief/Higuchi; HUZ/CLPW for QG | Codex's branch choice; Planck/COMPACT searches |
| GRUT-native? | no | no (a stated choice in 04c53f5) |
| Heavy window | no effect: homogeneous Gauss already imposed, no strings | same |
| Candidate 1 | no effect. The flat-FRW premise of its H₀ exclusion holds for any Ω_K > −0.563, which data already satisfy | supplies "flat" exactly, which was already satisfied |
| Dark-matter chain | mod-6 condition on the net X+3Y number; unobservable | same; 216 holonomy sectors; unobservable |
| **Verdict** | **MIXED**: the forced content is STANDARD; the sign is testable only if \|Ω_K\| ≳ 10⁻⁴ and the data are currently ambiguous | **UNTESTABLE**: standard content, no GRUT observable depends on L |

The **hidden-U(1) Gauss reading** (item 2.2.1–4) gets its own row in §5. Its verdict is STANDARD: homogeneity already
imposes it, and its remainder is unobservable.

## 3. P2: finite state capacity

### 3.1 Readings (closure motivates all of them and selects none)

- **P2-BF (Cosmological Complementarity).**
  - The postulate: "the physics of an asymptotically de Sitter space-time is completely described by a pure state in
    a Hilbert space of finite dimension, determined by the value of the cosmological constant" [Q, Banks–Fischler
    astro-ph/0307459 §1, full text read]. Then dim ℋ = exp(S_Λ).
  - Related: Banks hep-th/0007146 [R], Fischler 2000 [R], Bousso's N-bound hep-th/0010252 [Q abstract], Witten
    hep-th/0106109 [R].
- **P2-TCC.** For every expanding phase, a_f/a_i < M_P/H_f [Q, Bedroya–Vafa arXiv:1909.11063 abstract;
  Bedroya–Brandenberger–Loverde–Vafa arXiv:1909.11106]. This is a swampland conjecture, related in spirit to finite
  capacity but not derived from it.
- **P2-weak** [R]. These bounds are not binding here:
  - N ≲ S_dS(H_inf) (Arkani-Hamed et al., arXiv:0704.1814);
  - quantum break time N ≲ M_P²/H² (Dvali–Gomez–Zell, arXiv:1701.08776).

### 3.2 The heavy-window inflation numbers [C], `results/p2_capacity.json`

| Quantity | Value |
|---|---|
| H_i | 1.3335×10¹¹ GeV (H/M_P = 5.48×10⁻⁸) |
| V = 3M_P²H² | 3.163×10⁵⁹ GeV⁴, so **V^{1/4} = 7.50×10¹⁴ GeV** |
| r = 2H²/(π²M_P²A_s) | **2.90×10⁻⁷** (ε = 1.8×10⁻⁸), if H_i equals H at pivot exit. If H falls during inflation, r is larger. |
| Instant reheating T | 3.0–3.1×10¹⁴ GeV, consistent with the stated 2–5×10¹⁴ |
| N_min: today's Hubble scale inside the horizon at the start | **58.41** (58.27 with T_RH = 2×10¹⁴ and a matter era) |
| N_* at k = 0.05 Mpc⁻¹ | 53.0 |
| N for \|Ω_K,0\| ≤ 10⁻³ (10⁻⁴) from Ω_K,i = 1 | 61.9 (63.0) |
| Comparison with BK18 r < 0.036 (2110.00483) [Q] | allows H ≤ 4.70×10¹³ GeV. The heavy window passes by 10⁵ in r. The gauge point at H = 4.7×10¹³ gives r = 0.0360, exactly at the BK18 limit (retrospective). |

The heavy window's r = 2.9×10⁻⁷ cannot be tested by any planned B-mode programme. LiteBIRD's target is
δr < 10⁻³ (arXiv:2202.02773) [R].

### 3.3 Banks–Fischler: the e-fold cap

- **BF's counting [Q, eq. 3.2].** BF give N_e ≤ ¼ ln(Λ_I/λ) + (1/12) ln(M_P⁴/λ). We reproduce their illustrative
  numbers [C]:
  - 86.0 ("~85") for κ = 1;
  - 64.5 ("~65") for κ = 1/3;
  - 88.3 ("~88") for eq. 3.2.
  - Sympy check D3: residual 0.
- **The same counting with explicit O(1) factors [D, check D2].** Take the radiation entropy of the inflated patch
  after reheating, conserved adiabatically afterwards:
  - S_patch = S_H0 e^{3(N − N_min)}, where S_H0 = s₀(4π/3)H₀⁻³ = 3.14×10⁸⁸;
  - requiring S_patch ≤ S_Λ = 3.31×10¹²² gives **N − N_min ≤ Δ = ⅓ ln(S_Λ/S_H0) = 26.11**;
  - this holds **independently of the inflation scale** (any reheating that is adiabatic afterwards).
- **Heavy window [C].** N_max = 84.5 by our version, 84.1–84.7 by eq. 3.2 with the Planck-mass ambiguity. N_min is
  58.4. **The heavy window is compatible, with 26 e-folds to spare.** Candidate 1 has no inflation input and is
  unaffected.
- **Observable consequence [D].** The cap gives only an upper bound on N. Combined with compact space (§4) it gives
  a curvature floor |Ω_K,0| ≥ |Ω_K,i| e^{−2Δ} = 2×10⁻²³ |Ω_K,i|. That is 10¹⁸ below the cosmic-variance floor:
  **unobservable**.

### 3.4 TCC: the inflation-scale bound and the heavy window

- **The bound [D, check D1].** Horizon: e^N ≥ (g_s0/g)^{1/3}(T₀/T_rh)(H/H₀). TCC: e^N < M_P/H. With instant
  reheating,
  V^{3/4} < 3M_P³ (H₀/T₀) (30/π²g)^{1/4} (g/g_s0)^{1/3}.
- **Numbers [C].** **V^{1/4} < 6.89×10⁸ GeV** (6.91×10⁸ with GRUT dof), H < 0.113 GeV, **r < 2.1×10⁻³¹**,
  T_rh < 2.8×10⁸ GeV, N = 44.5 at the bound. BBLV quote "lower than 10⁹ GeV", ε < 10⁻³¹ and negligible r [Q].
- **Heavy window [C].**
  - TCC allows N ≤ ln(M_P/H_i) = **16.72** e-folds, against N_min = 58.41: a deficit of 41.7 e-folds.
  - H_i is 1.18×10¹² times the bound; V^{1/4} is 1.09×10⁶ times the bound.
  - The heavy window's own requirement, reheating within one Hubble time at 2–5×10¹⁴ GeV, puts it on the strict line
    (N_min = 58.3–58.4).
- **The most permissive post-inflationary history [C].** Minimize N_min over an intermediate phase with w ∈ [−1/3, 1]
  until T_rh ∈ [5 MeV, T_inst].
  - The optimum is w = −1/3 down to T_BBN. There aH is constant, so no mode exits.
  - It gives **H_end < 9.36×10⁹ GeV** and V^{1/4} < 2.0×10¹⁴ GeV.
  - r < 1.4×10⁻⁹ if H is constant. If H_* > H_end, r can be larger: Mizuno–Mukohyama–Pi–Zhang (arXiv:1910.02979)
    quote r ≲ 10⁻⁸ [Q abstract].
  - The heavy window's H_i = H_end still exceeds this by **14.2**. It is excluded even here, and this history
    contradicts its required fast reheating anyway.
- **Two-stage rescue [D, check D5].** Stage 1 is the heavy window, with N₁ ≤ 16.7; a later low-scale stage supplies
  the rest.
  - The horizon condition gives e^{N₁+N₂} ≥ (H₁/H₀)(T₀/T_rh1)(T_on2/T_rh2) ≥ e^{N_min,1}, so N₂ ≥ 41.7.
  - The X charge made at the end of stage 1 is then diluted by ≥ e^{3N₂} = **10^{54.3}**.
  - Recovering the yield through Y ∝ H^{3/2} would need H × 10³⁶, above the Planck scale. **Excluded.**
  - Torabian (arXiv:1910.06867) [Q, search level] notes that multi-stage inflation relaxes TCC for the CMB stage. It
    does not rescue a relic made in the first stage.
- **Every GRUT inflationary X-charge branch on record is excluded [C]:**

  | Branch | H (GeV) | N_TCC | N_min | H / H_TCC (standard) | H / H_TCC (relaxed) |
  |---|---|---|---|---|---|
  | heavy spot | 1.33×10¹¹ | 16.7 | 58.4 | 1.2×10¹² | 14 |
  | scan point | 3–6×10¹¹ | 15.2–15.9 | 58.8–59.2 | 2.7–5.3×10¹² | 32–64 |
  | gauge point | 4.7×10¹³ | 10.9 | 61.3 | 4.2×10¹⁴ | 5.0×10³ |

  **Inside the uniform scaling family**, survival needs c6 ≤ 1.486×10¹³, so H ≥ 3.9×10⁷ GeV at fixed ε. The yield
  with ε ≤ 1 needs H ≥ 3.7×10¹⁰ GeV. Reaching H ≤ 0.11 GeV would need ε c_H^{3/2}/λ_X × 1.3×10¹⁸. With λ_X fixed by
  Candidate 1 and ε ≤ 1, that means c_H ~ 10¹², outside the EFT [E].
- **TCC at late times [C].** The present acceleration must end within ln(M_P/H_Λ)/H_Λ = 138.9 × 17.5 Gyr =
  **2.4×10¹² yr**. There is no eternal de Sitter. This is untestable now.

### 3.5 Does finite capacity force Λ > 0? Only with "no crunch" added [D, C]

The countermodel is a k = +1 dust universe with Λ = 0 and Ω_m0 = 1 + δ.
- Its largest causal-diamond 2-sphere is the S³ equator at turnaround, r_max = (c/H₀)(1 + δ)/δ^{3/2}.
- Its covariant-entropy capacity 8π²M_P² r_max² is **finite**: 9.1×10¹²², 2.7×10¹²⁵ and 2.3×10¹²⁸ for
  δ = 1, 0.1 and 0.01.
- A flat universe with Λ < 0 also crunches and has finite capacity.

The forced statement is therefore "Λ > 0 **or** a future crunch". The retrodiction Λ > 0 follows only with BF's
explicit "asymptotically de Sitter". The observed acceleration (since 1998; DESI DR2 still has ρ_DE > 0) is
consistent with it, and is a retrodiction, not a test. The relation Λ = 3π/(G ln N) relabels Λ by N, and N is supplied
[C: S = 3π/(GΛ) = S_dS to 1 part in 10¹⁵]. This is the "N-bound / Banks–Bousso" row already credited in vacuum
SYNTHESIS §2 (class X).

### 3.6 The readings contradict each other and the q-theory branch [D/A]

| Reading | Asymptotic Λ | Eternal dS? | Inflation scale | Heavy window | q-theory P_q (final Λ_∞ = 0) |
|---|---|---|---|---|---|
| P2-BF | > 0, fixed by N | yes, required ("our bound would not apply to a meta-stable dS … vanishing cosmological constant" [Q]) | any; N ≤ N_min + 26.1 | allowed | **contradicts** (P_q's final state is Minkowski, with infinite capacity) |
| P2-TCC | → 0 or crunch | no: lifetime < 2.4×10¹² yr | V^{1/4} < 6.9×10⁸ GeV | **excluded** | compatible in asymptotics |
| GKS (hep-th/0212209) [Q abstract] | finite entropy and exact dS symmetry cannot both hold; a dS phase cannot outlive its recurrence time | no | — | — | — |

- **BF and TCC contradict each other** on late-time de Sitter.
- **BF contradicts P_q.**
- [A] A strict "asymptotically dS forever" reading of BF also forbids terminal decay of our vacuum (decay into a
  crunch ends the dS). That is the absolute-stability condition the q-vacuum packet found on the wrong side of the
  measured m_t(m_h, α_s) at about 3.1σ (HANDOFF item 3: unreviewed, retrospective, SM up to the Planck scale).
  - GRUT's own Higgs portals λ_Hi can move that boundary, so within GRUT this is not a fixed discriminator.
  - GKS instead requires dS to end before recurrence. A finite-capacity dS that is both eternal and exactly symmetric
    is not available.

"Closed" motivates all of these readings and selects none: the strengthening is a choice.

### 3.7 Candidate 1 (low scale, no inflation input)

- **P2-BF supplies the missing vacuum-sign premise [D/A].** The conditional H₀ exclusion of the 100 GeV, r = 1
  preparation (H₀ ≥ 116.66 km/s/Mpc against six masers; CAMPAIGN_LEDGER §13) assumed:
  - flat FRW (any Ω_K > −0.563 suffices);
  - a total vacuum energy ≥ −(2.28 meV)⁴, which Candidate 1's matching does not supply.

  Under P2-BF the asymptotic Λ > 0. If the dark sector rolls down monotonically, the present vacuum term is ≥ Λ_∞ > 0
  [A], so the premise holds by postulate.
  - The exclusion of that one preparation then needs no separate vacuum-matching rule, and the bound
    m_X|Y_Q0/Y_B0| < 63.93 GeV follows.
  - This uses already-exposed data and is **retrospective**; it is not a new test.
  - The one-loop scalar term −8.94×10⁻⁷ GeV⁴ = −3.6×10⁴⁰ ρ_Λ still has to be cancelled by a supplied counterterm [C].
- **P2-TCC: compatible [C].** Reheating is allowed up to 2.8×10⁸ GeV = 3×10⁹ m_X. But every GRUT inflationary
  production route is excluded (§3.4), so Candidate 1's relic charge would need a source after inflation. TCC fixes
  none of Candidate 1's couplings.
- **P1:** nothing new (§2.4).

### 3.8 Foundation gate for P2

| Item | P2-BF | P2-TCC |
|---|---|---|
| 1. Postulate | dim ℋ = e^{S_Λ}, pure state, asymptotically dS | a_f/a_i < M_P/H_f for every expanding phase |
| 2. Compatibility | consistent with semiclassical inflation if N ≤ N_min + 26; needs approximate dS symmetry (GKS) | consistent with low-scale inflation only; forbids eternal dS |
| 3. Forbidden baseline (derived) | eternal Λ ≤ 0 without a crunch; quintessence → 0; the P_q Minkowski final state; N > N_min + 26.1 | high-scale inflation: V^{1/4} ≥ 6.9×10⁸ GeV (standard) or H_end ≥ 9.4×10⁹ GeV (any history); r > 2×10⁻³¹ (H const.); dS older than 139/H |
| 4. Discriminator and comparator | sign of Λ: retrodiction (acceleration observed); e-fold cap unobservable (curvature floor 2×10⁻²³) | primordial B-modes. Comparator: BK18 r < 0.036 (not constraining); future δr ~ 10⁻³ (LiteBIRD) |
| 5. Failure condition | evidence that dark energy decays to ρ_DE ≤ 0 or ends in a crunch | any detection of primordial tensors with r ≳ 10⁻⁸ (MMPZ's most relaxed bound) |
| Supplied | N (equivalently Λ); state; couplings | V(φ), reheating history; couplings |
| Credit / GRUT-native? | Banks, Fischler, Bousso; relabelled ("N-bound" row of vacuum SYNTHESIS) | Bedroya–Vafa; BBLV; MMPZ; not GRUT-native |
| Heavy window | compatible (N_max 84.5 vs N_min 58.4) | **excluded** by 41.7 e-folds; relaxed history ×14 in H; two-stage dilution 10⁵⁴ |
| Candidate 1 | supplies the vacuum-sign premise of the 100 GeV, r = 1 H₀ exclusion (retrospective) | compatible; needs a non-inflationary relic source |
| Dark-matter chain | none | all on-record inflationary X-charge branches (H = 1.3×10¹¹ to 4.7×10¹³) excluded |
| **Verdict** | **MIXED**: new content UNTESTABLE; Λ > 0 a satisfied retrodiction; **inconsistent with the q-theory P_q branch** | **INCONSISTENT WITH GRUT BRANCH** (heavy window and every inflationary DM branch); the conjecture itself is refutable by a B-mode detection, but that is not a GRUT discriminator |

P2-weak gives N_max ≈ 2.6×10¹⁶ (S_dS at H_i) or 3.3×10¹⁴ (M_P²/H²) at the heavy point. It is not binding, and its
verdict is STANDARD.

## 4. P1 + P2 together: curvature and topology from a capped inflation [D, C]

**Exact relations (check D4).** For an initial S³ with Ω_K,i = 1/ℓ², or a T³ of side ℓ/H_i:
- |Ω_K,0| = |Ω_K,i| e^{−2(N − N_min)};
- L₀H₀ = ℓ e^{N − N_min}.

A capacity cap N − N_min ≤ f therefore gives a curvature **floor** and a topology **ceiling**. [A] Applying P2 to the
whole closed universe at the start of inflation (dim ℋ ~ e^{S_dS(H_i)}) motivates ℓ ~ 1. That value is still
supplied.

- **BF.** f = 26.1, so the floor is 2×10⁻²³ |Ω_K,i| and the ceiling is 2×10¹¹ ℓ. Both are **unobservable**.
- **TCC.** f_TCC(V) = ln(M_P/H) − N_min ≈ (3/4) ln(V_TCC/V). With ℓ = 1:
  - S³: |Ω_K,0| ≥ 10⁻³ for V^{1/4} ≥ 2.19×10⁸ GeV, and ≥ 10⁻⁴ for V^{1/4} ≥ 1.49×10⁸ GeV.
    - DESI DR2+CMB would put Ω_K = −10⁻³ at 3.0σ and −3×10⁻³ at 4.8σ. Planck+lensing+BAO gives 0.9σ and 1.9σ.
    - The curvature becomes observable only in a band V^{1/4} ≈ (1.5–2.2)×10⁸ GeV, i.e. 0.22–0.32 of the TCC limit.
  - T³: a null topology search (L₀ > 6.14 c/H₀) requires V^{1/4} ≤ 3.78×10⁸ GeV (0.55 of the limit).
  - **No fixed discriminator.** V is free, and TCC itself makes it unmeasurable through r. ℓ is supplied.
  - **Verdict: UNTESTABLE.** It is a correlated statement in an unobservable parameter.

## 5. Verdict table

| ID | Postulate | Verdict |
|---|---|---|
| P1-S3 | spatial slices S³ | MIXED (STANDARD content; sign weakly testable; data ambiguous, retrospective) |
| P1-T3 | spatial slices T³ (Codex branch) | UNTESTABLE |
| P1-Gauss | compact-slice Gauss law for the hidden U(1) (1, 3, 6) | STANDARD (homogeneity already imposes it; mod-6 remainder unobservable) |
| P1-QG | compact slice + gravitational path integral ⇒ one-dimensional closed-universe ℋ, observer-relative | UNTESTABLE (relabelling of HUZ/CLPW) |
| P2-BF | finite Hilbert space e^{S_Λ}, asymptotically dS | MIXED (UNTESTABLE new content; inconsistent with the q-theory P_q branch) |
| P2-TCC | Trans-Planckian Censorship | INCONSISTENT WITH GRUT BRANCH (heavy window and all inflationary DM branches) |
| P2-weak | N ≲ S_dS(H_inf), quantum break time | STANDARD (not binding) |
| P1+P2 | compact space + capacity-capped inflation | UNTESTABLE |

## 6. Limits

- The BF full text was read. These were read at abstract or search level only: TCC, BBLV, MMPZ, Torabian, GKS,
  Bousso, CLPW, HUZ, COMPACT, Planck 2015 XVIII, Caprini–Ferreira, the curvature-floor papers and BK18.
- Recalled: Banks hep-th/0007146, Witten hep-th/0106109, Moncrief, Higuchi, ADNTV, Dvali–Gomez–Zell, LiteBIRD.
- N_min uses H constant during inflation and instant (or one-Hubble-time) reheating. Changing either moves N_min by
  ≲ 1, and no verdict depends on that.
- The heavy-window exclusion under TCC holds for H_end. If the release occurred well before the end of inflation, the
  relevant H_end is lower, but the yield scaling would then change too (not computed).
- The [A] links are argued, not derived:
  - BF ⇒ no terminal decay;
  - P2 ⇒ ℓ ~ 1;
  - a monotone dark sector ⇒ present vacuum term ≥ Λ_∞.
- Reproduce with `cd code && python3 -I p2_capacity.py && python3 -I p1_space.py && python3 -I derivations.py`. Each
  takes a few seconds.
