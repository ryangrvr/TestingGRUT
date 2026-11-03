# lf_fermion: the reheating closure from fermions on the lattice

Claude, lattice-fermion closure pass (BRIEF_LF.md; GitHub issue #13 context). Same-session AI work on one supplied
model. Every method, cutoff, reduction, seed and closure here is a model choice; none follows from closure. Nothing
here is a validation, an exclusion, a new theory or a GRUT prediction. Labels: [D] derived, [C] computed, [Q] quoted
from a track, [E] estimate, [A] argued / assumption, [R] literature (arXiv numbers verified by search this session).

## 0. Q4, written before any closure number or window was computed

Question (BRIEF_LF Q4): does any principle become testable only once the reheating closure is known?

Rule used: a principle counts here only if its *own prediction* (a value of c0, δ or s_*) needs the closure-dependent
history (K, f_rad(t), the lattice-stage m_s²(t)) as an input. A principle whose value is closure-independent but whose
*comparison* with the windows moves with the closure does not count; that is every candidate already in
`tracks/pr_principle` (B, D, E, F, and the pre-registered H1).

Listed now, before any window is computed:

1. **Q4-a, branch-edge criticality: "the declared e sits at the branch edge, e = e_crit(c0, δ)".** It gives a line in
   (c0, δ), which needs e_crit, which needs the closure. **Warning, stated in advance:** e_crit is defined by the
   S/Q = 3.47 threshold, so this principle is branch-defining in kind, like candidate F; a hit would not count as a
   non-tuned hit.
2. **Q4-b, unique self-consistency: "the realised (c0, δ) is where the branch-aware fixpoint changes from two roots
   to one" (the edge of bistability, as at P0 c0 = 1.40).** It needs the closure through H_i(c0). It is again built on
   the relic rule's step at S/Q = 3.47, so it is branch-defining too.
3. **Q4-c, release synchrony: "the hidden release ends when the coherent stage ends" (t_release = t_coh,end).** Both
   times need the closure (the coherent-stage end is what this track computes). It gives a c0 only if t_release
   depends on c0 monotonically; that is not known here. It is not defined by S/Q.

No other principle was found whose own prediction needs the closure. **K itself does not depend on c0** (the holder
sector does not back-react on reheating), so "minimal K" or "maximal reheating efficiency" selects nothing in
(c0, δ) [A]. The closure is not a principle; Fermi statistics decides it (BRIEF_LF).

Who evaluates them: lf_chain, with unchanged definitions. This track only supplies the closure.

## 1. Regression R1, homogeneous (done first) [C]

`code/lf1d.py` (line solver) and `code/regress_homog.py` → `results/regress_homog.json`. Constant amplitude, no
expansion, units m = 1, Δk = 0.04, dt = 0.004, spectral free streaming (no doublers), exact 2×2 mass kicks.

| case | ε(1 half-period): line / mf_hybrid table | F_sat (mean of half-periods 7–16): line / skeptic |
|---|---|---|
| χ→N2, A = 42.3 | 1.0007 | 0.1310 / 0.1310 |
| χ→N3, A = 42.3 | 1.0002 | 0.0840 / 0.0840 |
| χ→N4, A = 42.3 | 1.0001 | 0.0582 / 0.0582 |
| χ→N2/N3/N4, A = 20 | 1.0034 / 1.0010 / 1.0004 | 0.0368 / 0.0274 / 0.0185 against 0.0369 / 0.0274 / 0.0185 |
| χ→N2/N3/N4, A = 10 | 1.0015 / 1.0005 / 1.0003 | 0.0094 / 0.0102 / 0.0070 against 0.0094 / 0.0102 / 0.0070 |
| σ→N1, A = 12.34 | 1.0002 | 0.1358 / 0.1358 |

No normalisation was fitted: the k > 0 and k < 0 halves of the line are the two helicities of the 3D radial problem,
and ε_3D = (1/2π²)∫k²E n(k) dk over the whole line. The whole F_tot(n) sequence agrees too (file).

## 2. Method for the inhomogeneous background (what was actually computed)

**Scales first [D/C].** The non-perturbative production fills the band |k| < μ5,max = A·m, with A = 42 at the slice
and local A ≈ 5–40 later. Along the recorded lines the comoving a·μ5 is ≈ 10⁴ (rms 3–8×10³ after fragmentation, hot
spots up to 5×10⁴); the χ-energy correlation length is ≈ 0.01–0.03 H⁻¹ (physical). So the band fermions have
k·ℓ_frag ≈ 30–300: their wavelength (≈ 2π/μ5 ≈ 10⁻³ H⁻¹) is far *below* the fragment size, not above it. The brief's
"0.05 H⁻¹" is the perturbative χ → N₄N₄ momentum (m/2·β), not the band. The P0-vs-GM question is therefore not
wavelength against fragment, but (i) how the band reacts when the local amplitude changes, and (ii) how far band
fermions stream (they move at ≈ c with a zig-zag) during the band memory of a few half-periods (≈ 0.04–0.08 H⁻¹),
which is longer than ℓ_frag [E].

**Method [C on stated choices].** `code/lf1d.py`, `code/lfline.py`.
- Unrotated basis: free massless chiral components (exact spectral streaming on a periodic line: no doublers, the
  SLAC-type derivative is harmless for a classical background), mixed locally by the complex mass a M e^{iθ(s,η)},
  θ = vφ/f, exact 2×2 kicks, Strang splitting in conformal time; comoving grid 4096 → 8192 sites (k_Nyq = 6.1×10⁴)
  at t = 0.34; time step from the no-aliasing condition dη(2k_Nyq + max aθ̇) ≤ π.
- **Fermi statistics:** the many-body state is the evolved Dirac sea; it is sampled by n_r = 4 random-phase
  superpositions per line (low-cost stochastic fermions, Borsanyi & Hindmarsh, arXiv:0809.4711 [R]); the estimator of
  C(k,k) is exact for a homogeneous background and has O(n_k) noise per mode otherwise.
- **Background:** σ, χ, σ̇, χ̇ on 12 lattice lines (4 per axis) of a new 64³ run (`code/lf3d.c`, a copy of
  e2e_lattice's lat3d.c, which is unchanged; it reproduces e2e_lattice P0's f_rad to ≤ 0.3%), recorded every lattice
  step; cubic Hermite in t and band-limited interpolation in s. Non-back-reacting.
- **Ray reduction [A]:** each line is a 1+1D problem (fermions moving along it, both helicities through k ≷ 0); the
  transverse gradient ∇⊥θ (helicity mixing, |b⊥|/k ≈ 0.1–0.3 at the band) and transverse streaming are dropped.
  The 3D energy is the k²-weighted line occupation, exact for a homogeneous background (R1).
- **Closures compared on the same lines and the same background:** passive tracers inside lf3d of the family
  bl = 1 − [(1−λ)B(x) + λ⟨B⟩]/[(1−λ)F_sat(E(x))E(x) + λF_sat(⟨E⟩)⟨E⟩] at λ = 0 (P0), ¼, ½, ¾, 1 (GM), with no N decay
  and B redshifting as a⁻⁶ (as the fermions have no decay either).

## 3. Non-back-reacting result on the P0 background (`results/cmp_REC_P0.json`) [C]

Energy in each species (line mean over 8 χ lines / 4 σ lines, ρ/v² units) against the closures on the same sites:

| t | Q/P0: N2 / N3 / N4 / N1 | Q/GM: N2 / N3 / N4 / N1 | homogeneous control Q_hom/closure_hom: N2 / N3 / N4 / N1 |
|---|---|---|---|
| 0.30 (before fragmentation) | 1.61 / 2.14 / 1.90 / 1.97 | 1.61 / 2.13 / 1.91 / 2.01 | 1.54 / 2.15 / 1.93 / 1.92 |
| 0.50 | 2.07 / 2.01 / 2.20 / 2.15 | 4.8 / 3.9 / 3.8 / 6.0 | 1.58 / 2.10 / 2.48 / 1.72 |
| 0.70 | 3.37 / 2.64 / 2.42 / 3.18 | 13.4 / 8.9 / 6.9 / 23.8 | 1.58 / 2.10 / 2.48 / 1.70 |
| 1.10 | 3.64 / 2.86 / 2.58 / 4.20 | 15.3 / 10.2 / 7.7 / 27.1 | 1.58 / 2.11 / 2.48 / 1.67 |

Line-to-line standard error of Q: ≤ 1% before t = 0.4, ≈ 10% after (N1: 30–40%, only 4 lines and hot spots).
The unblocked tracking closure (TR) lies above Q by ×1.3–1.9.

**Reading.**
1. **Even with no inhomogeneity the closures under-produce by ×1.5–2.5** (homogeneous control on the line-mean field
   of the same history). At constant amplitude they are exact (R1), so the deficit is the changing amplitude: when A
   falls, the fermions above the new band edge no longer block, but the closure keeps counting their energy in B
   against a capacity that shrinks as A⁴ [A, mechanism]. This is an error common to P0 and GM.
2. **Fragmentation adds a further ×1.5–2.5 over P0** (×2–4 in total), and ×4–27 over GM. The fermions are produced
   in the local band, but the local amplitude changes on the fragment time scale and the band fermions stream away
   (§2), so blocking bites even less than in the local closure.
3. **So the computed production is not between P0 and GM.** In the output-level interpolation
   λ_out = (E_P0 − E_Q)/(E_P0 − E_GM) it is **λ_PB ≈ −1.8 … −3.8** after fragmentation (−2.2 to −3.5 for the χ
   species at t ≥ 0.6; N1 −1.8 → −3.8), on the far side of P0, towards unblocked tracking. Inside the closure-level
   family above, λ = 0 (P0) is the extreme and no λ ∈ [0, 1] reproduces it.
4. **It is not unique; it depends on the history.** Before fragmentation λ is undefined (P0 = GM) while Q/closure is
   already 1.5–2.2; after it λ_out moves from about −2 (t = 0.5) to −3.5 (t = 1.1), and it differs by species
   (N1 most, N4 least). H_i dependence: see the 8.4×10¹⁰ paragraph below.

**A better description than λ_PB [C on a model choice]:** local blocking with the Pauli capacity multiplied by κ_j(t)
(`code/calib.py`, `calib_seq.py`, `kappa_table.py`). κ = 1 is P0, κ → ∞ unblocked tracking. Fitted window by window:
κ ≈ 2–4 before fragmentation (homogeneous control: 1.7 / 2.8 / 4.1 / 3.7 for N2 / N3 / N4 / N1) and
κ ≈ 10–90 in the fragmented stage. Adopted table (`results/kappa_REC_P0.json`): κ_hom = 1.71 / 2.81 / 4.10 / 3.67
(t ≤ 0.45), κ_frag = 24 / 16 / 18 / 29 (t ≥ 0.55, geometric mean of the window fits 0.55–1.1, held after 1.1 [A]);
brackets from Q ± 1 s.e.: κ_frag = 17 / 10 / 10 / 7 (lo) and 48 / 30 / 60 / ≥ 10² (hi). An effective coarse-graining
length is *not* the better description: no smoothing between the local and the box-mean capacity can exceed P0.

**Convergence of the fermion side [C]** (`results/conv_fermion.json`, line 1 or 0, t ≤ 0.8): time step cfl 0.5 → 0.25:
≤ 0.4%; fermion grid 8192 → 16384: +2–5% (+12% for N2 at t = 0.5, transient); realisations 4 → 8 (new seed):
−8% … +2% per line, consistent with the line-to-line scatter; low-k cut and projection vs chiral measure: ≤ 3%.

**H_i dependence [C, 2 χ + 2 σ lines only]** (`results/cmp_REC_H84_P0.json`, `calib_seq_REC_H84_P0.json`; 8.4×10¹⁰
GeV background, fermion grid 16384 because a·μ5 reaches 1.2×10⁵ in σ hot spots): Q/P0 = 1.7–2.0 at t = 0.3,
2.5–3.3 (χ species) at t = 1.1; λ_out = −1.8 … −3.3 (N1 −7.4 at t = 1.1, ±60%, one hot line). In the window fits κ
rises later (≈ 2–4 until t ≈ 0.65, 8–30 from t ≈ 0.75), following the later fragmentation (e2e_lattice: f_coh < 0.1
at 0.74 instead of 0.49). So κ depends on the time *since fragmentation*, roughly not on H_i.

## 4. The calibrated closure in the drain, rerun to the end of the coherent stage (t = 3) [C on a model choice]

Closure **LF** = local blocking with capacity κ_j(t)·F_sat(E(x))E(x) (`lf3d.c` option `kappa`, table in
`results/kappa_REC_P0.json`, linear switch 0.45 → 0.55), Γ_th = 0, perturbative tails free, as P0 otherwise. 64³,
L = 0.42, same seed and spectra as e2e_lattice P0.

| run | f_rad 0.3 / 0.5 / 1.0 / 3.0 | f_grad max | t(f_coh < 0.1) | f_−v max | K_cont (stage2, κ = 1) | K estimate [E] |
|---|---|---|---|---|---|---|
| e2e_lattice P0 [Q] | 0.17 / 0.28 / 0.60 / 0.604 | 0.148 | 0.49 | 0.018 | 3.59 | 8.77 (e2e_chain) [Q] |
| e2e_lattice GM [Q] | 0.17 / 0.17 / 0.17 / 0.248 | 0.221 | 0.47 | 0.031 | 5.65 | 15.28 (e2e_chain) [Q] |
| e2e_lattice TR (unblocked) [Q] | 0.54 / 0.73 / 0.83 / 0.857 | 0 | never | 0 | 2.05 | 4.3 |
| **LF (KFmid)** | **0.28 / 0.38 / 0.81 / 0.811** | 0.042 | 0.85 | 0.003 | 2.64 | **5.4** |
| LF lo (κ from Q − s.e.) | 0.27 / 0.36 / 0.78 / 0.783 | 0.053 | 0.75 | – | 2.72 | 5.7 |
| LF hi (κ from Q + s.e.) | 0.28 / 0.42 / 0.84 / 0.845 | 0.032 | never | – | 2.50 | 4.9 |
| LF, second seed | 0.28 / 0.38 / 0.81 / 0.811 (+0.07%) | 0.043 | 0.82 | – | 2.63 | 5.4 |
| LF at 8.4×10¹⁰ (κ from the H84 fit, switch 0.65–0.75) | 0.42 / 0.44 / 0.65 / 0.832 | 0.004 | never | – | 2.69 | 6.1 |

- The K estimate is e2e_chain's χ-weighted rule K ∝ (comoving χ + gradient energy at t_sw = 3)^{2/3}, anchored on its
  P0 (8.77) and GM (15.28) hand-overs; it reproduces GM from P0 to 2% (15.6 vs 15.28) [E]. **The real K must come
  from lf_chain's kinetic hand-over**, because the χ → N₄N₄ self-blocking of the remnant lives in the network.
- Stronger drain → weaker fragmentation (f_grad ≤ 0.05 instead of 0.15), fewer −v overshoots (0.3%), no walls.
- Lattice-stage spectator mass m_s² (HFa, 0.05-window mean): 8.42 / 8.93 / 8.98 / 8.99 at t = 0.5 / 1 / 2 / 3, against
  P0 7.88 / 8.86 / 8.95 / 8.98 and GM 7.93 / 8.65 / 8.88 / 8.95: +7% at t = 0.5, converged by t ≈ 2. Since a ±3% shift
  of the lattice-stage m_s² moves M2-X's e_crit by ×5 (e2e_chain), lf_chain must use this column, not P0's.
- ⟨R⟩_cyc = 0.047 / 0.015 at t = 1.5 / 3 (P0 0.105 / 0.032), R_rms max 1.9 (P0 3.5). Energy budget |Δ| ≤ 0.6%.

**Self-consistency on its own background [C]** (`results/cmp_KFmid.json`, `check_KFmid.json`). The fermions were
rerun on 8 χ + 4 σ lines of the LF background itself. There the LF closure reproduces the fermion production to
Q/LF = 0.86–0.95 for t ≥ 0.5 (all species), and 1.05–1.45 at t = 0.3–0.4 (N2 worst: its κ_hom is low by ≈ 1.3–1.5
on this background). On this background Q/P0 = 2.0–3.3 and λ_out = −2.4 … −4.5 again. So the calibration is good to
≈ 10% after fragmentation and ≈ 30% (N2) before it; it was not iterated further (budget).
Caveat [A]: the 8.4×10¹⁰ LF run does not fragment at all (the early κ_hom ≈ 2–4 drain empties χ first, as in TR),
so its κ_frag after 0.75 is applied to a state it was not calibrated on; it was not checked with fermions.

## 5. Exports

- `results/hist_HFa_LF.json`: e2e_lattice hist format (lattice block of box means on the 0.004 / 0.02 grid to t = 3,
  with per-species f_N added), primary KFmid with its HB.stage2 fluid part (κ = 1, only as e2e_lattice did; do not use
  its K), alternates KFlo, KFhi, second seed, 8.4×10¹⁰.
- `results/fluid_moments_HFa_LF.json`: hand-over moments at t_sw = 1.5, 2, 2.5, 3 for every LF run: a, H, f_σ, f_χ,
  f_grad, f_s, f_rad, f_SM, per-species E_N, n_N (so ⟨E⟩/M), B_N, ρ_SM, ρ_tot, m_s², s. These are the inputs of
  e2e_chain's `handover.py` (it reads f_σ, f_χ + f_grad, f_rad, f_SM at t_sw and splits N by species).

## 6. Convergence summary

| item | check | effect |
|---|---|---|
| fermion time step | cfl 0.5 → 0.25 | ≤ 0.4% |
| fermion grid / k_Nyq | 8192 → 16384 | +2–5% (+12% transient, N2, t = 0.5) |
| realisations | n_r 4 → 8, new seed | −8 … +2% per line; line-to-line s.e. ≈ 10% (χ), 30–40% (N1) |
| lines (sampling) | 8 χ / 4 σ lines, 3 axes | dominant statistical error; κ brackets lo/hi carry it to f_rad ±0.03 |
| lattice seed (LF drain) | seed 1 → 2 | f_rad(3) +0.07%, K est. 0% |
| lattice time step, box, spacing | **not redone**; e2e_lattice P0 [Q]: dt/2 ≤ 0.5%, box ×1.5 ≤ 1%, dx/1.5 +1.5% f_rad | – |
| ray reduction (1+1D) | not testable here [A] | unknown; the largest untested systematic, with the hot spots |

**Hot spots.** Most of the produced energy after fragmentation sits at k ≈ 1.5–4.5×10⁴ (comoving), from rare sites
with a·μ5 up to 5×10⁴. Their statistics depend on the 3D lattice spacing, which was not varied for the fermions.

## 7. Answers to the pre-registered questions (this track's part)

- **Q1.** The computed closure lies **outside** [P0, GM], beyond P0: λ_PB ≈ −2 … −4 (output-level) after
  fragmentation, ×1.5–2.5 above both closures even before it. It is **not unique**: it depends on time (it grows with
  fragmentation and saturates around λ_out ≈ −2.4 to −3.6 by t ≈ 1), on species (N1 > N2 > N3 ≈ N4), and on H_i only
  through the fragmentation time. Better description: local blocking with a capacity factor κ ≈ 2–4 before and
  ≈ 10–30 after fragmentation (between P0 and unblocked tracking, nearer the latter). Coarse-graining between local
  and box-mean cannot describe it.
- **Q2, Q3, H1.** Not evaluated here (lf_chain). Input: K ≈ 5.4 (4.9–5.7) [E] at the design scale and ≈ 6.1 at
  8.4×10¹⁰ [E], below the old P0 value of 8.8. In e2e_chain the S/Q floor drops below 3.47 for K ≲ 12, so the direction
  is towards P0's situation (M2-Q not marginal) [A, to be computed]. Use the m_s² column of hist_HFa_LF.json.
- **Q4.** §0, written before any window: Q4-a (e = e_crit), Q4-b (edge of bistability), Q4-c (release synchrony);
  a and b are branch-defining in kind.

## 8. Not done, and what would change the answer

- 3D fermions (helicity mixing by ∇⊥θ), or at least rays along lattice diagonals; a 96³ recording for the hot-spot
  spacing dependence; more lines for N1; an iterated calibration (Q on the LF background → new κ → rerun); HFb; Γ_th;
  N decays during the coherent stage in the fermion runs (they are absent on both sides of the comparison).
- Back-reaction is only through the calibrated closure; the fermions themselves never act on σ, χ.

## 9. Literature (arXiv numbers verified by web search this session)

- Sz. Borsanyi, M. Hindmarsh, *Low-cost fermions in classical field simulations*, PRD 79, 065010 (2009),
  arXiv:0809.4711 (stochastic fermion ensembles).
- P. M. Saffin, A. Tranberg, *Real-time fermions for baryogenesis simulations*, JHEP 07 (2011) 066, arXiv:1105.5546
  (low-cost fermions with axial couplings in 1+1D).
- G. Aarts, J. Smit, *Real-time dynamics with fermions on a lattice*, NPB 555 (1999) 355, arXiv:hep-ph/9812413
  (classical bosons with quantum fermions; doubler treatment).
- Quoted from earlier tracks, not re-verified here: P. Adshead, E. Sfakianakis, arXiv:1508.00891 (mf_hybrid).

## 10. Files

Code (`code/`): `lf1d.py` (line solver), `regress_homog.py` (R1), `lf3d.c` (copy of lat3d.c + mixed closure, κ,
tracers, line dumps; built into scratch), `lfcommon.py` (copy of latcommon.py, scratch redirected), `jobs.py`,
`lines_io.py`, `lfline.py`, `lfqueue.py` (2-worker queue), `closure_cmp.py`, `calib.py`, `calib_seq.py`,
`kappa_table.py`, `check_kf.py`, `conv_fermion.py`, `export.py`. Raw lines, fermion states and lattice restarts are in
scratch (`.../scratchpad/wt2/lf_fermion/`). Python ran with `-B`; no binaries or `__pycache__` in the tree.
