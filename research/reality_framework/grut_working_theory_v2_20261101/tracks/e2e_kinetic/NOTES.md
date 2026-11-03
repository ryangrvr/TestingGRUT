# e2e_kinetic: thermalisation and K for HFa / HFb

Claude, end-to-end pass ("M2"), GitHub issue #13 context. This track computes, with a kinetic network, how the energy of
HFa/HFb's oscillating σ and χ condensates reaches a Standard Model (SM) thermal bath. The prior passes only bracketed
this step. Brief: `BRIEF_E2E.md`, section e2e_kinetic.

Labels: [D] derived, [C] computed here, [Q] quoted, [E] estimate, [A] argued, [R] retrospective comparison. Every rate,
width, start state and switch time is a model choice; none follows from closure. Nothing here is a validation, an
exclusion, a new theory or a GRUT prediction.

The previous run of this track was interrupted by a container restart. This resumed run reused its code, p-tables,
coherent-stage checkpoint and finished runs. It added the natural shell width, found and fixed a factor-2 error in the
perturbative shell rate, re-ran the decisive cases and wrote the exports. Section 5 says which numbers are superseded.

## 0. Result in one paragraph

**[C on the stated rates]** At the design point (H_i = 1.30×10¹¹ GeV, m~ = 0.05 eV for all four singlets, mf_hybrid's
couplings), with the self-consistent natural width of the perturbative shell:
- **HFa: K = 16.35, T_RH = 3.06×10¹² GeV.** The SM holds 50% of the energy at t = 4420 H_i⁻¹ and 99% at 7280 H_i⁻¹.
- **HFb: K = 16.32, T_RH = 3.06×10¹² GeV**, with t₅₀ = 4434 H_i⁻¹.

The SM part thermalises among itself fast: its thermal share passes 90% at t ≈ 6 H_i⁻¹. The universe as a whole is
slow. Non-perturbative production moves at most 36% of the energy into singlets N (at t ≈ 0.3–0.45). The N decay to
LH at Γ_th ≈ 0.05–0.14 H_i, and the radiation then redshifts below the surviving χ condensate. The remaining χ decays
perturbatively through χ → N₄N₄ (2M₄ < m_χ) only. That channel is **Pauli-blocked by its own products**: the shell at
k_res fills to n ≈ 1/2, and N₄ empties it only at Γ_N₄ ≈ 0.12 H_i. Reheating therefore completes at t ~ 10⁴ H_i⁻¹.

**K = 16.3 lies outside the mf_critic / mf_chain bracket K = 2.5–10.** That bracket put no self-blocking on the late
perturbative products. The computed Γ_th is also far below the bracket's fast end (10 H_i).

## 1. Model (code: `code/`)

### Sectors
**Condensates.**
- σ (m = 600 H_i) and χ (m = 250 H_i): harmonic, w = 0, cycle-averaged energies.
- Back-reaction: the energy produced as N.
- Start: the mf_hybrid bench state at t₀ = 0.19 H_i⁻¹, the first χ turning point after the waterfall [Q mf_hybrid].
  The 0.9% that the bench's unblocked friction had already moved to radiation is returned to χ.

**Axial singlets N₁…N₄** [Q mf_hybrid params]:

| Singlet | Coupled to | M (H_i) | M (GeV) | f (H_i) |
|---|---|---|---|---|
| N₁ | σ | 250 | 3.3×10¹³ | 1200 |
| N₂ | χ | 400 | 5.2×10¹³ | 2000 |
| N₃ | χ | 200 | 2.6×10¹³ | 2000 |
| N₄ | χ | 120 | 1.6×10¹³ | 2000 |

Coupling: −(M/2) N̄ e^{iγ5φ/f} N. Only N₁ (M/m = 0.417) and N₄ (M/m = 0.48) have an open perturbative channel 2M < m.

**SM.**
- A hard, non-thermal component (ρ_h, n_h), injected by the N removal.
- A thermal component ρ_th.

### Stages
**1. Coherent stage, t₀ → t_coh = 3 H_i⁻¹** (`ek_modes.c`, COHERENT).
- Exact two-level (Bloch) evolution of every comoving mode and helicity; Pauli blocking is automatic (n ≤ 1).
- Collisions are applied at the parent's turning points, plus back-reaction and expansion [C].
- Regression against the mf_hybrid skeptic's exact multi-period runs (constant amplitude, no expansion,
  `results/check_constant_amplitude.json`): N₂, N₃ and N₄ reproduce F(n half-periods) to 3–4 digits (e.g. N₄ 0.0506
  vs 0.0506 at 16). N₁ is compared at χ time and sits at the skeptic's saturated value 0.13–0.16 [C].

**2. Dephased stage, t_coh → 10⁸** (rate equations).
- Non-perturbative part: dn/dt = R(1 − 2n) − Γn + Γ_eq f_eq, with R taken from the exact single-half-period transition
  tables (`ek_ptable.py`; Stueckelberg phases between half-periods randomised).
- Perturbative part, for A = μ5,max/m < 1 in the open channels: the golden-rule shell at k_res = (m/2)β.
- Switch-time test: t_coh = 3 → 10 changes K by +0.2% (fx_natw_tcoh10). Pre-fix t_coh = 1/3/10/30 stayed within 1% [C].

### Rates
**N removal** [Q/C]:
- **Vacuum decay N → LH** with the seesaw Yukawas of mf_hybrid's completion: h_j² = m~_j M_j / v², with v = 174.1 GeV
  and Γ_rest = h_j² M_j / 8π. Thermal-mass kinematics Φ(m_L, m_H) apply, and the decay closes when M < m_L + m_H.
- **In-medium N + SM → SM + SM**: the Besak–Bödeker leading-order total (arXiv:1202.1288, verified) in damping-rate form
  c_BB h² 𝒯²/ω, with detailed balance against the thermal bath.

**Medium scale.** 𝒯² = T_th² + (12/73) n_h²/ρ_h. The asymptotic masses are m_H² = k_H 𝒯² and m_L² = k_L 𝒯².

**SM couplings.** Taken at μ = 10¹⁵ GeV from higgs_reheat's 3-loop SM RGE (read-only import).

**SM thermalisation.** The hard SM energy relaxes into the thermal component at the Kurkela–Lu under-occupied
equilibration time:

  t_eq = (34 + 21 ln(Q/T)) (Q/T)^{1/2} / [(1 + 0.037 ln 1/λ) λ² T]

This is an effective kinetic theory with LPM-suppressed splittings: Kurkela–Lu arXiv:1405.6318 and Kurkela–Moore
arXiv:1107.5050, both verified. Harigaya–Mukaida (arXiv:1312.3097, verified) give the same parametric
(Q/T)^{1/2}/(α² T) scaling for reheating. The model takes λ = 2g₂² (SU(2), the charge of the primaries L and H) [Q/A].
Variants: SU(3) gives −2% on K, and t_eq × 10 gives +0.2% (pre-fix runs).

### Perturbative shell
The shell is a model choice, and the decisive one [A/D]. With fast N removal the shell occupation sits at the steady
state n = R/(2R + Γ), so its width matters. The physical width is the line of the damped final states.
- Each N has removal rate Γ(k_res), so the pair-energy Lorentzian has HWHM Γ, or Γ/2 in ω.
- Reference ("natw"): the Lorentzian in ω is spread over the whole kinematic range of χ → N + N*(→ LH), with
  M < ω < (m² + M² − (m_L + m_H)²)/2m.
- The width is recomputed per species and per step from the current medium (`ek_modes.c`, dshell ≤ −10⁵).
- **The rate in the reheating era is set by the N₄ drain of a nearly filled Fermi sea over that kinematic range.** It
  depends only weakly on the line shape (section 3).

## 2. Deliverables

### Γ_th(t), f_rad(t), K(t)
Reference run, HFa (`results/reference.json` → tables; full rows in `results/fluid_HFa.json`).
- Γ_th: N-removal rate, energy-weighted over all N.
- Γ_cond = −d ln(ρ_cond a³)/dt: the condensate → radiation transfer rate.

| t (H_i⁻¹) | f_rad | f_SM | f_N | f_σ | f_χ | Γ_th (H_i) | Γ_cond (H_i) | K = a²H | T_SM (GeV) |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.320 | 5×10⁻⁴ | 0.319 | 0.049 | 0.631 | 0.046 | 5.6 | 1.50 | 4.4×10¹³ |
| 0.45 | 0.364 | 0.003 | 0.361 | 0.052 | 0.584 | 0.051 | 0.16 | 1.57 | 6.1×10¹³ |
| 1.0 | 0.313 | 0.012 | 0.301 | 0.056 | 0.631 | 0.077 | −0.09 (blocked back-flow) | 1.74 | 6.6×10¹³ |
| 3.1 | 0.218 | 0.045 | 0.173 | 0.063 | 0.719 | 0.127 | 4.6×10⁻³ | 2.12 | 5.8×10¹³ |
| 10 | 0.142 | 0.086 | 0.056 | 0.070 | 0.788 | 0.127 | 2.1×10⁻³ | 2.82 | 4.0×10¹³ |
| 100 | 0.080 | 0.075 | 0.005 | 0.069 | 0.851 | 0.275 | 7.6×10⁻⁴ | 5.41 | 1.3×10¹³ |
| 305 | 0.100 | 0.097 | 0.003 | 0.001 | 0.899 | 0.240 | 1.4×10⁻⁴ | 7.50 | 7.8×10¹² |
| 1009 | 0.070 | 0.069 | 9×10⁻⁴ | 0 | 0.930 | 0.128 | 6.4×10⁻⁵ | 10.7 | 4.0×10¹² |
| 3017 | 0.257 | 0.255 | 0.002 | 0 | 0.743 | 0.121 | 2.6×10⁻⁴ | 14.6 | 3.2×10¹² |
| 5087 | 0.660 | 0.656 | 0.003 | 0 | 0.340 | 0.120 | 7.1×10⁻⁴ | 16.1 | 3.0×10¹² |
| 10⁴ | 1.000 | 0.999 | 8×10⁻⁴ | 0 | 10⁻⁸ | 0.108 | — | 16.35 | 2.3×10¹² |

### Headline numbers (HFa / HFb) [C]

| Quantity | HFa | HFb |
|---|---|---|
| K (a²H/H_i in radiation; a = 1, H = H_i at σ = σ_c) | 16.35 | 16.32 |
| T_RH (T of the SM energy when ρ_SM = ρ_tot/2) | 3.06×10¹² GeV | 3.06×10¹² GeV |
| SM reaches 50% | t = 4420 H_i⁻¹ | t = 4434 H_i⁻¹ |
| SM reaches 99% | t = 7279 H_i⁻¹ | t = 7286 H_i⁻¹ |
| SM bath thermal (thermal share ≥ 90% of SM energy) | t = 6.02 H_i⁻¹ | t = 6.02 H_i⁻¹ |
| Thermal equilibrium of the universe (SM ≥ 50% and thermal ≥ 90%) | t = 4420 H_i⁻¹ | t = 4434 H_i⁻¹ |
| T_max of the thermal part | 5.5×10¹³ GeV at t = 2.24 | 5.5×10¹³ GeV at t = 2.24 |
| T_K (radiation temperature extrapolated to a = 1) | 1.23×10¹⁵ GeV | 1.23×10¹⁵ GeV |
| SM ≥ 10⁻³ | t = 0.278 | t = 0.278 |

In seconds, 1/H_i = 5.1×10⁻³⁶ s, so t₅₀ ≈ 2.2×10⁻³² s. After t ≈ 300 the σ condensate is gone and only χ is left.

### Dependence on e2e_sm's neutrino completion
The rates are linear in h_j² = m~_j M_j/v², so the completion enters only through m~_j = Σ_i|h_ji|² v²/M_j.

**All four singlets at a common m~, H_i = 1.3×10¹¹ GeV** [C; power laws fitted over 0.015–0.15 eV]:

| m~ | K | T_RH (GeV) | t₅₀ (H_i⁻¹) |
|---|---|---|---|
| 0.015 eV | 20.53 | 2.30×10¹² | 7774 |
| 0.05 eV | 16.35 | 3.06×10¹² | 4420 |
| 0.15 eV | 12.78 | 3.92×10¹² | 2674 |

**Scaling: K ∝ m~^−0.21, T_RH ∝ m~^0.23, t₅₀ ∝ m~^−0.46.**

**Per species** [C, pre-fix one-bin runs, as relative shifts]: doubling m~₄ alone changes K by −5.4%. Doubling m~₂ or
m~₃ changes it by about −2%, and m~₁ by −0.1%. **N₄, the sink of the χ → N₄N₄ channel, is the lever.** With the
natural width the N₄ lever carries most of the m~^−0.21. If e2e_sm's completion changes m~₄ by a factor x, use
K ≈ 16.35 x^−0.2 as the first estimate, then rerun `ek_extra.py`.

### Dependence on H_i
Rescaling rule of mf_chain: masses and f fixed in H_i units, amplitude fixed in M_P units [C]:

| H_i (GeV) | K | T_RH (GeV) | t₅₀ (H_i⁻¹) |
|---|---|---|---|
| 5.0×10¹⁰ | 9.19 | 1.21×10¹² | 10594 |
| 1.3×10¹¹ | 16.35 | 3.06×10¹² | 4420 |
| 2.4×10¹¹ | 17.49 | 5.59×10¹² | 2427 |

At 5.0×10¹⁰ GeV the radiation briefly exceeds 50% at t = 0.22, before the χ redshift takes over. K is not monotonic
in the decay rates: lower H_i means larger amplitudes in H_i units and a longer non-perturbative stage.

## 3. Numerics and model sensitivity [C]

**Grid.** 6400 vs 12800 log bins in the dephased stage: K = 16.3485 vs 16.3485, t₅₀ 4420 vs 4427, T_RH 3.062 vs
3.061×10¹². The pre-fix one-bin shell was grid-dependent (K = 20.7 / 22.4 / 24.3 at 350 / 700 / 1400 bins). That was
the bin width acting as an unphysical line width, and it is what this resumed run fixed.

**Line shape and width** (all at nq ≥ 6400, pre-fix rate except where stated):

| Line shape | K |
|---|---|
| Lorentzian in ln k at fixed HWHM 0.0063, 6400 / 12800 bins | 17.613 / 17.614 |
| Self-consistent Lorentzian in ln k, cut at 40 HWHM | 17.70 |
| Same, cut at 80 HWHM | 17.37 |
| ω-Lorentzian over the kinematic range | 17.20 |
| Width × 0.5 | 18.87 |
| Width × 2 | 16.62 |
| After the fix: width × 2 | 15.61 (reference 16.35) |

So K ∝ width^−0.07 to −0.09, with a 3% spread from the tail treatment.

**Moving or fragmented parent** (lattice input): the products are spread as a Gaussian in ln k. σ = 0.1 gives
K = 15.52 (T_RH 3.38×10¹²); σ = 0.3 gives K = 13.47 (T_RH 4.17×10¹²).

**SM rates** (pre-fix, relative):
- c_BB × 2 / × 0.5: −2.6% / +1.0%.
- Kurkela–Lu t_eq × 10: +0.2% on K. The bath-thermal time moves from 6 to 770 H_i⁻¹.
- λ = 3g₃²: −2%.

**Coherent-stage numerics** (pre-fix): Δk 0.05 → 0.035 m and collision step 0.04 → 0.02 change K by < 10⁻⁴.

**Energy.** Bookkeeping is closed by construction (`results/check_energy_conservation.json`, no expansion, 5% initial
bath). Re-binning to the log grid conserves number to 10⁻⁸ and energy to 5×10⁻⁶.

**Lattice-interface test.** The run was stopped at t_sw = 1.0 and restarted from **moments only** (per-species N
energy and number, condensates, SM). K = 22.34 vs 22.43 uninterrupted (0.4%) and t₅₀ to 0.3% (pre-fix shell; the
restart affects only the early stage).

## 4. Interface for the chain and the lattice (`results/fluid_HFa.json`, `results/fluid_HFb.json`)

**Format.**
- The mf_hybrid bench history columns (t, N = ln a, a, H, ε_H, R, R_osc_σ, R_osc_χ, m_s, f_σ, f_χ, f_rad), with
  cycle-mean R and w = f_rad/3.
- Kinetic columns: f_SM, f_N, f_SM_hard, f_SM_thermal, T_SM, T_thermal, 𝒯², m_H²(medium) = k_H 𝒯², τ_KL, K_running,
  Γ_th, Γ_cond, and per species f_Nj, Γ_band, n̄_band, A_j.
- 689 rows, log-spaced in t, from t₀ = 0.19 to 1.07×10⁵ H_i⁻¹. Column definitions are in `history_columns`.
- For the bench (t < 0.19) use mf_hybrid's `bench_<name>.json`.
- Spectator mass m_s is in mf_chain's convention.

**Switch to a lattice coherent-stage history.** The lattice (e2e_lattice) runs from the waterfall to t_sw. Take
t_sw ∈ [0.45, 3] H_i⁻¹: after the first χ half-periods, while the kinetic code's coherent stage is still running.
Two ways:
- **(a) Moments hand-over (preferred).** Write a JSON with t, a, H; ρ_σ and ρ_χ (zero-mode, cycle-mean); ρ_h, n_h
  (SM not yet thermal; zero if the lattice has no SM); ρ_th; and lists E_N[4], n_N[4] (per-species N energy and number
  densities from the lattice's Pauli drain). All are physical, in H_i units. Then run:

  `python3 -B code/ek_run.py export HFa TAG --state lattice.json --out <path>`

  The defaults are the reference shell (dshell −1e5, 6400 bins). The occupations are rebuilt as uniform up to
  k_max = 4E/3n with n̄ ≤ 1/2, and the coherent stage continues to max(3, t_sw + 2). Tested at t_sw = 1 (0.4%).
- **(b) Splice.** Use the lattice to t_sw and these rows at t ≥ t_sw. This holds only if the lattice's f_σ, f_χ and
  f_rad at t_sw agree with the rows to a few %. Otherwise use (a).

**What the lattice must also supply.** If the condensate fragments, give the fraction of χ energy in non-zero modes
and the ln k spread σ of the χ quanta. Rerun with `--dshell σ`: σ = 0.1 → K 15.5; σ = 0.3 → K 13.5. **This is the
largest open input to K**, larger than the line shape.

## 5. What is superseded

- `scan.json` entries without the `fx_` prefix were made with the perturbative shell rate a factor 2 too low. Each
  decay fills the pair (k, −k), so the per-entry rate is g_i N_tot/(w_i/a³), not N_tot/(2w_i/a³). The fix is in
  `ek_modes.c` at line ~434.
- The one-bin entries (dshell = 0, among them ref, H*, mnu*, cBB*, KL*, lamSU3, tcoh*, restart*) are also
  grid-dependent.
- The fluid files of the interrupted run (one-bin, K ≈ 22.4) are overwritten.
- The pre-fix entries are kept only for relative sensitivities (SM rates, per-species m~, t_coh, the restart test),
  and are used that way in section 3.

## 6. Limitations [A]

- **Homogeneous parents.** Fragmentation is the lattice's job (section 4).
- **The off-shell tail is treated as pair production.** It uses Pauli factor 1 − 2n and 2N per event. Physically the
  far tail is χ → N + N*(→ LH), with factor 1 − n and only one N. This bears on the 3% tail spread and may lower K
  somewhat. Not computed.
- **No elastic N–SM momentum diffusion out of the shell.** It is O(h²) like the in-medium damping, and small here:
  Γ_in-medium ≪ Γ_rest at T ≲ 50 H_i.
- **SM cascade.** One relaxation rate (Kurkela–Lu, under-occupied) for the whole hard SM component. No flavour or
  charge structure, no soft-sector evolution. The KL fit is for pure SU(N) gluodynamics.
- **G12 = 4.3×10⁻⁴ is extrapolated** from Besak–Bödeker's table (±10%).
- **The SU(2) λ at μ = 10¹⁵ GeV** is a choice of renormalisation scale.
- **The N's leptogenesis CP asymmetry is not tracked.** That belongs to e2e_sm: T_RH = 3×10¹² GeV < M_N, so the N
  population is non-thermal.
- **The spectator s** is carried only through m_s(t). Its own decay (Γ_s = 1.3×10⁻⁶ H_i assumed, bench) is not part of
  this network.

## 7. Files

**`code/`** (python3 -B; the binary compiles into scratch):
- `ek_common.py`: constants, rates, I/O.
- `ek_modes.c`: the engine. Coherent, dephased and table modes, with the natural-width shell and the rate fix added in
  the resumed run.
- `ek_ptable.py`: transition tables.
- `ek_run.py`: driver, regression checks and fluid export.
- `ek_scan.py`: pre-fix scans.
- `ek_extra.py`: this run's cases; `fx_*` are authoritative.
- `ek_post.py`: Γ_cond column, scan.json, reference.json.

**`results/`:** fluid_HFa.json, fluid_HFb.json, reference.json (headline, tables, scalings), scan.json (all runs),
check_constant_amplitude.json, check_energy_conservation.json, ptable_summary.json.

**Reused read-only.**
- mf_hybrid `results/bench_HFa.json` and `bench_HFb.json` (parameters, start state).
- The mf_hybrid skeptic's `sk_axial.json`: the regression target for the coherent stage. It holds the exact
  multi-period output of the axial and Pauli codes of mf_hybrid and its skeptic. This track's engine is its own C code
  and does not import `axial.py` or `sk_pauli.py`.
- higgs_reheat's `sm_rge.py`, imported for the couplings.

## 8. Literature

Verified in this resumed run by search of the arXiv listing:
- Kurkela & Lu, arXiv:1405.6318 (PRL 113, 182301).
- Kurkela & Moore, arXiv:1107.5050.
- Harigaya & Mukaida, arXiv:1312.3097.
- Besak & Bödeker, arXiv:1202.1288 (JCAP 1203:029).

Cited through other tracks and not re-verified here:
- Buttazzo et al., arXiv:1307.3536 (SM RGE matching, via higgs_reheat).
- Arnold, Moore & Yaffe, hep-ph/0209353 (the effective kinetic theory under Kurkela–Lu).
