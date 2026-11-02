# mf_ideal code: the M1r chain on an external (idealised or benchmark) history with a spectator holder

Claude, multi-field inflaton pass (issue #13 context). Same-session AI work on one supplied model. Every history, ramp,
width and slice convention is a model choice, not a consequence of closure.

## What is here

| file | role |
|---|---|
| `rpmfi.c` | the engine: copy of `tracks/inf_chain/code/rpchain.c` (lineage rpinf.c, rpdq.c, rptri.c) plus **MODE 5** (header, marks MF0–MF2). MODE 0–4 are byte-identical to rpchain.c / rpdq.c (checked by `regress.py`) |
| `vmodes.c` | longitudinal-vector mode integrator, copied unchanged from `tracks/inf_chain/code/vmodes.c` |
| `mfi_common.py` | build, driver, parser, preparation (dq_heat skeptic's `dq_model_copy.py`, read-only), fluid K(Γ) |
| `release_mfi.py` | one history: dense run (τ ≤ 100) + long run (τ ≤ 10⁴, conformal bins) + a⁻² freeze-out fits → `results/release/` |
| `vectors_mfi.py` | vectors on the history's own recorded background → `results/vectors/` |
| `chain_mfi.py` | ledger, S/Q(e), e_crit, branch, relic, H_i, isocurvature, Hartree check → `results/chain/` |
| `grid.py` | the idealised grids (P primary, D small decay, F fine K = 1, X ξ_s tracking, I isocurvature) |
| `osc_bg.py`, `smooth_end.py` | sensitivity backgrounds (oscillating matter of mass m_bg; finite-duration end) |
| `bench_run.py` | **runs the whole chain on a `bench_*.json` file** (below) |
| `regress.py`, `summary.py` | regression (R1 MODE 0 byte check, R2 table bridge); the map `results/sq_map.json` and `results/summary.json` |

Binaries and run series live in the scratch area
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_ideal/` (`build/`, `runs/`); set
`MFI_SCRATCH` to move it. Run Python with `-B` and one thread (`OMP_NUM_THREADS=1`).

## MODE 5 in one paragraph

The hidden sector (M1r at λ_Y = λ_Σ = 0.5) is integrated on an external background (a, H, R). The holder is
u = μ + G_eff c with c = s/s_*. The engine's `gre gim` arguments carry G_eff = u_slice − μ (= 10i − μ for M1r's
calibration). The preparation is the held stationary point at u_slice and R = 12 (constant-H inflation, s frozen).
The background is either the built-in **fluid** (dust decaying into radiation at the `gamma` argument, R = ρ_m) or a
**table** (`MFI_TABLE`). The spectator either **evolves**, with c'' + (3H + Γ_s)c' + m_s²(t)c + (H_i/s_*)² ∂V_h/∂c = 0,
or **follows the table** (`MFI_STAB=1`). The environment variables are listed in the header of `rpmfi.c`.

Table format (text): the first line is `NROWS NCOLS` (4, 6 or 7). Each row is `t lna H R [c cdot [ms2]]`, in H_i units,
with a = 1 and H = 1 at t = 0. ln a and c are cubic-Hermite with their exact slopes (H, cdot). H, R and cdot use
finite-difference slopes. Beyond the last row the engine continues with exact radiation (a²H fixed, R = 0). So **the table
must end in radiation** (the bench files run to K convergence).

## Running on a bench_*.json file (mf_chain)

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_ideal/code; export OMP_NUM_THREADS=1
python3 -B bench_run.py ../../mf_hybrid/results/bench_HFa.json --tag HFa            # ≈ 1.5 min, one history
python3 -B bench_run.py ../../mf_hybrid/results/bench_HFa.json --tag HFa_mean --R trace --win 0.05   # cycle-mean R check
```

**Fields read from the bench file:**
- `history.t`: H_i⁻¹, with the slice at t = 0. Rows with t < 0 are used only for the exit-factor check.
- `history.N` (ln a), `H`, `R`, `m_s` (H_i), `s_over_sstar`, `s_is_envelope`, `R_osc_sig` (null marks the homogeneous part).
- `params.H_i_GeV`, `params.s_star_GeV`.
- `spectator.Gamma_s_assumed`, if present.

The script inserts a t = 0 row if needed, renormalises a = 1 at t = 0 and writes the MODE 5 table (7 columns). It then
runs release → vectors → chain.

**Options and what they mean:**
- `--R raw` (default): use the bench's R. This is valid if the 2m_χ, 2m_σ oscillation of R is resolved by the sampling;
  bench_HF* sample at Δt = 4×10⁻⁴ there. The script prints the median Δt of the homogeneous part.
- `--R trace`: replace R in the homogeneous part by the cycle mean 6H̄²(2 − ε̄), with H̄ a running mean over `--win`. The
  fluid part is already a cycle mean.
- `--s evolve` (default): the engine evolves s with m_s²(t) from the bench's column (`MFI_SM2TAB=1`), the back-reaction
  and Γ_s. The bench's exact s(t) is compared (`s_check` in `results/bench/bench_run_<tag>.json`).
- `--s table`: s follows `s_over_sstar` inside the table. The table is then cut where `s_is_envelope` starts. This is valid
  only if the background has reached radiation by then.
- `--gs` overrides Γ_s. `--tlong` sets the freeze-out run length (default 10⁴).

**Outputs:**
- `results/release/release_<tag>.json`: N_X, E_X, ε, E/E_c, N_YL, N_YH, N_σ, K, Gauss/work/H checks.
- `results/vectors/vectors_<tag>.json`: N_A(e), R_vec(e).
- `results/chain/chain_<tag>.json`: S/Q(e), the branch at 10⁻⁶ and 3×10⁻⁶, e_crit (bulk, Q = 2×10⁵), the relic (M, p),
  H_i,req and r.
- `results/bench/bench_run_<tag>.json`: the interface diagnostics, including the pre-slice exit factor
  f = [h_r F_s/F(R = 6(2 − ε)h_r²)]² at 1–4 e-folds.

**What bench_run.py does not do (mf_chain's items):**
- the self-consistent H_i iteration (rescale V0, re-export, re-run until H_i,req = H_i);
- hidden heat from s and σ decays;
- gravitational production at the abrupt end;
- heavy remnants beyond the homogeneous engine;
- Q-matter ⟨Q⟩ and M beyond the quoted qm_micro table.

The isocurvature derivative needs two extra releases at c0 = 1 ± 0.01 (see `chain_mfi.isocurvature`). For a bench,
call `release_mfi.release(tag, K, m_s, 0, bg='table:<file>', c0=...)` with `MFI_SM2TAB` set, as `bench_run.py` does.

## Reproduce this track

```bash
python3 -B regress.py                                          # R1 + R2 (~30 s)
python3 -B grid.py P 0 2 & python3 -B grid.py P 1 2 & wait      # primary map (~8 min)
#   (D, I, F, X: same with the set name; grid.py skips finished tags)
for m in 3 6 10 30; do python3 -B osc_bg.py $m 2.0 0.1 1; done; python3 -B osc_bg.py 10 1.3 0.1 1; python3 -B osc_bg.py 30 1.3 0.1 1
for e in 0.03 0.1 0.3 1.0; do python3 -B smooth_end.py $e; done
python3 -B vectors_mfi.py <tags>; python3 -B chain_mfi.py P D F X I <tags>
python3 -B summary.py
```

The isocurvature (`chain_mfi.isocurvature`; plus `extras.py iso05` and `extras.py isoupdate`), Hartree
(`extras.py hartree`), holder-amplitude (`extras.py holder`), s_* / mirror / fast-decay (`extras.py sstar`) and
tolerance (`extras.py rtol 1.03`) checks are in `extras.py`. Their results are recorded in `../NOTES.md` §6.
