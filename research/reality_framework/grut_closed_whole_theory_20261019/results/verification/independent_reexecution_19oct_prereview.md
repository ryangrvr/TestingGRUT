# Fresh re-execution of the closed-whole calculations (workflow 2, rerun step)

Date: 2026-10-19. Repo: /home/user/TestingGRUT, branch codex/grut-closure-fork-20261019, HEAD 96a7e7e
(the packet directory research/reality_framework/grut_closed_whole_theory_20261019/ is untracked at this commit).
Nothing under /home/user/TestingGRUT was modified. All files written under
/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/ (abbreviated `$W` below).

Reference inputs (read-only):

- `$R` = /home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019
- scripts: `$R/code/closed_whole_demo.py` (sha256 ba50f60e58a0219a228ea5829b2116c4e1ca3f2c7008e4ff2eb8f51de7090ed2),
  `$R/code/marginal_compatibility.py` (sha256 1aad4b064ec518d7928ec454c2941267a9b19df84e1dab1d5f6731b85beb2cec)
- reference outputs: `$R/results/closed_whole_demo.json` (100188 bytes, sha256 f7eab04b3d606561949fa2d7280764d315c80fdb22d6f5455561220bd15e47d5),
  `$R/results/marginal_compatibility.json` (3250 bytes, sha256 1dbc40df2dfea69a5691862f9bdb12a645862df42bbfe03423eadc4364feaac8)

Environment (identical to the reference provenance strings): Python 3.13.16 (main, Oct 1 2026) GCC 13.3.0,
numpy 2.5.3, scipy 1.18.1, cvxpy 1.9.3, CLARABEL available, OpenBLAS 0.3.34 (scipy-openblas64, Haswell kernels),
Linux-6.18.44-fc-v139-x86_64-with-glibc2.39, 4 CPUs. `/usr/bin/time` is not installed; process wall time was
measured with `date +%s.%N` around the python3 call.

## Verdict

1. Every numeric value agrees: **maximum relative deviation 0.0** over all numeric leaves of both JSON files
   (5030 numeric leaves in closed_whole_demo.json, 24 in marginal_compatibility.json; every one bit-exact equal).
2. **No key is present in one file and not the other** (5058/5058 common leaves for the demo, 46/46 for the marginal file).
3. **All SDP/LP statuses match**: C2 `infeasible`, C4 chain `optimal`, C4 triangle `optimal`, C5 odd cycle `infeasible`,
   C6 `infeasible`. (closed_whole_demo.json contains no status keys; it is numpy/scipy exact diagonalization only.)
4. The two reruns are **bitwise identical to each other and to the reference once `provenance.wall_seconds` is removed**;
   with nothing removed, `provenance.wall_seconds` is the only differing line (the marginal file's run 1 even reproduced
   the reference byte-for-byte, wall 0.29 s in both).
5. Wall time: closed_whole_demo.py 50.75 s / 50.05 s process wall (in-script 50.22 s / 49.51 s; reference 53.03 s);
   marginal_compatibility.py 1.55 s / 1.44 s process wall (in-script 0.29 s / 0.26 s; reference 0.29 s).

## Commands

```bash
W=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2
R=/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019

# layout: scripts in rerun/code, outputs land in rerun/results via the scripts' own HERE.parent/"results"
mkdir -p $W/rerun/code $W/rerun/results $W/run1 $W/run2
cp $R/code/closed_whole_demo.py $R/code/marginal_compatibility.py $W/rerun/code/
sha256sum $W/rerun/code/*.py      # ba50f60e... and 1aad4b06... : identical to the originals

# run 1
cd $W/rerun/code
S=$(date +%s.%N); python3 closed_whole_demo.py       > $W/run1/closed_whole_demo.stdout;       E=$(date +%s.%N)   # 50.75 s
S=$(date +%s.%N); python3 marginal_compatibility.py  > $W/run1/marginal_compatibility.stdout;  E=$(date +%s.%N)   # 1.55 s
cp $W/rerun/results/*.json $W/run1/

# run 2 (same commands) -> $W/run2/                                                                # 50.05 s, 1.44 s

# key-by-key comparison (flattens to dotted paths; rel dev = |a-b|/max(|a|,|b|))
python3 -I $W/compare_json.py $R/results/closed_whole_demo.json      $W/run1/closed_whole_demo.json      --exclude-prefix provenance
python3 -I $W/compare_json.py $R/results/marginal_compatibility.json $W/run1/marginal_compatibility.json --exclude-prefix provenance
# same with only wall time excluded (so provenance strings are checked too)
python3 -I $W/compare_json.py $R/results/closed_whole_demo.json $W/run1/closed_whole_demo.json --exclude-key provenance.wall_seconds

# bitwise check: re-serialize with provenance.wall_seconds popped, then cmp / sha256
for run in run1 run2; do for f in closed_whole_demo marginal_compatibility; do
  python3 -I $W/strip_wall.py $W/$run/$f.json $W/stripped/${run}_$f.nowall.json; done; done
for f in closed_whole_demo marginal_compatibility; do python3 -I $W/strip_wall.py $R/results/$f.json $W/stripped/ref_$f.nowall.json; done
sha256sum $W/stripped/*.nowall.json
cmp $W/stripped/run1_closed_whole_demo.nowall.json      $W/stripped/run2_closed_whole_demo.nowall.json
cmp $W/stripped/run1_marginal_compatibility.nowall.json $W/stripped/run2_marginal_compatibility.nowall.json

# raw diffs (no stripping)
diff $R/results/closed_whole_demo.json $W/run1/closed_whole_demo.json
diff $W/run1/closed_whole_demo.json    $W/run2/closed_whole_demo.json
diff $R/results/marginal_compatibility.json $W/run1/marginal_compatibility.json
diff $W/run1/marginal_compatibility.json    $W/run2/marginal_compatibility.json
```

## Numbers

### Key-by-key comparison, reference vs run 1 (compare_json.py, `--exclude-prefix provenance`)

| file | leaves ref / new / common | keys only in ref | keys only in new | numeric leaves | bit-exact | max rel dev | max abs dev | non-numeric leaves | non-numeric diffs |
|---|---|---|---|---|---|---|---|---|---|
| closed_whole_demo.json | 5058 / 5058 / 5058 | 0 | 0 | 5030 | 5030 | 0.000000e+00 | 0.000000e+00 | 22 | 0 |
| marginal_compatibility.json | 46 / 46 / 46 | 0 | 0 | 24 | 24 | 0.000000e+00 | 0.000000e+00 | 14 | 0 |

With only `provenance.wall_seconds` excluded (provenance strings included): demo 27 non-numeric leaves, 0 differing;
`script_sha256`, `python`, `numpy`, `scipy`, (`cvxpy`, `solver`), `platform` all SAME in both files.
Exit status of compare_json.py: 0 in every case.

### SDP / LP statuses (marginal_compatibility.json), reference = run 1 = run 2

| key | reference | rerun |
|---|---|---|
| C2_bell_chain.sdp_status_common_global_state | infeasible | infeasible |
| C4_werner_thresholds.chain_status | optimal | optimal |
| C4_werner_thresholds.triangle_status | optimal | optimal |
| C5_classical_contrast.odd_cycle_anticorrelated_status | infeasible | infeasible |
| C6_inconsistent_overlap.sdp_status | infeasible | infeasible |

Solver values, identical to 12 significant digits in all three files: C3 max F_AB+F_BC = 1.49999999899,
max symmetric min fidelity = 0.750000000106, optimum pair (0.749999999496, 0.749999999496); C4 chain p* = 0.666666673771,
triangle p* = 0.33333333563; C1 nested dev 0.0 and 3.46944695195e-18; C5 chain extension error 5.55111512313e-17.

### Bitwise identity of the two reruns (and of the reference), `provenance.wall_seconds` removed

| file | sha256 (wall stripped) ref | run 1 | run 2 | cmp run1 vs run2 |
|---|---|---|---|---|
| closed_whole_demo.json | 17e1fd6da29a02dea10bbf215b86eb1612fbe10934aaa9927aca32e27a3a654f | same | same | identical |
| marginal_compatibility.json | 1ef9002357fdae754b42557ea4eda121f51b31f5939fd15feedb0c92cb4c3e55 | same | same | identical |

Raw files (nothing removed): all 100188 / 3250 bytes. `diff` output is exactly one line per pair, the wall-time line:

```
closed_whole_demo.json      ref vs run1:  5210c5210  "wall_seconds": 53.03 -> 50.22
closed_whole_demo.json      run1 vs run2: 5210c5210  "wall_seconds": 50.22 -> 49.51
marginal_compatibility.json ref vs run1:  (no diff; sha256 1dbc40df... identical, wall_seconds 0.29 in both)
marginal_compatibility.json run1 vs run2: 78c78      "wall_seconds": 0.29 -> 0.26
```

Raw sha256: run1 demo 76505d12bce5d5e17cfb229dcc665cb1be891f63d84fe678138d47ae82da03b0, run2 demo a2025e35726e6b443db62e8b5f4005d0852fd810f73b0edea6d9a003f01a7457,
run1 marginal 1dbc40df2dfea69a5691862f9bdb12a645862df42bbfe03423eadc4364feaac8 (== reference), run2 marginal f9726958c5e9d37248fb9f6c2c6002ad8417ca560e1e804fbc667fa384558e83.

stdout: run 1 vs run 2 differ only in the demo's final `wall 50.2s` / `wall 49.5s` line; marginal stdout identical.
stderr empty (0 bytes) for all four runs.

### Wall time

| script | run 1 process wall | run 2 process wall | run 1 in-script | run 2 in-script | reference in-script |
|---|---|---|---|---|---|
| closed_whole_demo.py | 50.75 s | 50.05 s | 50.22 s | 49.51 s | 53.03 s |
| marginal_compatibility.py | 1.55 s | 1.44 s | 0.29 s | 0.26 s | 0.29 s |

(process wall includes interpreter + numpy/scipy/cvxpy import; in-script = `provenance.wall_seconds`, measured from `t0` after
model setup in the demo and from the start of the computations in the marginal script.)

### Run-1 stdout of closed_whole_demo.py (identical in run 2 except the wall line)

```
S1 invariants: max|dE_tot| = 7.11e-15, max|dS_tot| = 4.00e-15, max|d purity| = 6.66e-16; Born local-vs-global max dev = 2.2e-16
S1 ranges: E_A in [-0.1693, 0.0000], E_int in [-0.2229, 0.0000], S_A in [0.0000, 0.4893] (ln2 = 0.6931)
S1 reduced-equation check at t=3: |lhs-rhs| = 1.9e-10 (full), 1.9e-10 (without H_Abar term); |Tr_Abar[H_Abar,rho]| = 3.0e-16
S2 S_A: t0 0.0000, max 0.4893 at t=37.5, min after max 0.4042 at t=39.5, sign changes of dS/dt: 26; BLP backflow sum 1.2035 (D: 1.000 -> min 0.542)
S3 CPTP: min Choi eig 2.96e-11, TP dev 2.7e-15; CP-divisibility most negative intermediate Choi eig -0.1590
S3 TCL generator |K(t)-K(0.5)| thermal: {'0.5': 0.0, '1.0': 0.1598, '2.0': 0.5202, '5.0': 0.9689, '10.0': 0.6789, '20.0': 0.601}
S3 TCL generator |K_thermal-K_ground|: {'0.5': 0.0172, '1.0': 0.0758, '2.0': 0.1001, '5.0': 0.1084, '10.0': 0.0688, '20.0': 0.0651}
S4 inf_T_bath_quantum_correlation_Y0X1: eps=0.001953 (max 0.003906), marginal dev 0.0e+00, D(rhoA) t0 1.6e-18, max 0.1025 at t=0.7, |inhom(0)|=0.3536
S4 inf_T_bath_classical_correlation_Z0Z1: eps=0.001953 (max 0.003906), marginal dev 0.0e+00, D(rhoA) t0 4.4e-16, max 0.0097 at t=14.1, |inhom(0)|=0.0000
S4 thermal_beta1_bath_dressed_quantum_correlation: eps=0.25 (max 0.5), marginal dev 2.2e-16, D(rhoA) t0 6.9e-18, max 0.0544 at t=0.7, |inhom(0)|=0.1925
S5 b' exact recurrence: D(T)=7.6e-15, S_A(T/2)=0.4487, S_A(T)=0.0e+00
S5 trace distances rhoA at t=5: {'a_vs_b': 0.6046, 'a_vs_c': 0.2765, 'b_vs_c': 0.7107} | model b recurrence min S_A after max = 7.56e-02 at t=31.4
S6 Born p(Z+,t1)=0.354910 register pop=0.354910; p(Z+,t2) with record 0.338531, without 0.449147, dephased 0.338531
S6 record drift at t=10: kappa=0 7.1e-15, kappa=0.3 0.2505
S7 forward-backward trace distance 5.1e-15; reversed preparation S_A: {'0.0': 0.4198, '5.0': 0.4005, '10.0': 0.4484, '20.0': 0.3337, '30.0': 0.3929, '40.0': 0.0}
S8 pure-global checks: |S_A-S_Abar| max 3.3e-16, subadditivity violations 0, Araki-Lieb violations 0, closest return D=0.4012 at t=11.1
wall 50.2s; sha256 ba50f60e58a0219a
```

## Caveats (scope of this check)

- This is a reproduction check only: it establishes that the recorded JSON files are what the recorded scripts produce in
  this environment, deterministically (two runs bit-identical). It does not audit the physics or the correctness of the
  scripts' claims, and it claims no new physics.
- Bit-exactness was obtained on the same machine image, same numpy/scipy/OpenBLAS build and same CPU kernel path as the
  reference. On a different BLAS/CPU the floating-point leaves (already rounded to 12 significant digits by the scripts'
  `rnd`) could differ in the last digits; the SDP values (about 1e-9 from the analytic 3/2, 3/4, 2/3, 1/3) are solver
  tolerances, not exact numbers.
- `compare_json.py` treats a JSON leaf as numeric if it is int or float (bool excluded) and uses rel dev = |a-b|/max(|a|,|b|)
  with 0/0 := 0.

## Files written

- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/rerun/code/closed_whole_demo.py and marginal_compatibility.py (verbatim copies)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/rerun/results/closed_whole_demo.json and marginal_compatibility.json (run 2 outputs, last written)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/run1/ and run2/ (per-run JSON, stdout, stderr, .time)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/stripped/ (wall-stripped JSONs used for cmp)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/compare_json.py and strip_wall.py (helpers)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/run1/compare_demo_vs_ref.txt and compare_marg_vs_ref.txt (full comparison logs)
- this file: /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/rerun/REPORT.md
