# Track `independent_solver`: an independent radial solver for the prepared formation run

Claude review of Codex packets ff49e5d (causal receiver) and bb211d1 (core conversion), 25 October 2026.

**Independence.** This track wrote its own solver from the equations in Codex's MATHEMATICS.md files (c4b99f3
section 4 and bb211d1 section 2). It never ran Codex code.

**Sources of numbers.**
- Codex's numbers are quoted [Q] from the packet JSON files, read through `git archive` into the scratchpad.
- One configuration item was read from Codex code: the old cubic damper profile in bb211d1
  `radial_nodal.cpp`, gamma = 0.03 ((r-270)/90)^3. It was used only to diagnose the old boundary (section 3.2).
  It is labelled [Q-code].

**Claim labels.** [D] derived, [C] computed here, [Q] quoted, [E] estimate, [A] argued.

Everything here concerns one classical spherical prototype with supplied fields. Nothing here is evidence
that the fields exist, a cosmological population, or a quantum particle count.

## 1. What was solved

The equations, from MATHEMATICS.md [Q], with a = 0.1 and b = 3:

    psi_tt = lap psi - (S^2 + a|psi|^2) psi,      S_tt = lap S - [b(S^2-1) + |psi|^2] S,
    Qbar = -int Im(psi* psi_t) d^3rho.

The energy Ebar is defined as in c4b99f3. The variables are U = r psi (complex, mass 1) and Z = r(S-1)
(real, mass sqrt 6).

The prepared data are:
- psi = 0.18 exp(-rho^2/7200) and psi_t = -i psi;
- S from the static neutral equation, with S_t = 0.

My Ebar0 = 38951.47049696 and Qbar0 = 38969.38665303. Codex's are 38951.470496956 and 38969.386653027 [Q].
They agree to 11 digits [C].

## 2. Method, chosen to differ from Codex

| | Codex ff49e5d | This track |
|---|---|---|
| Space | 4th-order staggered finite differences | **Fourier-sine (DST-I) pseudo-spectral**, h = 1, 0.5, 0.25 (h = 0.125 only to tau = 1956) |
| Time | RK4, dt = 0.015625 | **Symplectic splitting**: exact free Klein-Gordon rotation of each sine mode, a nonlinear position kick, Yoshida 4th-order composition; dt = 0.0625, 0.03125, 0.015625 |
| Boundary | CFS-PML at r = 360 to 450 | **(a) Box too big for reflections to return**: hard wall at L = 6400, so nothing leaving r < 600 can return to r = 60 before tau = 12140. **(b) Own absorber**: smooth velocity sponge gamma = 0.03 ((r-R_s)/3000)^5, applied exactly in the kick with an energy/charge receiver; R_s = 1200 for the long runs |
| Static S | continuum solve_bvp | Spectral fixed point on the grid. It agrees with an independent solve_bvp in Z to 7e-14 |
| Q60 | nodal sum | Trapezoid plus Euler-Maclaurin end term (h^2/12)q'(60). Without this term the observable converges only as h^2 |

Conservation [C]:
- Charge is conserved to round-off: |dQ/Q0| <= 2e-10 with the sponge receiver and <= 1.5e-12 in the box.
- Maximum relative energy error: 7.7e-6 at dt = 0.0625, 4.8e-7 at dt = 0.03125, 3.0e-8 at dt = 0.015625.

## 3. Results

### 3.1 Reflection at p/m_X = 0.02, 0.1, 0.5 (energy fraction)

| | 0.02 charged | 0.1 charged | 0.5 charged | 0.02 neutral | 0.1 neutral | 0.5 neutral |
|---|---|---|---|---|---|---|
| Big box (production) [C] | 0 before causal return | 0 | 0 | 0 | 0 | 0 |
| Own sponge, continuum ODE [C] | 5.8e-13 | 1.4e-24 | 1.1e-28 | 6.4e-11 | 8.1e-24 | 8.7e-30 |
| Own sponge, time-domain packets [C] | 1.4e-11 | 6.1e-14 | 4.0e-25 | 8.5e-10 | 5.1e-14 | 5.3e-24 |
| Codex CFS lattice [Q] | 3.3e-16 | <=1.8e-18 | <=4.0e-17 | 4.2e-7 | <=3.0e-19 | <=8.1e-18 |

The packet errors are measured in the physical region, relative to the same packet propagated exactly in a box
twice as large. Field energy plus receiver closes to <= 2e-10.

For the box, the measured error before the return time is <= 6.4e-25. The sponge changes Q60 relative to the box
by at most 1.8e-9 (absolute fraction) through tau = 12000.

Two numbers were rederived exactly (CONFIRMED):
- Codex's ideal-continuum CFS wall return. Analytically, |R|^2 = exp(-4k * 450 * w/(0.01+w^2)), giving
  1.8217e-8 (charged, p = 0.01), 4.2465e-7 (neutral, p = 0.02) and 6.515e-4 (neutral, p = 0.01).
- 9636950's old-damper table, from my own ODE code: 31.73% at p = 0.02, 4.92% at p = 0.5, 14.83% at p = 1.

### 3.2 Retained core charge, Q60 mean over tau in (9600, 12000]

| Run | Q60 mean |
|---|---|
| h = 0.5, dt = 0.0625 / 0.03125 / 0.015625 [C] | 35.142 / 35.346 / **35.357%** |
| Fourth-order Richardson in dt | 35.357% (the differences 0.2035 and 0.0111 give a ratio of 18) |
| h = 0.25, dt = 0.0625 [C] | 35.144% (h = 0.5 to 0.25 changes it by +0.002 pp) |
| h = 1 [C] | 34.58% (under-resolved after tau of about 840) |
| Codex P00625 / P0125 / P025_time [Q] | 35.363 / 35.357 / 35.188% |
| bb211d1 N00625 [Q] | 35.153% |
| Same mean with Codex's 40-tau sampling | 35.353% |

Pointwise spatial convergence is spectral [C]. The maximum |dQ60| is 2e-7 between h = 0.25 and h = 0.125, and 7e-5
between h = 0.5 and h = 0.25.

**Sensitivity [C].** The trajectory is strongly sensitive. Scaling the initial amplitude by (1 + 1e-6) shifts the
window mean by -0.022 pp, and the pointwise difference saturates near 2.7e-4 after tau of about 4000. Differences
below about 0.03 pp (about 0.06% relative) therefore cannot count as convergence evidence. That includes Codex's
0.019% agreement between its two finest grids. It is consistent with convergence, but it does not certify accuracy
at that level.

**Diagnosis of bb211d1's 35.15% [C].** My solver with the old cubic damper (L = 360, start at 270, eta = 0.03
[Q-code]) gives 35.168% at dt = 0.015625. That reproduces N00625: its absorbed energy and charge, 23998 and 22023,
compare with 24002 and 22027 [Q]. The window slope also flattens, from -0.16 to -0.009 pp per 1000 tau.

The bb211d1 to ff49e5d shift is therefore attributable to soft-wave reflection from the old boundary, not to the
spatial Hamiltonian. Codex attributed it to both.

### 3.3 Beyond Codex's window

| Window | Mean, dt = 0.03125 | Mean, dt = 0.0625 (long sponge run) | Slope (pp per 1000 tau) |
|---|---|---|---|
| 9600 to 12000 | 35.346% | 35.142% | -0.162 |
| 19200 to 24000 | **34.729%** | 34.532% | -0.023 |
| 28800 to 36000 | | 34.354% (about 34.56% after the dt offset) [E] | -0.012 |
| 38400 to 48000 | | 34.243% (about 34.45% after the dt offset) [E] | -0.009 |

The dt offset is set during the violent phase (tau of about 1600 to 4000) and then stays constant: 0.204 pp at
12000 and 0.196 pp at 24000.

Q60 is still decreasing at tau = 48000. Power-law fits Q_inf + A tau^-beta are unstable: Q_inf falls from 34.0% to
33.55% and beta from 1.07 to 0.51 as the fit start moves from 10000 to 16000. **No asymptotic retention can be
certified.**

Codex's P025_remote value, 34.574% in 19200 to 24000 [Q], carries its h = 0.25 bias, which is about -0.17 pp. The
corrected value is 34.73% at dt = 0.03125, or about 34.74% extrapolated in dt.

### 3.4 Escape through r = 240 by tau = 12000 (h = 0.5, dt = 0.015625, box) [C]

| | This track | Codex I025 [Q] |
|---|---|---|
| Q_esc | 23526.9 | 23680.9 (+0.65%) |
| E_esc | 25486.1 | 25618.6 (+0.52%) |
| m Q_esc/E_esc | 0.92313 | 0.92436 |
| m Q_esc/E_charged | 0.93617 | |
| Neutral-s energy E_s | 354.7 (1.39% of E_esc) | |

Codex's I025 is the P025_time trajectory, whose retention is 35.19%, so its escaping totals inherit that
coarse-grid bias.

Escape is time-dependent and not asymptotic:
- m Q/E at r = 240 is 0.811 at tau = 2000, 0.865 at 4000 and 0.923 at 12000.
- Increments after tau = 4000 have dQ/dE of about 0.99.
- In the sponge runs, charge passing r = 240 (field plus receiver) reaches 24606 by tau = 24000
  (dt = 0.03125) and 25503 by tau = 48000 (dt = 0.0625), 65.4% of Q0.

**Free-field projection of the exterior field** (r > 240, smooth mask). The box keeps all escaped radiation, so
nothing is missing to absorption.

| Quantity | Value |
|---|---|
| N+ | 23750.6 |
| N- | 255.0 |
| E- | 432.3, a mean frequency of 1.70 per negative-frequency quantum |
| N_s, E_s | 110.2, 354.7 |
| N-/(N+ + N-) | **1.06%** |

This is stable to within 1% across masks at 240, 480 and 1000, and across h and dt. It is well below bb211d1's
positivity bound of 2.69%, and below the bounds from 0.92313 (3.84%) and from 0.93617 (3.19%).

The negative-frequency waves are emitted before tau of about 4000 and are fast (N- = 251 already at tau = 4000).
The late leakage is almost purely slow X: in the long run, N-/N+ is below 1e-4 and mean k is 0.03 to 0.09.

This is classical negative-frequency content [C]. Calling it an antiparticle number is [A]. Modes with p below
about 0.015 had not reached r = 240 by tau = 12000.

## 4. Verdicts

| Codex claim | Codex value | My value | Verdict |
|---|---|---|---|
| CFS reflection (ideal continuum) | 1.82e-8 / 6.5e-4 at p = 0.01 | identical, analytic [D] | CONFIRMED |
| CFS lattice reflection at 0.02 / 0.1 / 0.5 | as in 3.1 | own absorbers below 1e-9 at all three momenta; Codex's lattice solve not re-implemented | NOT CHECKED (lattice); comparable performance |
| Old damper reflection (9636950) | 31.7% at 0.02, 14.8% at 1 | 31.73%, 14.83% [C] | CONFIRMED |
| Retention, ff49e5d finest grid | 35.36% | 35.357 +- 0.02% [C/E] | CONFIRMED |
| Retention, bb211d1 | 35.15% | 35.36% (35.17% only with the old damper) | CORRECTED; cause: old-boundary reflection |
| Two finest grids agree to 0.019% | 0.019% | below the trajectory-sensitivity floor of about 0.06% relative | CONFIRMED as arithmetic; not a precision certificate |
| Moving the layer 360 to 720 changes retention by 4e-7 relative | 4e-7 | box vs sponge 1.8e-9 absolute [C] | CONFIRMED (absorber-insensitive) |
| Retention later, 19200 to 24000 (P025_remote) | 34.57% | 34.73% | CORRECTED |
| Still evolving, no asymptote | stated | slope -0.16 to -0.009 pp per 1000 tau through tau = 48000 | CONFIRMED and quantified |
| I025 escaping E, Q through 240 | 25618.6, 23680.9 | 25486.1, 23526.9 | CORRECTED (-0.5% / -0.65%); ratio 0.9231 vs 0.9244 |
| bb211d1 antiparticle fraction <= 2.69% (bound) | bound | 1.06% classical negative-frequency content | CONFIRMED (bound respected) |

## 5. Limits

- One spherical classical prototype; one prepared Gaussian and one parameter point.
- Q60 is a finite-time, window-dependent observable. Its intrinsic sensitivity is about +-0.02 pp.
- The dt = 0.015625 run exists only at h = 0.5. The h-convergence evidence is at dt = 0.0625 (h = 1, 0.5, 0.25;
  h = 0.125 to tau = 1956).
- The long-run levels beyond 24000 come from dt = 0.0625 plus a measured constant offset [E].
- The exterior projection is a free-field, equal-time decomposition; slow modes (p < 0.015) are not yet outside
  r = 240.
- No quantum particle spectrum, cosmological source, hot-core leakage or vacuum density is addressed.

Method credits: PML/CFS absorbers (Berenger; Kuzuoglu-Mittra; Roden-Gedney); symplectic splitting and Yoshida
composition; Q-ball formation and fragmentation (Kusenko-Shaposhnikov; Enqvist-McDonald; Kasuya-Kawasaki).

## 6. Reproduce

```bash
D=research/reality_framework/grut_formation_review_20261025/tracks/independent_solver/code
gcc -O3 -march=native -o <scratch>/radial_spec $D/radial_spec.c /lib/x86_64-linux-gnu/libfftw3.so.3 -lm
# (set BIN in drive_c.py to the binary path)
python3 $D/drive_c.py run --h 0.5 --L 6400 --dt 0.015625 --tmax 12000 --snaps 12000 --out <scratch>/A_h05_dt2
python3 $D/drive_c.py run --h 0.25 --L 6400 --dt 0.0625 --tmax 12000 --out <scratch>/A_h025
python3 $D/drive_c.py run --h 0.5 --L 4200 --dt 0.03125 --tmax 24000 --sponge 1200,3000,0.03,5 --out <scratch>/S_h05_long_dt
python3 $D/drive_c.py run --h 0.5 --L 360 --dt 0.015625 --tmax 12000 --sponge 270,90,0.03,3 --out <scratch>/D_olddamper_dt2
python3 $D/drive_c.py run --h 0.5 --L 6400 --dt 0.03125 --tmax 12000 --amp 0.18000018 --out <scratch>/A_h05_dt_pert
python3 $D/sponge_tests.py out.json ; python3 $D/sponge_tests.py out_soft.json soft
python3 $D/collect.py
```

The other runs in `results/formation_runs.json` vary --h, --dt and --tmax as listed there. The h = 0.5 runs take
about 0.5 to 2 hours on a loaded 4-core machine.

`spectral_radial.py` is the Python reference. It runs the same scheme directly:
`python3 spectral_radial.py --h 0.5 --L 6400 --dt 0.0625 --scheme yoshida4 --out X`.

Results are in `results/summary.json`, `results/formation_runs.json` and `results/absorber_tests.json`. Raw
histories and snapshots are in the scratchpad (`review/independent_solver/runs/`).

Caveat on the absorber file: the first two p = 0.02 packet rows in `absorber_tests.json`, from
`sponge_tests.json`, used a packet that was not odd-symmetrized; their error of about 1e-5 is very-low-k content.
They are superseded by the `soft` rows.
