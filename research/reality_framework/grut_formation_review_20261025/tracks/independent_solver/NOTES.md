# Track `independent_solver`: an independent radial solver for the prepared formation run

Claude review of Codex packets ff49e5d (causal receiver) and bb211d1 (core conversion), 25 October 2026.
This track wrote its own solver from the equations in Codex's MATHEMATICS.md files (c4b99f3 section 4 and
bb211d1 section 2). It did not run Codex's code. Codex's numbers below are quoted [Q] from the packet JSON
files, which were read through `git archive` into the scratchpad. Labels: [D] derived, [C] computed here,
[Q] quoted, [E] estimate, [A] argued.

STATUS: DRAFT, numbers being filled in as runs finish.

## 1. What was solved

The equations, from MATHEMATICS.md [Q]: tau = g_X v t, rho = g_X v r, X = v psi/sqrt2, phi = v S, with

    psi_tt = lap psi - (S^2 + 0.1|psi|^2) psi,      S_tt = lap S - [3(S^2 - 1) + |psi|^2] S,
    Qbar = -int Im(psi* psi_t) d^3rho,
    Ebar = int [ |psi_t|^2/2 + S_t^2/2 + |grad psi|^2/2 + |grad S|^2/2 + (3/4)(S^2-1)^2 + S^2|psi|^2/2 + (0.1/4)|psi|^4 ].

Prepared data: psi = 0.18 exp(-rho^2/7200), psi_t = -i psi, S_t = 0, with S taken from the static neutral
equation, S'(0) = 0 and S(inf) = 1. My initial data give Ebar0 = 38951.47049696 and Qbar0 = 38969.38665303.
Codex records 38951.470496956 and 38969.386653027 [Q]. The two agree to 11 significant digits [C].

## 2. Method, chosen to differ from Codex

| | Codex ff49e5d | This track |
|---|---|---|
| Variables | U = r psi, Z = r(S-1) | Same (standard) |
| Space | 4th-order staggered finite differences, h = 0.25 to 0.0625 | **Fourier-sine (DST-I) pseudo-spectral**, h = 1, 0.5, 0.25 (0.125 partial) |
| Time | RK4, dt = 0.015625 (refined) or 0.03125 | **Symplectic splitting**: exact free Klein-Gordon rotation of each sine mode plus a nonlinear position kick; Yoshida 4th-order composition; dt = 0.0625, 0.03125, 0.015625 |
| Boundary | CFS-PML with memory, from r = 360 to 450 | **(a) A box too big for any reflection to return** (hard wall at L = 6400; a signal leaving r < 600 at tau = 0 cannot return to r = 60 before tau = 12140). **(b) Smooth velocity sponge** gamma = 0.03 ((r-R_s)/3000)^5 on [R_s, R_s+3000], R_s = 1200, used only for the long runs |
| Static S | solve_bvp on the continuum | Spectral fixed point on the actual grid (agrees with an independent solve_bvp of the Z form to 7e-14) |
| Core charge Q60 | nodal sum | Trapezoid plus Euler-Maclaurin end correction (h^2/12)q'(60), using spectral q' |

Charge is conserved to round-off: the kick multiplies by a real factor, and the free rotation is unitary.
Measured max |dQ/Q0| < 2e-12 [C]. Energy error is bounded because the scheme is symplectic.

Code (all in `code/`):
- `spectral_radial.py` is the reference implementation (Python).
- `radial_spec.c` is a fast C/FFTW port of the same algorithm. It agrees with the Python version to 1e-14
  relative over tau = 400, with and without the sponge [C].
- `drive_c.py` prepares the data, runs the C solver and post-processes the output.
- `sponge_tests.py` holds the absorber tests; `collect.py` writes `results/`.

## 3. Results

(filled in below)
