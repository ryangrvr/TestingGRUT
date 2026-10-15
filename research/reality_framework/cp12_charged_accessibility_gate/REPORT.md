# CP-12 charged-accessibility gate

Parent: [a2f447d](https://github.com/ryangrvr/TestingGRUT/commit/a2f447dc37ea0a71f62c6e990a4630d727ebd85e).

**The charged route passes a strict finite-system accessibility test. Appreciable macroscopic weight, persistence, dimension selection, and local propagation remain open.**

The initial state is the original coherent neutral ray. The frozen Hamiltonian and all its parameters are unchanged. There is no source kick, initial loop filter, spatial host, or phase mask.

## 1. A precisely defined target

For N >= 5 and Q=2N, let G_N contain every physical configuration satisfying:

- The graph is connected and four-regular, with no preferred labelled host.
- Every active edge belongs to the union of directed electric circuits.
- n_v >= ceil(N/4).

The graph has N+1 independent cycles. This is a label-invariant, extensively cyclic charged target; it is not a definition of three-dimensional space. All upper bounds below also hold for the larger event that drops the graph-degree and connectivity requirements.

## 2. Exact access from the actual neutral ray

At the smallest allowed size N=5,Q=10, every graph is K5. Complete shortest-path amplitude summation gives

\[
\boxed{
p_{G_5}(T)=
\frac{1326415312}{848779513441154588027109375}\,T^{44}
+O(T^{46}).
}
\]

The coefficient is approximately 1.5627324776282548 x 10^-18. It is strictly positive.

All target amplitudes vanish before 22 hops. At order 22 there are exactly 960 configurations: 24 regular tournament orientations and 40 admissible pairs of edge holes per orientation. They form eight signed-permutation orbits of size 120. The two absolute amplitude values are 49,583,226,880 and 40,677,670,912, each occurring in four orbits. The exact norm square is

\[
\|P_{G_5}H^{22}\psi_0\|^2
=1974321263356092214149120.
\]

Independent forward and reverse BigInt calculations agree for every orbit. All 960 labelled forward amplitudes were also calculated directly, confirming the orbit sum. The same individual target was also checked in Python and JavaScript.

This is a sum of **signed quantum amplitudes**, not a count of allowed paths. It proves that the original ray's cyclic subspace has nonzero target overlap. A deliberately conservative explicit bound is

\[
\|P_{G_5}\Pi_{\mathrm{cyc}(\psi_0)}\|^2
\ge 3.2723\times10^{-42}.
\]

That number bounds an operator norm; it is not a transition probability or a lower bound on a dominant phase. The short-time coefficient establishes positivity for sufficiently small nonzero T. It is not extrapolated to T=1 or long times.

K5 is already a fixed complete graph because Q=M. This exact test cannot demonstrate graph selection, emergent spatial dimension, or the action of rewiring on a choice of graphs.

## 3. Collective probability, rather than hopping order alone

Write e_C for the number of edges belonging to directed electric circuits. Every configuration with e_C >= L and n_v >= k has electric L1 norm at least

\[
d=2L+k.
\]

Rewiring and the diagonal Hamiltonian preserve that norm; a hopping event changes it by one. The full, unbounded-electric-field Hamiltonian obeys

\[
b=|t|\sqrt{2Q\min(N,Q)},\qquad \|H_{\rm hop}\|\le b.
\]

For the specified neutral ray the first hopping norm is exactly sqrt(2Q)|t|. A unitary Duhamel remainder then gives, for d >= 1,

\[
\boxed{
P_T(e_C\ge L,n_v\ge k)
\le
\min\!\left\{1,\,
\frac{1}{\min(N,Q)}
\left[\frac{(b|T|)^d}{d!}\right]^2
\right\}.
}
\]

The bound includes interference and every possible transition path, and has no electric cutoff. It remains valid for arbitrary T, although it becomes vacuous outside its useful early-time regime.

For G_N, d=4N+ceil(N/4) and b=2N|t|. At t=1 there is exponential suppression in N for each fixed

\[
|T|<\frac{4.25}{2e}=0.781743812489315\ldots .
\]

Examples at T=0.5:

| N | Analytic probability upper bound |
|---:|---:|
| 16 | 5.852 x 10^-31 |
| 64 | 2.363 x 10^-111 |
| 1024 | log10 upper bound = -1696.845262 |

Values below floating-point range are stored logarithmically, never as a zero probability.

An extensive minimum hop count therefore **does not establish a formation time proportional to N**. The collective hopping norm also grows with N. This is a bound on distance through configuration space, not a causal cone in physical space. At T=1 the general bound is vacuous.

## 4. A spectral bound using the initial ray

A separate N=4,Q=4 rewiring benchmark uses all 577 retained signed S4 modes, with coefficients taken from the actual neutral ray and no mode pruning. Its target is full circulation with n_v >= 1, a charged spanning ring. This is a separate benchmark, not a second size of the degree-four family.

For physical trial modes v_a, real frequencies lambda_a, and their actual initial coefficients c_a, define

\[
A=\sum_{a,b}|c_a^*c_b\langle v_a|P_G|v_b\rangle|,
\quad
\epsilon(T)=\left\|\psi_0-\sum_a c_av_a\right\|
+|T|\sum_a|c_a|\|(H-\lambda_a)v_a\|.
\]

The residual uses the full unbounded H, including every outgoing transition. Then

\[
p_G(T)\le \min\{1,(\sqrt A+\epsilon(T))^2\}.
\]

The executed numerical calculation obtained

- A = 0.14752894104361644.
- Initial reconstruction norm = 4.100982338017993 x 10^-14.
- Weighted full-H residual = 1.4352906680026263 x 10^-13.

Thus the evaluated full-law upper bound is below 0.147528942 throughout 0 <= T <= 10. At T=5 the target probability is approximately 0.055022107666, with a truncation/residual interval [0.05502210766561084, 0.05502210766632265].

This exploits the initial ray's occupied trial directions and target Gram matrix. It does not replace them with an unrestricted energy-window norm. The finite trial signal has an all-time envelope; the full-H estimate retains its explicit time-dependent error.

## 5. Capacity-valid five-vertex upper bound

The N=5,Q=10 retained-space calculation has electric-square budget 32, dimension 39,969, and all 960 minimum-order charged targets within it. Their electric-square requirement is 26. Higher charged targets outside that budget are included in the full-law error, not declared impossible.

| T | Retained target probability | Full-state error upper bound | Full target probability upper bound |
|---:|---:|---:|---:|
| 1 | 2.126 x 10^-19 | 8.540 x 10^-7 | 7.300 x 10^-13 |
| 2 | 2.259 x 10^-9 | 0.003164 | 1.032 x 10^-5 |
| 3 | 3.745 x 10^-7 | 0.04977 | 0.002539 |

Every numerical lower interval here includes zero. The exact onset proof supplies nonzero accessibility; these finite-time numerical estimates supply upper bounds and do not independently resolve a positive probability.

The error bounds in Sections 4-5 are analytically derived and evaluated in ordinary floating point. They are not formal outward-rounded interval-arithmetic certificates.

## Verification and next gate

The exact core, all eight forward/reverse orbit calculations, 20 analytic-bound evaluations, and 680 independent finite-unitary sanity cases were executed in the orchestration JavaScript runtime. The numerical session commands completed before the container disconnected. The packaged Python reproduction wrapper has not been rerun; this packet does not claim a full reproduction-suite PASS or external mathematical certification. Details are in VERIFICATION.json.

The next useful calculation is **charged target weight and persistence at Q=2N with N>5**, where graph choice exists, using the same unfiltered state and full-law error control. A growing charged network would then need a joint test of intrinsic volume growth and spatial response. The present result establishes neither a preferred dimension nor extended causal propagation.

The original coherent neutral preparation remains a specified input; a physical mechanism that autonomously prepares it has not been derived.
