# CP-09: spectral clock budget and why the autonomous apparatus does not yet prepare ideal relational time

**Date:** 2026-10-11. Independent bridge between the [frozen CP-09 law](https://github.com/ryangrvr/TestingGRUT/blob/1dadbd60ae5719508ac40ce7633768549052691a/research/constructive_physics/attempt09/LAW.md), [autonomous clock–pointer apparatus](https://github.com/ryangrvr/TestingGRUT/blob/99cd69493b3c11220c6f56bd19902d560b110d4a/research/reality_framework/cp09_autonomous_apparatus/REPORT.md), and [exact clock-conditioned response](https://github.com/ryangrvr/TestingGRUT/blob/3020ab9c2a28434a7e3c6584a17e04cd42b88fd8/research/reality_framework/cp09_clock_conditioned_response/REPORT.md). **No change to CP-09, no CP-10 or Stage-3 Card 3, no attempt to merge frozen packets.** The supplied Page–Wootters-type clock state is tested as a resource, not proposed as new physics.

## 1. Exact spectral-support theorem for the *noninteracting stationary constraint*

Let a finite system have Hamiltonian H_s and a q-level clock have Hamiltonian H_c. Suppose a pure joint state satisfies the exact stationary constraint

\[
(H_c\otimes I+I\otimes H_s)|\Psi\rangle=E_\ast|\Psi\rangle.
\]

Expanding in energy eigenspaces, a system energy E_n can have nonzero amplitude only if the clock has an eigenvalue \(\epsilon_a=E_\ast-E_n\). **Therefore at most q distinct system energy eigenvalues can appear** in any exact stationary state; degeneracies within one energy do not count separately. A clock with q<k distinct energies cannot reproduce for an interval of clock readings a specified conditional history whose unperturbed system state has k nonzero, distinct energy components. This assumes the clock and system are noninteracting in the constraint, and the conditional clock time states generate energy phases as in the published construction.

This is elementary finite-dimensional spectral matching. It is **not** a no-go for q-level clocks interacting with a system and pointer under a different Hamiltonian, for nonstationary autonomous evolution, or for approximate/time-limited conditional behavior.

For a given ground-state preparation |g> and a source B, the source-pulse history \(\exp(i\varepsilon B)|g\rangle\) has first-order tangent \(iB|g\rangle\). To reproduce the *entire first-order conditioned trajectory*, the stationary clock must be able to correlate all excited energies for which \(\langle n|B|g\rangle\ne0\), in addition to the ground energy. A single scalar response sampled at one time could accidentally match with a smaller clock; the theorem concerns the whole prepared conditional history.

## 2. Exact N=2 CP-09 spectral budget

The supplied CP-09 N=2 Hamiltonian has 12 eigenstates. For the frozen disjoint minimal-source B and maximal-detector A, an **independent** 12x12 construction finds six excited eigenenergies with nonzero source-pulse amplitudes, plus the ground level. The energy gaps and squared B matrix elements are:

| Excited eigenstate index (zero-based) | E_n-E_0 | \(|\langle n|B|g\rangle|^2\) |
|---:|---:|---:|
| 2 | 0.3881916254 | 0.1500014541 |
| 3 | 0.7216272304 | 0.0131772514 |
| 6 | 1.4708319081 | 0.1783721202 |
| 8 | 2.9961143834 | 0.0319341483 |
| 9 | 3.4183737513 | 0.0662245298 |
| 11 | 3.8415741156 | 0.0009087024 |

The other five excited eigenstates have vanishing source overlaps within 1e-11 in the numerical control. The six nonzero weights sum to the source variance

\[
{\rm Var}_g(B)=0.4406182061422456.
\]

If a **two-energy stationary clock must include the ground energy** so it can also reproduce the initial, unperturbed ground state, it can include at most one of the six excited energies. Choosing the largest overlap still omits first-order source tangent norm

\[
\boxed{\sqrt{{\rm Var}_g(B)-\max_{n>0}|\langle n|B|g\rangle|^2}=0.5120996835836915.}
\]

This is 77.15% of the excited tangent norm; it is a **state-level** lower bound for any such two-energy exact spectral restriction, **not** a lower bound on one scalar readout's absolute response error. A specific 2-level truncation with ground plus highest-weight excitation gives response at t=.4 of 0.129836435962, versus ideal 0.052588513539. That example is diagnostic, not an optimization of the scalar response.

Because B is a block-diagonal Pauli source that is zero on the empty geometry and squares to the nonempty-geometry projector, B²|g> lies in precisely the same seven-energy sector. Therefore the **entire finite source-pulsed state** exp(±i ε B)|g> stays in this seven-dimensional spectral support, not only its tangent. The independent test checks this at ε=1e-5 and also verifies B²|g>.

## 3. Exact seven-level clock reproduces the finite response

For the seven occupied eigenstates \(S=\{0,2,3,6,8,9,11\}\), supply a seven-level clock with energies \(\epsilon_n=-E_n\) and a normalized stationary history \(\sum_{n\in S} a_n |n\rangle_c\otimes|n\rangle_s\). Clock test kets \(|t\rangle=7^{-1/2}\sum_{n\in S}e^{-i\epsilon_n t}|n\rangle\) yield exactly \(e^{-iH_st}|\psi(0)\rangle\) after conditioning, with rank-one test probability 1/7. The independent implementation reproduces the established CP-09 disjoint-probe impulse response:

| Clock reading t | q=7 conditioned impulse response | Ideal Kubo response |
|---:|---:|---:|
| 0.1 | 0.004121694443 | 0.004121694452 |
| 0.4 | 0.052588513531 | 0.052588513539 |
| 1.0 | 0.029091102669 | 0.029091102668 |

State-vector residuals are below 1e-11 and the stationary constraint is satisfied within numerical tolerance. This improves the **resource count of the existing conditional representation**, from 12 clock energies to 7 for this one prepared experiment. It is *not* an experimentally accessible physical seven-state clock, and it still requires complete knowledge of the selected system spectrum and correlated initial state. It also does not turn separate pre-prepared \(\pm\varepsilon\) histories into an intervention performed dynamically at a clock reading.

## 4. Why the actual autonomous qubit apparatus is different

The [autonomous apparatus](https://github.com/ryangrvr/TestingGRUT/blob/99cd69493b3c11220c6f56bd19902d560b110d4a/research/reality_framework/cp09_autonomous_apparatus/REPORT.md) has

\[
H_{\rm tot}=H_s+H_c+V,\quad
H_c=(I-X_c)/2,\quad V=P_{1,c}\otimes F_s\otimes(I-X_r).
\]

It uses **a two-level clock interacting with both geometry and a record qubit**, and its initial state is not the stationary spectrally matched Page–Wootters state. The above noninteracting spectral theorem is NOT a theorem that this 48-dimensional autonomous apparatus cannot generate partial, imperfect, or even eventually accurate conditional records under a separately specified protocol. Its clock is dynamically driven: \([H_c,V]\neq0\), and the relative clock/system/pointer energies actually change in the published experiment. The global F_s remains a supplied geometry-controlled coupling.

For a truly *stationary interacting* constraint \((H_c+H_s+V)|\Psi\rangle=E_\ast|\Psi\rangle\), with covariant clock bras \(\langle t|=\langle0|e^{iH_ct}\) and unnormalized relative states \(\phi(t)=\langle t|\Psi\rangle\), direct differentiation gives

\[
\boxed{ i\partial_t\phi(t)=(H_s-E_\ast)\phi(t)+\langle t|V|\Psi\rangle.}
\]

That extra term is a **physical interaction source**, not automatically zero or a scalar phase. Its time-nonlocal interpretation is established in the interacting relational-clock literature, e.g. Smith & Ahmadi, *Quantizing time: Interacting clocks and systems*, Quantum 3, 160 (2019), https://doi.org/10.22331/q-2019-07-08-160 . The equation does not by itself solve for a physical clock POVM, a preparation mechanism or the conditional multi-time detector protocol.

The new result therefore separates an exact **representation of a preselected history** from an apparatus that *prepares, triggers and reads* histories by its own energy-conserving interactions.

## 5. What must be constructed next

A serious physical bridge would specify a fixed finite-resource clock/reference/pointer Hamiltonian with independently calibrated initial preparation and **no matched lookup of all CP-09 eigenenergies**. It must define clock readings as a genuinely normalized instrument over an operational window, implement the source trigger at a clock reading within one joint history, and evaluate the response conditioned on the actual pointer and clock records. It should report the total closed-system energy ledger, disturbance, clock accuracy, usable duration, and error versus a predeclared ideal response for N=2,3,... under a controlled scaling.

Crucially, merely adding more clock levels or a spectrally tuned history state does not repair CP-09's all-pairs escape scale, emergent localization, manifold selection, massless spin-2 dynamics or the physical meaning of its labels. No new CP candidate is introduced here.

## Reproduction

Run [check.py](check.py) with NumPy and SciPy. It independently constructs the frozen CP-09 N=2 Hamiltonian, computes the ground state and relational source/detector, verifies the six excited overlaps and exact seven-level support, constructs q=7 stationary conditional histories, and checks three conditioned impulse responses against ideal Kubo values and the q=2 example. Executed independently before publication. Expected terminal: \`PASS: seven matched clock energies recover exact response; two cannot reproduce the full first-order clock history\`.

Other contextual prior art: Page–Wootters relational quantum dynamics reviewed in Höhn, Smith & Lock, *Trinity of relational quantum dynamics*, Phys. Rev. D 104 066001 (2021), https://doi.org/10.1103/PhysRevD.104.066001 . These are previously known mathematical mechanisms; the present spectrum count and exact control refer only to this supplied CP-09 toy.
