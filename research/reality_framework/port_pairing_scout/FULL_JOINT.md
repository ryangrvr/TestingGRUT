# SCOUT-RG02 full quantum matter–geometry extension: no Born–Oppenheimer shortcut

**2026-10-11**, separate continuation of port-pairing SCOUT-RG02, not a modification of frozen CP-11 or any Stage-3 Card. The companion [full_joint.py](full_joint.py) solves the complete finite, exactly number-conserving Hamiltonian using the fermion canonical anticommutation relations and all connection matchings.

## 1. Fixed fermionic Hilbert space while geometry changes

The earlier matching calculation treated matter through its exactly derived ground energy, then used a quantum graph-flip Hamiltonian. That was an effective Born–Oppenheimer approximation. To remove it, place **two persistent fermion modes on the two ports of each vertex**. With 2N physical modes and N spinless fermions, the fermion Hilbert dimension is \(\binom{2N}{N}\), the same for every graph. Each perfect port pairing P still labels one geometry state, but graph reconnection can now act as identity on the *same* port-fermion space without inventing an identification of edge-fermion modes after a graph changes.

Let \(c_{i,a}\) be canonical fermionic annihilation operators on port a=0,1 of site i. Supply the time-independent joint Hamiltonian
\[
\boxed{
H_{\rm joint}=
\sum_P |P\rangle\langle P|\otimes H_f(P)
+\frac{\Gamma}{N}\sum_{\{P,P'\}\in\mathcal E}
(|P\rangle-|P'\rangle)(\langle P|-\langle P'|)\otimes I_F
+\frac{4tN}{\pi}I,
}
\]
with \(\Gamma=.1,t=1\), total fermion number N, and
\[
\boxed{H_f(P)=t\sum_{i=1}^N(c_{i,0}^\dagger c_{i,1}+{\rm h.c.})
+t\,s\sum_{(p,q)\in P}(c_p^\dagger c_q+{\rm h.c.}),\quad s=\pm1.}
\]
Here the reconnection links are *pairs of matching configurations*, not physical spatial bonds. The relative sign s is a **supplied microscopic parameter with a physically nontrivial flux effect** on odd-length graph cycles. We test both signs adversarially; no after-the-fact choice is hidden.

For the s=-1 sector, each connected l-cycle of vertices becomes a 2l-site alternating port chain with one-body bands
\[
\omega_\pm(q)=\pm t|1-e^{iq}|=\pm2t|\sin(q/2)|.
\]
Thus the half-filled N-fermion ground energy on each matching is **exactly** the massless incidence model's
\[
E_{\rm fill}(P)=-2t\sum_{l\in {\rm cycles}(P)}\cot(\pi/(2l)).
\]
In particular the positive Hamiltonian lower bound follows WITHOUT adiabatic elimination:
\[
H_{\rm joint}\ge0
\]
at s=-1 after the displayed common +4tN/pi shift, because every conditional fermion block has lowest eigenvalue \(-4tN/\pi+\sum_l\delta_l\ge -4tN/\pi\) and the graph reconnection Laplacian is nonnegative. This establishes finite-N ground-energy stability with no negative O(N²) collapse.

For s=+1 the band has \(\pm t|1+e^{iq}|\). Self-loop vertices permit intra-port and inter-port hoppings to add constructively, sharply favoring disconnected l=1 components rather than long loops. The two signs differ in gauge-invariant product around an odd-l cycle, so one cannot gauge them into one another throughout the full matching ensemble. This is a **physical flux-sector selection assumption, not its dynamical explanation**.

## 2. Full finite quantum calculation

All matching configurations and all \(\binom{2N}{N}\) N-fermion Fock states are retained. The hopping matrix is built with exact Jordan–Wigner/wedge anticommutation signs. There is no external source schedule, premanufactured history state, selected ring, Born–Oppenheimer kinetic projection, or manual matching graph postselection.

| N | Graph matchings | Fermion Fock states | Full joint dimension | P(connected), s=-1 | P(connected), s=+1 |
|---:|---:|---:|---:|---:|---:|
| 3 | 15 | 20 | 300 | **0.995020829** | 0.000004128 |
| 4 | 105 | 70 | 7,350 | **0.985452621** | 0.000000009997 |
| 5 | 945 | 252 | 238,140 | **0.968940897** | not executed |

The sign s=-1 is declared before comparing to these outputs, using the same \(\Gamma=.1,t=1\) for every N. The s=+1 counterexample is not concealed. Negative sign results are close but NOT identical to the effective-geometry quantum probabilities 0.994388878 (N=3), 0.982037359 (N=4), 0.959993501 (N=5), which demonstrates that the effective model did omit coherent matter dynamics.

Lowest full joint energy after the common shift:
- N=3: **0.456749763** (s=-1), **-1.983674562** (s=+1).
- N=4: **0.419356126** (s=-1), **-2.610844672** (s=+1).
- N=5: **0.417896147** (s=-1).

The s=+1 energies are allowed negative because the +4tN/pi shift is the lower-bound normalization proved only for s=-1.

## 3. Genuine fermion–geometry quantum entanglement

For the normalized eigenstate \(|g\rangle\in\mathcal H_{\rm geometry}\otimes\mathcal H_F\), calculate the Schmidt coefficients and reduced geometry purity:
\[
\rho_G=\mathrm{Tr}_F|g\rangle\langle g|,\qquad
\mathcal P_G=\mathrm{Tr}(\rho_G^2),\qquad
S_G=-\mathrm{Tr}(\rho_G\log\rho_G).
\]

| N | Geometry purity \(\mathcal P_G\) | Entanglement entropy S_G (nats) | Schmidt rank |
|---:|---:|---:|---:|
| 3 | 0.772249 | 0.529023 | 7 |
| 4 | 0.561881 | 1.173879 | 27 |
| 5 | 0.404441 | 1.869061 | 111 |

The state is **not** a product of a chosen geometry with a passive matter state. Its connected probabilities reflect the actual entangled joint ground state.

We tested signed vertex transposition covariance at N=3,4. The permutation of fermion port modes acts with the correct wedge sign on each occupied Fock string, and the exact **full-Hamiltonian commutator** vanishes within numerical precision (Frobenius norm 0). Eigenstate residual norms are below 6e-15 for the above tests.

## 4. Complete quantum finite-temperature and relational detector controls

The full-Hilbert-space Gibbs trace at **N=3** can be computed exactly from all 300 joint eigenstates, with no BO elimination. [full_joint_thermal_n3.py](full_joint_thermal_n3.py) produces:

| Inverse temperature beta | P(connected), s=-1 | P(connected), s=+1 |
|---:|---:|---:|
| 0 | 0.533333 | 0.533333 |
| 1 | 0.692101 | 0.369106 |
| 5 | 0.997022 | 0.000380 |
| 20 | 0.997740 | 0.000004128 |

These are actual canonical quantum thermal probabilities for **one small N**. There is no controlled large-N joint-Gibbs extrapolation and no proof that finite positive temperature sustains global connectivity.

A more operational test is available without projecting the graph onto a chosen classical ring. Anchor a source to an actual physical fermion port 0: \(B=n_0\). In each matching P, define the detector \(A_1(P)=n_{\mathrm{mate}_P(0)}\). For the two-link detector, follow the first port bond, cross its vertex to the other port, and follow the second port bond; \(A_2(P)=n_{\mathrm{mate}_P(\mathrm{mate}_P(0)\oplus1)}\). These probes are diagonal in graph and fermion occupation, commute at equal time with B, and depend on **local bond traversals**, not a postselected classical graph. Their realization by a physical moving/readout apparatus is not derived—following an actual edge still requires a prepared local reference.

The full joint ground-state Kubo response \(R_{AB}(t)=i\langle[A(t),B]\rangle\) was calculated with the entire matrix exponential and checked by independent source kicks \(\exp(\pm i10^{-5}B)\):

| N | t | one-link response | two-link response | fixed other-vertex port response |
|---:|---:|---:|---:|---:|
| 3 | 0.4 | -0.038040 | -0.000829 | -0.017842 |
| 4 | 0.1 | -0.083989 | -0.0000592 | -0.014136 |
| 4 | 0.2 | -0.160504 | -0.0004745 | -0.027532 |
| 4 | 0.4 | **-0.266133** | **-0.00381282** | -0.049455 |

At N=4, an independently evaluated **exact first-time derivative** gives one-link \(R_1'(0)=-0.85262727\), while the two-link derivative vanishes numerically below \(10^{-9}\); the two-link response is consistent with leading **cubic** time dependence, \(R_2(t)/t^3\simeq-0.0591\) at t=.02,.04,.08. This is a finite-ground-state selection result, **not** a universal graph light cone: short loops/shortcuts, reconnection, changing preparations and longer times need not preserve that hierarchy.

The calculation demonstrates an operationally meaningful distinction between one- and two-step **relational** detector observables while every geometry fluctuates and is coupled to fermions. It improves on CP-09's label-only and globally determined role probes. It does NOT prove that a laboratory observer can implement these graph-controlled readouts without a supplied access mechanism, and does not fix any physical distance/time units.

Executable: [relational_response.py](relational_response.py); source-sign finite impulse comparisons pass within \(10^{-8}\).

## 5. What this construction genuinely resolves—and what it still does not

**New mathematical progress relative to the earlier SCOUT:** the mechanism now supports a fixed physical matter register and a fully coherent quantum geometry changing in time, with the same port pairing controlling fermionic hopping and the energy that favors connected cycles. It has an extensive lower bound and robust small finite-N connectedness, and explicitly exhibits matter–geometry entanglement. This addresses the Born–Oppenheimer and edge-register identification concerns left by the earlier conditional and quantum-effective benchmarks.

**Not resolved:**
- The preferred \(\mathbb Z_2\) sign/flux sector **is an input**, not dynamically selected. The opposite sign selects fragmented loops almost perfectly.
- The microscopic two-port rule builds in degree two, so this is at most an emergent **one-dimensional** ring, not derivation of three space dimensions.
- Reconnection acts on arbitrary pairs of bonds irrespective of their separation in a candidate emergent geometry. Operational interaction locality has not been earned.
- A coherent state over relabelled rings gives no separation profile to identically defined fixed-label detectors without a physical reference. Actual relational source and detector interactions, clocks and energy accounting remain unconstructed.
- The full quantum **finite-temperature** fermion-plus-geometry Gibbs trace and its large-N limit are NOT yet calculated. The earlier conditional Gibbs theorem and the zero-mode ground-manifold control warn of fragmentation but cannot simply be advertised as rigorous theorems for this complete Hamiltonian.
- There is no macroscopic 3+1-dimensional Lorentzian state, graviton, universal metric coupling, Einstein equation, or new empirical prediction. Finite N=3..5 does not establish a surviving thermodynamic phase.

Do not infer that one sign is selected by reality merely because its small systems look better. A gauge/flux dynamical principle must earn that sign, select dimensionality, survive finite temperature, and restore operational locality **within the same law**. The comparison is a reproducible falsifiable hypothesis, not a completed Theory of Everything.

### Reproduction and verification

Run \`OPENBLAS_NUM_THREADS=1 python research/reality_framework/port_pairing_scout/full_joint.py\`. It constructs the 300-, 7,350- and 238,140-dimensional joint sparse Hamiltonians, diagonalizes the lowest modes, checks covariance with signed fermion permutations, and outputs connectedness and Schmidt entanglement. Equivalent independent implementation and controls were executed locally. The exact committed script is still available for a verbatim partner rerun.

Prior-art context: ordinary fermionic tight-binding models, gauge-invariant cycle flux, CAR Fock-space signs and quantum graph geometry are established ingredients. No independent priority claim is made for this assembly.
