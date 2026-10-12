# SCOUT-RG02: a ring is energetically selected without supplying its ordering—and why the thermodynamic universe still fails

2026-10-11. Independent sustained theory-construction result, under the user's Einstein challenge. New branch from prior SCOUT-RG01 benchmark (8813abc); frozen CP-11, GRUT-RAI, and Stage-3 Card 3 untouched. **This is not a finished theory of reality, CP-12, or a claim of scientific priority.** Full reproducible [check.py](check.py) accompanies this note.

## 1. The physical hypothesis removes an assumed circle

Instead of assigning each particle a position *on an already-existing S1*, take N vertices with exactly TWO connection ports each. Allow every possible perfect matching of the 2N ports. These matchings make degree-two multigraphs: any partition into disconnected cycles of sizes l_1+...+l_k=N. Even self-loops (l=1) and double edges (l=2) are included.

There are \((2N-1)!!\) matching basis states; a connected N-cycle is **not preselected**. Port number two, the permutation-invariant matching measure and elementary quantum degrees of freedom remain inputs. This does not derive physical dimensionality in general: the possible graphs are essentially unions of one-dimensional circles by kinematics. The full model is permutation covariant in the labelled vertices and ports, but whether relabellings should be gauged requires a supplied reference interpretation.

For each configuration G, let its vertex–edge incidence fermions determine energy **and** conditional propagation by massless Dirac hopping D_G=t[[0,B_G],[B_G^dagger,0]], with N fermions filling negative levels and suitable zero modes. At t=1, with common edge cost κN,
\[
U(G)=\kappa N-\mathrm{Tr}\sqrt{L_G}.
\]
The **massless** choice m=0 is a NEW physical assumption (frozen CP-11 had m=1). The zero modes imply degeneracies that would matter in a full finite-temperature fermion partition function. We only use the *filled-band minimum* to construct the conditional geometry energy.

## 2. A genuine exact zero-temperature geometry-selection result

For an l-cycle (including allowed self-/double-edge cases), the graph Laplacian has
\[
\lambda_k=4\sin^2(\pi k/l).
\]
Therefore its fermionic filled-band energy is **exactly**
\[
\boxed{U_l=\kappa l-2t\cot(\pi/(2l)).}
\]
Define the *excess above the common bulk energy*
\[
\boxed{\delta_l=\frac{4t}{\pi}l-2t\cot(\pi/(2l)).}
\]
For continuous l>=1, let x=π/(2l). Differentiation gives
\[
\delta'_l=\frac{4t}{\pi}\left[1-\left(\frac{x}{\sin x}\right)^2\right]<0,
\qquad \delta_l>0,\quad \delta_l\to0.
\]
Any partition l_1+...+l_k=N with k>=2 obeys
\[
\sum_{j=1}^k\delta_{l_j}>\delta_N.
\]
This is strict because every l_j<N and δ_l is positive decreasing. Hence, **at every finite N>=2, the conditional fermionic ground-energy minimizers among ALL 2-port pairings are precisely the connected N-cycles**, all relabellings included. We did NOT insert a global connectivity term, a specific cyclic order, an embedding, or a lattice-energy reward.

The long-wavelength one-body bands on the *dynamically preferred fixed-graph class* are
\[
\omega_\pm(q)=\pm2t|\sin(q/2)|=\pm t|q|+O(q^3).
\]
Thus *conditional* matter excitations on the resulting ring see the same neighbor links that its binding energy favored, with maximum group velocity t in unit edge spacing. This is a conventional 1D Dirac/phonon-like band, **not** a physical spacetime light cone or spin-2 gravity. Full joint graph-matter causal propagation has not been derived.

### The hidden gap failure

The selection becomes arbitrarily soft for a split into two macroscopic cycles:
\[
\delta_l=\frac{\pi t}{3l}+O(l^{-3}),\qquad
2\delta_{N/2}-\delta_N=\frac{\pi t}{N}+O(N^{-3})
\]
for even N. So the energy difference between one ring and two equally macroscopic rings **vanishes like 1/N**. This is much weaker than a size-extensive energetic preference for a connected universe. It warns against declaring thermodynamic spontaneous spatial selection from finite-N minimizers alone.

## 3. Exact finite-temperature failure for the conditional geometry Gibbs ensemble

The combinatorial count for a matching with c_l cycles of each length l is
\[
\#(c_1,c_2,\ldots)=
N!\prod_{l\ge1}\frac{2^{(l-1)c_l}}{l^{c_l}c_l!},
\qquad \sum_l l c_l=N.
\]
Consequently, after removing a common factor 2^N N! and the bulk fermionic energy, the conditional-geometry partition function is
\[
\boxed{Z_N(\beta)=[z^N]\exp\left(
  \frac12\sum_{l\ge1}\frac{e^{-\beta\delta_l}}{l}z^l
  \right).}
\]
The exact connected probability is
\[
\boxed{P_\beta(\text{one N-cycle})=
\frac{e^{-\beta\delta_N}}{2N\,Z_N(\beta)}.}
\]
The constant-beta recursive evaluator \(N Z_N=\frac12\sum_{l=1}^N e^{-\beta\delta_l}Z_{N-l}\) is independently implemented and tested at beta=0,1,5,20.

At beta=0 these are the standard random-pairing / Ewens(theta=1/2) cycle counts. This is known mathematical physics, not claimed novelty. At any **fixed finite beta**, δ_l=πt/(3l)+O(l^-3), so
\[
F_\beta(1)=\exp\Big[
 \frac12\sum_{l=1}^\infty \frac{e^{-\beta\delta_l}-1}{l}
 \Big]
\]
is finite and strictly positive. Standard coefficient asymptotics give
\[
Z_N(\beta)\sim F_\beta(1)/\sqrt{\pi N},\quad
\boxed{P_\beta(\text{one connected ring})
\sim\frac{\sqrt{\pi}}{2F_\beta(1)}N^{-1/2}\longrightarrow0.}
\]
The thermal suppression is **polynomial**, not the exponentially strong fixed-bond graph-entropy no-go, but is still decisive for a single connected universe.

| N | beta=0 | beta=1 | beta=5 | beta=20 |
|---:|---:|---:|---:|---:|
| 16 | 0.223294 | 0.423471 | 0.792106 | 0.992893 |
| 64 | 0.110995 | 0.213897 | 0.448140 | 0.787488 |
| 256 | 0.055416 | 0.106785 | 0.227865 | 0.445593 |
| 1,024 | 0.027698 | 0.053314 | 0.113845 | 0.226588 |

This explicitly demonstrates a potentially misleading finite-size success: even at inverse temperature beta=20, where N=16 appears >99% connected, the probability eventually falls with N. The statement concerns the **Born–Oppenheimer conditional filled-band geometry Gibbs ensemble**, not the full fermionic Gibbs ensemble. Extra fermion zero-mode multiplicities and any dynamically quantized connection field must be included in a complete finite-T test.

Typical configurations retain cycles whose sizes are O(N): a giant component can exist without the complete geometry being one connected component. Do not equate these claims.

### An additional adversarial factor: exact massless-fermion zero-mode degeneracy

The conditional energy calculation fills all negative one-body levels. For a graph with c disjoint cycles, the incidence matrix rank is N-c. Its 2N-dimensional massless Dirac operator therefore has N-c negative eigenvalues and **2c exact zero eigenvalues**. At fixed N fermions, filling the negative levels requires filling c of these 2c zero-energy states. Consequently the minimum-energy fermionic sector has
\[
\boxed{d_{\mathrm{zero}}(c)=\binom{2c}{c}}
\]
distinct Fock-state fillings, even before dynamical graph reconnection. This degeneracy grows rapidly with c; **an energy-only geometry partition undercounts fragmented states**.

As a deliberately restricted **ground-manifold thermal control**, multiply every graph partition's statistical weight by \(\binom{2c}{c}\), retaining its exact filled-band energy but ignoring all nonzero-energy fermionic excitations. Define \(Q_{n,c}\) by
\[
\sum_{n,c}Q_{n,c}z^nu^c=
\exp\left[\frac u2\sum_{l\ge1}e^{-\beta\delta_l}z^l/l\right].
\]
The recurrence
\[
nQ_{n,c}=\frac12\sum_{l=1}^n e^{-\beta\delta_l}Q_{n-l,c-1}
\]
gives the exact low-energy manifold weight \(Z_n^{\rm zero}=\sum_c\binom{2c}{c}Q_{n,c}\). The one-cycle probability is \(e^{-\beta\delta_n}/(nZ_n^{\rm zero})\), since \(\binom21=2\). Independent direct fermionic incidence-rank controls verify the zero-mode count.

| N | beta=0, ground manifold | beta=1 | beta=5 | beta=20 |
|---:|---:|---:|---:|---:|
| 16 | 0.014170 | 0.096522 | 0.531319 | 0.978969 |
| 32 | 0.004098 | 0.034228 | 0.267228 | 0.820577 |
| 64 | 0.001147 | 0.011030 | 0.112192 | 0.522418 |
| 128 | 0.000314 | 0.003320 | 0.041008 | 0.262486 |

These are **NOT probabilities from the complete canonical fermion Gibbs trace**; nonzero-energy excitations can modify them. The important physical lesson is unambiguous: exactly degenerate matter fillings are another state-counting obligation that a model claiming thermal connectivity cannot hide. The strict finite-N energy selection remains true but need not imply a unique ground *state* or thermodynamic dominance.

The independent verification is [zero_mode_control.py](zero_mode_control.py).

### Robustness check: positive fermion mass

The exact **massless** cotangent theorem does not require retuning any parameters. However CP-11's original fermion mass was m=1, so a separate, explicitly **non-frozen diagnostic** checks whether a positive mass invalidates the loop-energy mechanism.

For an l-cycle define
\[
\delta_l(m)=l\int_{0}^{2\pi}\frac{dq}{2\pi}
\sqrt{m^2+4t^2\sin^2(q/2)}
-\sum_{k=0}^{l-1}\sqrt{m^2+4t^2\sin^2(\pi k/l)}.
\]
A Bernstein/Laplace integral representation of the square root together with the **positive** modified-Bessel coefficients of the ring's heat kernel establishes \(\delta_l(m)>0\) for finite l and m>=0: its Fourier coefficients away from the zeroth mode are negative and the finite ring samples extra negative harmonics. For m>0 the analyticity strip makes \(\delta_l(m)\) **exponentially small** for large l; hence the same fixed-beta cycle-counting argument still drives one-cycle probability to zero in this *conditional* geometry ensemble.

An independent numerical Laplacian-integral control (not a new complete m>0 quantum simulation) gives:

| l | delta_l(m=0) | delta_l(m=0.25) | delta_l(m=1) | delta_l(m=2) |
|---:|---:|---:|---:|---:|
| 2 | .546479 | .359705 | .119152 | .035600 |
| 4 | .264531 | .112694 | .010270 | .000647 |
| 8 | .131237 | .023720 | .00014395 | .000000375 |
| 12 | .087366 | .006605 | .000002449 | .000000000261 |

This is a robustness check of the **finite-size penalty and thermal asymptotic**, not an analytical proof that single-ring minimization under arbitrary m>0 or all fluctuating quantum states follows from the exact m=0 theorem. In fact, gapping the fermions makes the energetic distinction between large competing loops decay even faster.

## 4. Quantum reconnection test through N=36

Supply, without any favored cycle or global connectivity projector, a **positive** geometry reconnection Hamiltonian on perfect matchings:
\[
\boxed{H_{\rm eff}=\sum_G\left(\sum_{\ell\in{\rm cycles}(G)}
\delta_\ell\right)|G\rangle\langle G|
+\frac{\Gamma}{N}\sum_{\{G,G'\}\in\mathcal E}
(|G\rangle-|G'\rangle)(\langle G|-\langle G'|),\quad\Gamma=.1.}
\]
The edges of the *configuration* transition graph are obtained by choosing two existing bonds and replacing them by either alternative pairing of their four endpoints. There are N(N−1) transitions out of every matching. This kinetic term is Hermitian, nonnegative, globally permutation symmetric and O(N)-bounded in operator norm. It can reconnect distant emergent locations: physical locality has not been earned.

We computed the **complete** matching-basis Hamiltonian by independent sparse enumeration at N=3,4,5,6 (dimensions 15,105,945,10395). We also derived its exactly invariant **permutation/port orbit representation**, in which an orbit is specified by the integer partition of N formed by its cycle lengths, and each basis vector is the normalized sum of *all* matchings with those lengths.

Its adjacency matrix is constructed from exact combinatorial split/merge rates and \(\sqrt{s_\lambda/s_\mu}\) orbit-size weights, with explicit symmetry and constant-vector checks. Because this Hamiltonian is stoquastic and its reconnection graph connected, the unique positive ground state belongs to that invariant subspace: its ground energy and geometry-class probabilities are **exactly** recovered by orbit reduction. The complete N=3..6 calculation agrees with the reduced one to ~1e-10, **without excluding physical states before defining the Hamiltonian**. The gap from the reduced spectrum is the *invariant-sector* gap, not necessarily the full first excitation.

| N | Complete matching states | Exact orbit-sector states | Quantum P(one ring) | Mean largest-cycle fraction |
|---:|---:|---:|---:|---:|
| 3 | 15 | 3 | 0.994389 | 0.9981 |
| 4 | 105 | 5 | 0.982037 | 0.9928 |
| 5 | 945 | 7 | 0.959994 | 0.9855 |
| 6 | 10,395 | 11 | 0.929347 | 0.9734 |
| 12 | 316,234,143,225 | 77 | 0.719330 | 0.9034 |
| 20 | >10²³ | 627 | 0.553833 | 0.8543 |
| 24 | >10²⁹ | 1,575 | 0.502504 | 0.8400 |
| 30 | >10⁴⁰ | 5,604 | 0.445824 | 0.8250 |
| 36 | >10⁴⁹ | 17,977 | **0.404312** | **0.8146** |

These probability numbers belong to **the quantum ground of H_eff with Γ=.1**, not the beta>0 Gibbs state. The data show increasing fragmentation even in the quantum zero-temperature state, but the large-N limit of the *quantum* connected probability is NOT proved. It is also meaningful that the largest cycle remains macroscopic at the largest checked N; we should not conflate failing exact connectedness with total absence of organized macroscopic structure.

### Symmetry and operational locality remain open

The ground state is invariant under all simultaneous relabellings. Fixed labelled probes i,j therefore yield the same average response for every distinct pair, unless a physical identifying reference is included. Conditional graph-distance propagation on one ring is not, by itself, an operational separation-dependent law in the fully fluctuating ground state. The actual reference system and local coupling are still unconstructed. Port roles that carry true physical reference records have energetic and accessibility costs.

## 5. Stronger follow-up: complete joint fermions, thermal control and relational response

The initial \`H_eff\` above contains an explicitly adiabatic filled-band approximation. This limitation motivated a separate **fully non-adiabatic microscopic candidate** [FULL_JOINT.md](FULL_JOINT.md) with [full_joint.py](full_joint.py): fixed port-fermion modes (2N total, N occupied), graph-dependent hopping, **all pairings retained**, and the same positive quantum reconnection operator. The relative sign s of within-site and bond hoppings is supplied. The s=-1 sign has the same conditional filled-band energy **without needing to eliminate the fermions**, while the opposite sign s=+1 is a decisive short-loop counterexample.

| N | Complete geometry × fermion states | P(connected), s=-1 | P(connected), s=+1 | Geometry purity at s=-1 |
|---:|---:|---:|---:|---:|
| 3 | 300 | 0.995021 | 0.000004128 | 0.772249 |
| 4 | 7,350 | 0.985453 | 0.000000010 | 0.561881 |
| 5 | 238,140 | 0.968941 | not run | 0.404441 |

These are actual **joint quantum ground states**, not the effective \`H_eff\` predictions; the nonzero geometry–matter entanglement entropies are 0.529, 1.174, 1.869 nats. The exact signed fermion-wedge representation makes the full Hamiltonian commute with a vertex transposition at N=3,4. For s=-1, \(H_{\rm joint}\ge0\) after the constant \(4tN/\pi\) shift.

An exact 300-state **full canonical thermal trace** at N=3 gives P(connected)=0.53333 at beta=0, 0.69210 at beta=1, 0.99702 at beta=5, for s=-1. This is finite-system thermal evidence only; the **full** quantum thermodynamic limit is not proved by the earlier conditional Gibbs no-go, whose scope remains explicit. Reversing sign gives 0.000380 connected at beta=5.

The same ground state has a true *relational detector hierarchy*: source occupancy \(n_0\), one-link detector occupancy at \(\mathrm{mate}_P(0)\), and two-link detector occupation after following two bonds across the intervening vertex. At N=4, t=0.4 the Kubo responses are **-0.266133** and **-0.00381282**, respectively, and agree with independent finite source impulses. At t≈0 the one-link slope is -0.852627 while the two-link linear coefficient vanishes within numerical precision and its first visible term scales approximately as \(-0.0591t^3\). This is stronger than just finding a low energy gap; it links the same evolving connections to actual fermionic detector response. The readout operations remain mathematically bond-controlled and a physical apparatus must still earn access to the mate map.

The crucial **hostile locality control** is [shortcut_control.py](shortcut_control.py). Any two originally distant ring bonds can be reconnected in one quantum geometry operation. A fermion can then hop over the newly created bond in the next operation, giving a second-order matrix amplitude magnitude \(\Gamma t/N\) in \(H^2\) independent of the original ring distance; for N=4 this is exactly **0.025**. Consequently a robust dynamical physical light cone is **not derived**, notwithstanding the relational detector hierarchy and small-N connectedness.

This separate fully joint extension genuinely addresses the original *adiabatic fermion* limitation for N<=5; it does **not** establish a selected 3D spatial phase, a thermally stable complete quantum universe, a dynamically chosen sign/flux sector, local graph reconnection, Einstein's equations, or any empirical novelty.

## 6. What this actually contributes

**Positive:** one connected, gapless 1D geometry is favored by fermion binding across every finite-N pairing in the declared all-pairings state space. Unlike SCOUT-RG01, the cyclic ordering is NOT supplied by the initial kinematics; a single cycle competes with all partitions. Exact conditional energy, a thermal theorem and coherent reconnection ground states were tested, rather than assuming any one from pictures.

**Negative:** two ports per vertex build 1D **valence** into the microscopic law; zero-mode filled-band physics has degeneracy, the strict energetic preference softens as 1/N, fixed-beta conditional Gibbs connectedness vanishes as N^-1/2, quantum effective reconnection appears to fragment the ground with increasing N, labels are not physical localized references, there is no dynamically recovered 3D topology, Lorentz cone, spin-2 sector, universal gravity coupling or experimental prediction. The original `H_eff` fermionic potential is adiabatic; the separate full joint extension removes that shortcut only for the finite sizes tested, and its nonlocal reconnection and full thermodynamic limit remain unresolved.

These negative results are **part of the result**. It would be false to claim we have derived physical spacetime or a thermal quantum gravity theory.

## 7. Next physical obstruction and falsifier

The next successor must NOT fix this by commanding “one component” or by assuming an R³ embedding. It must justify microscopic interaction resources and local reconnection rules which constrain the real accessible matching space, and it must show robust spacetime connectivity at fixed physical temperature or explain its non-equilibrium/cosmological preparation without hiding a temperature that scales with system size. If it changes the number of ports per particle, prove dimensional geometry selection rather than choosing valence to fit a target dimension. An apparently connected small universe at very low temperature is not a valid substitute for those tests.

The strongest new mathematical question suggested by this scout is whether **gauge/charge constraints can make physically accessible reconnections local relative to a geometry generated in the same Hilbert space**, while avoiding the unrestricted graph-entropy no-go. This cannot be asserted from permutation symmetry alone.

## External research

Random pairing/multi-loop counts and \(\theta=1/2\) Ewens asymptotics have substantial prior art, e.g. [Lugo, Profiles of Permutations (2009)](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v16i1r99), [Tavaré, The Magical Ewens Sampling Formula (2021)](https://doi.org/10.1112/blms.12537), and the random configuration-model literature. Coherent quantum loop/reconnection ideas also have prior art, e.g. [Ioffe & Tóth, Split-and-Merge in Stationary Random Stirring (2020)](https://doi.org/10.1007/s10955-020-02487-2). [Konopka, Markopoulou, Severini, Quantum Graphity (2008)](https://doi.org/10.1103/PhysRevD.77.104029) is an essential broader predecessor. The exact free-fermion ring energy is standard spectral mechanics. **Novel mathematical/physical priority has not been assessed adequately to claim originality.**

[check.py](check.py) builds the all-pairing model independently, checks its full quantum reconnection and exact symmetry quotient, tests conditional ring energy and finite-temperature recurrence/asymptotic, and writes results.json. The equivalent standalone algorithms were executed locally. The committed exact script should be rerun verbatim by the partner.
