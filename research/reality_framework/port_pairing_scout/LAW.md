# SCOUT-RG02 — permutation-symmetric port pairings and matter-mediated cycle selection

2026-10-11. New **declared test hypothesis**, not CP-12, a law of reality, or a repair to frozen CP-11 / SCOUT-RG01. Forked separately from SCOUT-RG01 at 8813abc. Main GRUT-RAI and Stage-3 Card 3 untouched.

## Kinematics without a supplied circle

N identical vertices each have TWO distinguishable connection ports (a,b). The 2N ports are paired into N undirected bonds, in **every** possible perfect matching: exactly (2N-1)!! orthogonal matching states. The resulting degree-2 multigraph consists of one or more cycles. Self-loops (one vertex) and double edges (two vertices) are **included**; no single connected ring, order, metric, target N-cycle, valence selection or dimension selection is imposed after enumeration.

The two-port constraint, the pair-matching Hilbert-space measure, port distinguishability, ordinary quantum mechanics, external evolution time and fermion species ARE supplied. This is a strong one-dimensional-valence input even though no particular cyclic order is specified.

## Conditional free fermions: the same bonding determines propagation and binding

At each graph G, construct its oriented vertex-to-edge incidence B_G. A self-loop has a zero incidence column. Use **massless** off-diagonal vertex/edge Dirac operator (this differs explicitly from CP-11's frozen mass m=1):
\[
D_G=t\begin{pmatrix}0&B_G\\B_G^\dagger&0\end{pmatrix}.
\]
There are N vertex modes, N edge modes and N fermions. Filling all negative one-particle energies and then zero modes gives the conditional ground energy
\[
U(G)=\kappa N-t\,\mathrm{Tr}\sqrt{B_GB_G^\dagger}
=\kappa N-t\,\mathrm{Tr}\sqrt{L_G}.
\]
N occupied bonds makes the link-cost term common. A connected N-cycle has modes \(\omega_\pm(q)=\pm2t|\sin(q/2)|\), q=2πk/N, so conditional matter propagation has a finite-range graph structure and a gapless long-wavelength band. This **does not** constitute a relativistic gravitational field.

## Explicit geometry reconnection: positive quantum kinetic operator

For each matching G and each pair of distinct bonds, replace the two bonds by either of the two other possible matchings of their four endpoints, obtaining G' (two alternatives). These are all physically allowed reconnection moves. Let \(\mathcal E\) contain each undirected neighboring-pair transition once. Supply
\[
H_{\rm eff}=\sum_G[U(G)-N(\kappa-4t/\pi)]|G\rangle\langle G|
+\frac{\Gamma}{N}\sum_{\{G,G'\}\in\mathcal E}
(|G\rangle-|G'\rangle)(\langle G|-\langle G'|).
\]
Set t=1 and Gamma=.1 as a fixed scout benchmark, not an optimized physical coupling. Every matching has N(N-1) reconnection neighbors counting each move; the reconnection term is positive and has operator norm <=2 Gamma (N-1).

IMPORTANT: This is a **Born–Oppenheimer / integrated-out filled-band EFFECTIVE** geometry quantum Hamiltonian, *not* a complete joint graph+fermion quantum Hamiltonian including nonadiabatic fermionic backreaction. The conditional band-minimizing energy is exact, but replacing the original matter Hilbert space by U(G) is an additional approximation/assumption. Reconnection itself is nonlocal in the emergent graph and requires global access to port identities.

## Scientific tests

1. Prove exact finite-N cycle energy with no preferred graph: for a cycle of length l, U_l = kappa*l - 2t*cot(pi/(2l)).
2. Establish whether the STATIC full ensemble's minimum is one connected N-cycle, then price the gap to fragmentation at large N.
3. Enumerate the exact port-matching measure; derive the finite-temperature distribution of cycle partitions and connected probability for every fixed beta>0. No temperature tuning.
4. Build the full quantum reconnection H_eff for N=3..6 and an **exact permutation/orbit isometry** for larger N; compare energies and connectedness, not only an invariant-sector eigenvalue.
5. Distinguish all-matching support, giant components, the exactly one-ring probability, and operationally labelled vs relational matter propagation.
6. Preserve failures; do not equate a supplied valence-2 topology class with derivation of physical 3+1 spacetime, spin-2 gravity, gauge-invariant local observers or a new empirical prediction.

All calculations are mathematics of an explicitly supplied toy model. Priority over random pairing/Ewens and quantum random-loop statistical physics is NOT claimed.
