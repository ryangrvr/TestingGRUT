# CP-12 higher-density locality: exact flux/fermion protection tradeoff

**14 October 2026.** Independent analytic extension of 03256f9. Frozen Hamiltonian unchanged. This is a **necessary resource constraint**, not a demonstrated spontaneously formed spatial phase.

## Exact statewise identity

On every physical CP-12 basis state with N vertices and Q present connections, let

- \(n_v=\sum_i n_{a_i}\in[0,N]\) be occupied vertex modes;
- \(r\) be the number of **removable** present edges: edge fermion occupied and both half-leg electric fields zero;
- \(k\) be the number of **occupied-edge-fermion** present edges with nonzero electric flux;
- \(h\) be the number of present edges with **unoccupied** edge fermions.

The absent-edge fermions are filled, and the conserved total fermion count is M=binom(N,2). Therefore
\[
h=n_v,\qquad r+k+h=Q,\qquad k=Q-r-n_v\ge0.
\]
On an occupied-edge-fermion present edge, Gauss law gives \(E_{ie}+E_{je}=0\). If it is nonremovable, integer fields force \((E_{ie},E_{je})=(j,-j)\), \(|j|\ge1\), contributing **at least 2** to \(\sum E^2\). On a hole edge, Gauss law gives \(E_{ie}+E_{je}=1\), so integer fields contribute at least 1 to \(\sum E^2\). Thus the frozen electric Hamiltonian \(F=\frac{\epsilon_E}{2}\sum E^2\) obeys the **diagonal operator inequality**
\[
\boxed{F\ge \epsilon_E k+\frac{\epsilon_E}{2}n_v
=\epsilon_E(Q-r-\tfrac12 n_v)
\ge\epsilon_E(Q-r-\tfrac12 N).}
\]
The positive mass-plus-electric part obeys the **stronger** diagonal inequality (when \(2m\ge\epsilon_E/2\), true for frozen m=1, epsilon_E=1/4):
\[
\boxed{2m n_v+F\ge\epsilon_E(Q-r)
+(2m-\epsilon_E/2)n_v\ge\epsilon_E(Q-r).}
\]
For Q=3N and r=o(N), the electric energy is necessarily at least
\[
F\ge (5/2)\epsilon_E N-o(N)=0.625N-o(N)
\]
and mass plus electric energy at least
\[
2mn_v+F\ge 3\epsilon_E N-o(N)=0.75N-o(N).
\]
The parent paper's \(F\ge(\epsilon_E/2)(Q-r)\) remains valid but is substantially weaker because it counts only one half-leg unit on a nonremovable occupied edge.

These are *lower bounds on positive Hamiltonian components*, **not lower bounds on the total conserved energy**, because hopping can be negative. They do not by themselves forbid a prepared or dynamically formed low-r state.

## Exact dynamical attribution

The rewiring Hamiltonian R acts only by swapping one removable edge for an absent edge that is initially in the same zero-field, occupied-fermion configuration. The swapped edge is removable before the move, and the newly present edge is removable after it. All other edge occupancies and fields are unchanged. Thus r is constant across each off-diagonal rewiring matrix element. Its diagonal Laplacian terms commute with r as well. Mass and electric diagonal terms commute with r. Hence
\[
\boxed{[R,r]=0,\quad [H_{\rm diag},r]=0,\quad
\frac{d}{dt}\langle r\rangle=i\langle[H_{\rm hop},r]\rangle.}
\]
**Only matter hopping can change the total number of removable connections.** Rewiring can change *which* edges are removable, not their number. The source-free neutral sea starts with r=Q, n_v=0, F=0. Any autonomous transition toward r=o(N) requires matter hopping to create a macroscopic amount of electric/fermionic protection; rewiring alone cannot do it.

For a density matrix the diagonal inequalities remain valid on expectations, e.g.
\[
\langle F\rangle\ge\epsilon_E(Q-\langle r\rangle-\tfrac12\langle n_v\rangle).
\]
The commutator identity holds for finite-energy states in the unbounded-rotor model; bounded r avoids domain difficulties under the usual finite-energy dynamics.

## Scientific interpretation and next gate

The 03256f9 state-dependent bound shows r=o(N) is a sufficient condition to suppress coupling to an O(N) candidate shortcut pool at fixed times, **provided the small r persists over the time interval**. The present result shows the required protection is an *extensive gauge-flux/matter reorganization*, and identifies the only term that can generate it: \(H_{\rm hop}\). At Q=3N, occupancy of vertex fermions alone cannot protect most edges. This makes the correct constructive target a **self-generated flux-protected phase**, not merely an arbitrary pre-pinned graph.

Caveat: globally reducing r may freeze *all* rewiring, including potentially desirable local geometry relaxation; no locality, 3D dimension, propagation cone, state selection, or ground-state phase is established. Stronger work must demonstrate **selective** distant shortcut suppression with retained local transport, and account for energy/entropy/preparation. The operator-capacity bound itself is worst-case and may be avoided by destructive interference even with r=Theta(N); r=o(N) is sufficient, not necessary.

## Minimal independently executable identity check

The companion check.py enumerates all N=3,4 physical basis states at finite electric cutoffs using the frozen model, validates the statewise inequalities and verifies that every off-diagonal rewiring transition preserves r. The mathematical proof above does not rely on the cutoff. The checker does not establish a many-body phase.
