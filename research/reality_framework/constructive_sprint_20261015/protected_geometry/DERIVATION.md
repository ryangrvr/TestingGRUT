# Conditional flux lifetime and matter propagation under frozen CP-12

This calculation supplies an initial graph and a flux record. It establishes a
conditional dynamical mechanism, not the autonomous formation of either object.
All couplings retain their frozen values. The simulation includes every allowed
graph, and finite electric boxes are accompanied by full outgoing-error bounds.

## 1. A bounded local hopping term despite unbounded electric fields

For one active connection e=(i,j), write

\[
 K_e=t\{(a_i^\dagger U_{ie}-a_j^\dagger U_{je})b_e+\mathrm{h.c.}\}.
\]

The two rotor shifts are commuting unitaries and commute with fermions. Therefore
the bright mode

\[
 c_e=(U_{ie}^{\dagger}a_i-U_{je}^{\dagger}a_j)/\sqrt2
\]

obeys the canonical anticommutation relations, anticommutes with b_e, and
K_e=√2 t(c_e†b_e+b_e†c_e). The two-mode hopping operator has norm one. Thus

\[
 \boxed{\|K_e\|\le\sqrt2|t|.}
\]

The inactive-edge projector can only reduce this norm. This norm calculation does
not replace an electric shift by a truncated nonunitary shift. `run.py` also
checks an independent three-fermion CAR matrix at several rotor angles.

## 2. A local field-distance bound

Let H_0=H−K_e. Every term in H_0 preserves the pair (E_ie,E_je): neutral rewiring
changes an edge only at zero fields and never changes its field values; other
edge hopping terms change other fields. Let the initial state have

\[
 |E_{ie}|+|E_{je}|\ge d
\]

on every configuration in its support. A K_e event shifts just one half-leg by
one unit. Consequently, at least d events from K_e are needed to reach E_e=(0,0).
H_0 and H are self-adjoint on the electric-energy domain because the remaining
terms are bounded at fixed finite N. The interaction-picture Dyson expansion
through order d−1 has zero overlap with that target. Its unitary integral
remainder obeys

\[
 \boxed{\|P_{E_e=0}e^{-iHT}\psi\|
 \le \frac{(\sqrt2|t|\,|T|)^d}{d!}.}
\]

There is no exponential prefactor: the final exact evolution in the remainder
is unitary. Squaring the norm and capping at one gives a full-law, all-time
inequality, although it is useful only before its right-hand side becomes large.

For a supplied Euler-oriented graph with every edge carrying (±L,∓L), d=2L.
Finite supplied matter dressings which only increase positive half-legs preserve
this conservative distance. The simulation's single charged edge actually has
distance 2L+1, and the other three have distance 2L.

## 3. Final graph persistence

The active-edge count Q is conserved. An admissible absent edge has both fields
zero. If the final graph differs from the supplied Q-edge host, at least one of
its original Q edges is absent. Commuting diagonal projectors therefore give

\[
 \boxed{
 1-\langle P_g\rangle_T
 \le\min\left\{1,
 Q\left[\frac{(\sqrt2|t|\,|T|)^{2L}}{(2L)!}\right]^2\right\}.}
\]

This is a statement about the final graph readout. It is not the probability of
a continuously recorded history in which no rewire occurred. It does not assume
that individually labelled edges can become preferred from a permutation-neutral
preparation. The supplied host explicitly breaks that preparation symmetry.

Stirling's lower bound (2L)!≥(2L/e)^(2L) gives the convenient corollary

\[
1-\langle P_g\rangle_T
\le Q\left(\frac{e|t|\,|T|}{\sqrt2 L}\right)^{4L}.
\]

The guaranteed interval increases linearly with L, while the supplied vacuum
electric energy is ε_E QL². This is a resource-costed protection mechanism,
not free geometry stabilization. It says nothing by itself about the mechanism
that might generate these large fluxes from the neutral input.

## 4. Approximation by a fixed-host local Hamiltonian

Let H_g=P_gHP_g and evolve the same supplied initial state with H_g. Its
diagonal rewiring term is retained, Γ(M−Q)r/M. For each removable edge e and
absent destination f, the output graph g−e+f identifies e and f uniquely. Hence

\[
\|(1-P_g)HP_g v\|^2
=\frac{\Gamma^2}{M^2}(M-Q)\sum_{e\in g}\|r_e v\|^2.
\]

Applying the edge-distance estimate inside H_g, followed by Duhamel, yields

\[
\boxed{
\|e^{-iHT}\psi-e^{-iH_gT}\psi\|
\le\frac{\Gamma}{M}\sqrt{Q(M-Q)}
\frac{(\sqrt2|t|)^{2L}|T|^{2L+1}}{(2L+1)!}.}
\]

This coherent-state norm estimate controls every bounded final observable and
bridges the supplied local-host dynamics to the unchanged full Hamiltonian.
Its root-mean-square factor applies to one fixed host. It is not justified for
an unrestricted coherent superposition of input hosts with overlapping outputs.

There is also a less sharp, more economical bound requiring no initial flux at
all. Since each r_e is a projector, the same identity gives

\[
\boxed{\|e^{-iHT}\psi-e^{-iH_gT}\psi\|
\le C_N|T|,\quad
C_N=\frac{\Gamma}{M}\sqrt{Q(M-Q)}.}
\]

For Q=νN and fixed ν, C_N~Γ√(2ν/N). Thus every supplied sparse host has a
controlled fixed-host approximation for T=o(√N), independently of the electric
resource cost. The frozen 1/M rewiring normalization and the orthogonality of
outgoing graph records produce this result. The global rewiring operator norm
is not small, and coherent multihost initial states are not covered. In
particular, this is a preservation theorem, not a formation theorem. It can be
combined with locality of H_g on any bounded-degree supplied host; it does not
select which host should be supplied.

Within a fixed host, the electric energies and removable-edge penalty act on
individual edge Hilbert spaces. Vertex mass terms act on vertices. Hopping joins
only incident vertex and edge degrees of freedom and has bounded norm independent
of L. Interaction-picture locality estimates can therefore use the incidence
graph without truncating electric fields. A lattice and a preferred dimension
are not derived by this observation: the host remains an explicit input.

## 5. Same-state response experiment

The diagnostic input is a supplied directed four-cycle carrying flux L, with one
vertex-0 fermion made by the gauge-invariant edge-to-vertex hop on edge (0,1).
The matching edge hole and +1 electric shift are included. Its full-law moments
are

\[
 \langle H\rangle=L^2+\frac L4+\frac{33}{8},\qquad
 \sum E^2=8L^2+2L+1.
\]

After unforced waiting time w, apply exp(−iθ n_0), evolve for lag τ, and measure
receiver density n_j, host identity, full protection, and the original directed
circulation on the same evolved preparation. No circulation postselection is
used. With x=U(τ)U(w)ψ and y=U(τ)n_0U(w)ψ,

\[
 \chi_j=\left.\partial_\theta\langle n_j\rangle_\theta\right|_0
 =2\operatorname{Im}\langle x,n_jy\rangle.
\]

The receivers at vertices 1,3 have initial-host distance one; receiver 2 has
distance two. This tiny supplied cycle cannot establish a macroscopic causal
cone, but it directly tests nonzero matter response while the same graph record
is retained with high probability.

For errors ||δx||≤e_x and ||δy||≤e_y, the conservative detector error is

\[
 |\delta\chi_j|\le2(e_x+e_y+e_xe_y).
\]

All finite-box exits are explicitly constructed from an independent full-law
transition implementation. The imported leakage integrator bounds the error by
unitary Duhamel evolution and analytic midpoint remainder sums. Calculations use
ordinary floating-point arithmetic, so these are analytically derived numerical
error estimates, not formal interval-arithmetic certificates. Comparisons at two
electric cutoffs check their practical consistency. The saved probability and
response intervals add a separate numerical reporting allowance of 5×10⁻¹².
This exceeds observed cutoff discrepancies but is not a certified floating-point
roundoff bound; analytic electric-truncation estimates are saved separately.

The exact hopping-off control is χ_j=0: every vertex density commutes with the
diagonal Hamiltonian and neutral rewiring. Thus this source-density response
requires the matter-hopping mechanism.

## 6. Flux does not guarantee efficient transport

At a vacuum unit-occupied edge, creating matter on a half-leg of field ±L costs

\[
 \Delta_\pm=2m+\frac{\epsilon_E}{2}(\pm2L+1)
 =\frac{17}{8}\pm\frac L4.
\]

Thus stronger flux need not monotonically suppress motion. The negative-field
channel is nearly resonant around L=8.5 at the frozen couplings. At much larger L
both channels are strongly detuned. The protection theorem is independent of
these detunings; useful propagation is an additional dynamical requirement, which
the same-state response experiment measures rather than assumes.

## Scope

This construction demonstrates how supplied relational geometry can be protected
and can coexist with matter response under one unchanged physical Hamiltonian.
It does not select a host, a dimension, an initial flux record, an initial state,
a macroscopic metric, Einsteinian gravity, or the observed matter spectrum.

## Exact supplied-square winding certificates

The generic 2L local lower bound need not be tight once the other Gauss equations
are imposed. For the charged square used in the simulation, fixing any original
edge to zero fixes the only cycle-flow variable. The alternating vertex/edge
Gauss constraints then bound every target half-leg by four. Complete fermionic
enumeration yields 93 possible targets, whose actual maximum |E| is two.

For this finite complete target catalog, the minimum electric L1 distance from
the prepared input is

\[
\boxed{d_{\min}(L)=8L-5,\qquad L=1,2,\ldots.}
\]

L=1 is checked directly. For L≥2, every target's distance is an affine function
of L with slope eight because its fields have |E|≤2; checking L=2 therefore
establishes the entire tail. Before the first zero edge, no rewire is allowed.
Consequently, at least 8L−5 hopping events are required to expose a removable
edge, and at least one subsequent rewiring event is required to change the host.

The independent unregulated integer Hamiltonian implementation, retaining signed
fermionic amplitudes, gives the following exact *nonzero* first loss terms:

\[
\begin{aligned}
1-p_g(T)\big|_{L=1}
 &=\frac{T^8}{259200}+O(T^{10}),\\
1-p_g(T)\big|_{L=2}
 &=\frac{216017}{12906142470144000000}T^{24}+O(T^{26}).
\end{aligned}
\]

The exact code also verifies nonzero leading loss at T^40 for L=3 and T^56 for
L=4, with rational coefficients saved in `EXACT_ONSET.json`. These are local
finite-preparation certificates, not long-time Taylor extrapolations or claims
that matter response vanishes until those orders.
