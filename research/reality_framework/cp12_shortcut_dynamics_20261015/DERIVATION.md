# Physical shortcut operators and the scope of the controls

Parent: `aedf49249d911d938b8c0a0f73dd98f9804879fe`, descending from `69098c127e33ad545044e524fe57cf80339cca08`. All numerical full-law evolutions below use the frozen CP-12 Hamiltonian. A comparison Hamiltonian is always an explicitly identified diagnostic intervention.

## 1. Microscopic hopping, gauge invariance and through-flow

For an active edge e=(i<j), in frozen units t=1,

\[
K_e=n_{g,e}\{(a_i^\dagger U_{ie}-a_j^\dagger U_{je})b_e+
b_e^\dagger(U_{ie}^\dagger a_i-U_{je}^\dagger a_j)\}.
\]

Every term preserves both endpoint Gauss constraints, the edge Gauss constraint and total fermion number. Inactive edges contribute zero. The electric shift by one unit is retained; this is not an adjacency-matrix approximation.

An edge carries a directed through-flow only if E_ie E_je<0, directed from its positive endpoint to its negative endpoint. The diagonal indicator chi_e(s) is one precisely when that arc lies on a directed cycle, calculated by directed reachability. A charged hole may have fields (2,-1); requiring equal and opposite fields would exclude it incorrectly.

Electric fields and chi_e commute with local gauge transformations. They transform covariantly under vertex relabelling. Fermion signs still enter every K_e matrix element; signs are summed before any norm or probability is formed.

## 2. An exact Hermitian moving-core split

Define on the full physical Hilbert space

\[
K_{\rm in}=\frac12\sum_e\{\chi_e,K_e\},\qquad
V_{\rm out}=K-K_{\rm in},\qquad H_{\rm c}=H-V_{\rm out}.
\]

For a physical hop s→u along e, the matrix-element multipliers are

\[
w_{\rm in}=\frac{\chi_e(s)+\chi_e(u)}2,\qquad
w_{\rm out}=1-w_{\rm in}.
\]

They are 0, 1/2 or 1. A hop changing core membership receives half weight on each side. Merely multiplying K_e on one side by a source-state core indicator would generally violate Hermiticity. The symmetric definition does not have that defect.

The split is gauge invariant and permutation covariant, including transitions that change the core. All diagonal electric/mass/graph-cost terms and the full rewiring operator remain in H_c. In particular, its physical rewiring escape diagonal is not recomputed after deleting hopping. No graph projection is introduced.

There are at most 2Q nonzero unit-magnitude hopping transitions from a physical basis state. Multipliers in [0,1] preserve the physical Schur row bound, so ||V_out||≤2Q and ||K_in||≤2Q. Thus these are bounded self-adjoint perturbations of the electric operator; H_c is a well-defined self-adjoint comparison dynamics at each finite N,Q. The signed normalized permutation quotient is a restriction of that operator, not a replacement law.

The original coherent neutral ray psi_0 has no through-flow on any edge. One microscopic hop from it still has no directed cycle, so K_in psi_0=0. Neutral rewiring annihilates psi_0, hence

\[
H_c\psi_0=(Q/2)\psi_0,\qquad
\|V_{\rm out}\psi_0\|^2=\|K\psi_0\|^2=2Q.
\]

The moving-core-only comparison cannot form a core from that state. Its difference from full neutral-origin evolution therefore mixes formation and transport. It cannot by itself diagnose a spatial shortcut. We also start both operators from the **same full-law waiting state**, and from the same explicitly recorded charged-event branch, to test the approximation after circulation has formed.

## 3. A fixed recorded ring and an actual off-ring transfer

N=4,Q=5 is the smallest spanning ring plus one additional active edge. It is an auxiliary resource sector of the same frozen microscopic law; it is not the Q=2N sequence used for N5/N6 formation. Every one of its six labelled hosts is K4 minus one edge and has a unique four-cycle C. The fifth active edge is a diagonal of C. Three hosts put this diagonal between source0 and its opposite receiver; three put it between the other two vertices.

Evolve the original uniform coherent neutral ray over all six hosts under H to time t_w=5. Make the complete twelve-outcome record

\[
\{P_g G_g,\ P_g(1-G_g)\}_{g},
\]

where G_g requires a directed spanning four-cycle, n_v≥1 and a neutral fifth edge (b_e=1,E_ie=E_je=0). All outcomes remain in the complete channel. The selected branches, their unconditional probability, and the energy and entropy change of this external record are reported. No supplied flux or matter is substituted for the waiting state.

The recorded C and the opposite receiver label are fixed during the probe. On its full physical branch Hilbert space define

\[
V_C=\sum_{e\notin C}K_e,\qquad H_C=H-V_C.
\]

This removes hopping along **both possible diagonal labels**, including a diagonal activated later by rewiring. It retains all graphs, all gauge fields, all mass/electric terms and all rewiring. The pointer is fixed; it is not recalculated from the final graph. At most two labels contribute, so ||V_C||≤4 is a conservative physical row bound. The reference graph may change under either dynamics.

For a neutral active fifth edge, its first action has the exact norm

\[
\|K_e v\|^2=\langle v,(2-n_i-n_j)v\rangle,
\quad b_e v=v,\quad E_{ie}v=E_{je}v=0.
\]

There is no zero-field gauge selection rule suppressing this action. Pauli filling can block it when both endpoints are occupied; the actual endpoint occupations must be measured rather than assumed.

For a fixed recorded core C with all active non-C edges neutral in the initial
branch, first-hop outputs on different such edges are orthogonal electric
patterns. Thus the collective identity is

\[
\|V_Cv\|^2=\sum_{e=(i,j)\notin C\text{ active}}
\langle v,(2-n_i-n_j)v\rangle.
\]

This identity assumes that neutral-edge record; it is not asserted for a
general coherent moving-core superposition.

More directly, introduce the gauge-invariant matter transfer

\[
T_{j\leftarrow i}=a_j^\dagger U_{je}U_{ie}^\dagger a_i.
\]

It transfers one vertex fermion i→j, shifts E_ie by -1 and E_je by +1, and leaves b_e unchanged. Write d^dagger=a_i^dagger U_ie-a_j^dagger U_je. The CAR give {d,d^dagger}=2 and

\[
K_e^2=2n_{b,e}+(1-2n_{b,e})
\bigl(n_i+n_j-T_{j\leftarrow i}-T_{i\leftarrow j}\bigr)
\]

on an active edge. Therefore the transfer coefficient is +1 when b_e=1. An occupied source, empty receiver and neutral fifth edge admit one physical two-hop path: b_e→receiver, then source→b_e. Its amplitude in exp(-iH tau) is -tau^2 T/2. Its magnitude is not cancelled by the gauge constraint or Fock signs.

The final fifth-edge fields are (1,-1) or (-1,1); the original ring fields are unchanged. No diagonal or rewiring term can produce these two nonzero half-legs at second order. For a recorded branch X,

\[
\Pr(\text{new fifth-edge through-flow AND original ring at }\tau)
=\frac{\tau^4}{4}\langle X,
[n_i(1-n_j)+n_j(1-n_i)]X\rangle+O(\tau^5).
\]

The coefficient is unnormalized unless divided by ||X||². It is a finite derivative, not a Markov rate or a long-time lower bound. It quantifies allowed charge transfer in either direction. It is **not** the onset coefficient of the phase-response signal below. Under H_C both non-C electric fields stay identically zero, so this joint probability is exactly zero at all times.

The diagonal energy changes along this transfer are also exact. With an
occupied source, empty receiver and neutral filled edge, the first hop creates
one extra vertex particle and an edge hole. Its electric square increases by1,
and the number of removable edges decreases by1. The completed transfer
returns the vertex-particle count, increases electric square by2, and still
blocks removal of that edge. Therefore

\[
\Delta D_1=2+\frac18-\frac{\Gamma(M-Q)}M,\qquad
\Delta D_2=\frac14-\frac{\Gamma(M-Q)}M.
\]

The original core's flux magnitude does not enter either difference, since
its fields are unchanged. These are physical **basis-diagonal differences**,
not spectral gaps or a controlled perturbative effective Hamiltonian. Core
flux protection against rewiring does not itself supply a growing diagonal
detuning for the neutral shortcut. Interference or a spectral suppression
mechanism still requires a calculation in the actual state.

`transfer_algebra.py` checks the CAR identity in all eight local fermion states using integers, and independently enumerates the two-hop Gauss-admissible transfers from recorded physical states, retaining destinations outside its check basis. These are exact algebra checks, separate from the numerical waiting-state expectation.

## 4. Local source, fixed receiver, and record cost

After recording, apply the local gauge-invariant pulse exp(-i theta n_0). The source commutes with all record projectors, so the initial branch probabilities have zero source derivative. For each branch let

\[
x=U(\tau)X,\qquad y=U(\tau)n_0X.
\]

The receiver is the vertex opposite source0 on the **initial** recorded ring. Its density O=n_j remains fixed independently of the final graph. The unnormalized joint susceptibility is

\[
\chi=\left.\frac{d\langle O\rangle_\theta}{d\theta}\right|_0
=2\operatorname{Im}\langle x,Oy\rangle.
\]

Dividing the joint susceptibility by the initial record probability is valid because its derivative vanishes. Its sign is a coherent response to a dimensionless phase pulse; it is not an arrival probability or a transport rate. Both actual ±0.001 pulses and the exact ±pi/2 projector phase identity are evaluated as controls.

The twelve record projectors commute with the unbounded electric/mass terms and all diagonal terms. Consequently A=Phi(H)-H=Phi(K+R_off)-(K+R_off) is bounded:

\[
\|A\|\le2\left(\sqrt{2Q\min(N,Q)}+
\frac{\Gamma Q(M-Q)}M\right).
\]

A waiting-state norm error epsilon bounds the record energy error by ||A||(2epsilon+epsilon²). We never use a bare state-norm error to bound an unbounded energy expectation. For the pure waiting state, the complete recorded-system entropy is the Shannon entropy of twelve orthogonal branches. A conservative total variation bound 2epsilon+epsilon² enters the finite-outcome entropy continuity bound. An apparatus, its preparation and its reset are external resources, not derived dynamics.

The shortcut/no-shortcut categories have source and receiver degree3/degree2 respectively. Comparing those categories is not a degree-matched scaling test. In each full-versus-H_C comparison the initial branch is identical; the intervention isolates off-ring hopping. N4 distances only distinguish one from two active-graph edges and provide no macroscopic distance window.

The separate N5 fixed-receiver response uses the S_(N-1) source stabilizer. The density (n_v-n_0)/(N-1) then equals the expectation of any one fixed receiver j!=0, for both Hamiltonians and the pulse. It does not define a core-distance bin.

## 5. Full-law error accounting

Let S be the retained electric-square-budget space, B=(1-S)HS, and B_c=(1-S)H_cS. Every outgoing destination and amplitude is retained; the electric budget is a numerical approximation, not a modified physical law. For a retained vector v,

\[
\|U_H(t)v-e^{-iSHSt}v\|\le
\int_0^t\|B e^{-iSHSs}v\|\,ds.
\]

The same formula applies to H_c or H_C with its own outgoing matrix. The stream evaluates independent analytic upper sums for each fixed action L. On a cell of width h, midpoint w and A=SHS-energy_shift,

\[
\int_{-h/2}^{h/2}\|Le^{-iAs}w\|\,ds
\le\sum_{k=0}^{5}\frac{h^{k+1}\|LA^kw\|}{2^k(k+1)!}
+\frac{h^7\|L\|\|A^6w\|}{2^6 7!}.
\]

Unitarity makes the sixth moment constant along the retained trajectory. The smaller of this estimate and a first-order-plus-remainder estimate is used. Sparse Schur norms bound the fixed action norms.

For the shortcut Duhamel integral, L stacks retained V and outgoing (1-S)VS. If its retained-trajectory integral upper sum is I_V, the waiting-state error is epsilon_w and the core-trajectory leakage bound at t is epsilon_c(t), then

\[
\|U_H(t)v_{\rm true}-U_c(t)v_{\rm true}\|
\le I_V+\|V\|t[\epsilon_w+\epsilon_c(t)].
\]

For a computed state distance delta_hat an independent interval is delta_hat ±[2epsilon_w+epsilon_H+epsilon_c]. The phase-optimized distance uses the same Lipschitz radius; record branches are phase-optimized separately. Thus a large common dynamical phase is not mistaken for an observable difference.

Record-branch errors combine as a Frobenius norm in the direct sum of orthogonal classical records. The common waiting-state error enters that norm once. If x,y have errors epsilon_x,epsilon_y and ||O||≤1, the susceptibility radius is

\[
2[\epsilon_x\|Oy\|+\epsilon_y\|Ox\|+
\epsilon_x\epsilon_y].
\]

Probability intervals use the contraction amplitude bound (sqrt(p)±epsilon)². Conditional intervals divide joint intervals by a separately bounded positive initial probability. A 5e-10 floating-point reporting allowance is added separately. These are analytic truncation estimates evaluated with ordinary floating point, **not** outward-rounded machine proofs.
