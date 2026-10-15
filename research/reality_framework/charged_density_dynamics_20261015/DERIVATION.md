# Charged circulation density: finite-time dynamics and exact local boundary data

Parent: `1e68b1a7d29f3f7b36350aaf7d72700175ad6182`, following the published
memory bridge at `5d63406`. The frozen CP-12 Hamiltonian and original coherent
neutral ray are unchanged. Q=2N except the declared N4,Q4 validation cases.

## 1. Repairing the charge-compatible observable

The edge Gauss law and total fermion number imply

\[
E_{ie}+E_{je}=h_e=1-n_{b_e},\qquad
h_e\in\{0,1\},\qquad\sum_e h_e=n_v.
\]

Equal-and-opposite nonzero fields imply h_e=0. Hence circulation restricted
to such edges satisfies c_bal<=Q-n_v. In particular c_bal=Q and n_v>0 is
empty. This was an observable-definition error in `1e68b1a`, not a physical
obstruction discovered in the Hamiltonian. Its primary definition is repaired
in the accompanying edit, retaining the immutable original commit as history.

For the primary through-flow diagnostic use E_ie*E_je<0. Orient the edge from
its positive-field endpoint to its negative-field endpoint. Writing the fields
in that order, every such edge has the unique form

\[
(q_e+h_e,-q_e),\quad q_e\in\mathbb N,\qquad
E_{ie}^2+E_{je}^2=2q_e^2+(2q_e+1)h_e.
\]

The common through-flow is q_e=min(|E_ie|,|E_je|). This permits charged holes,
including (2,-1). Retain an edge in the circulation subgraph if its directed
endpoints belong to the same nontrivial strongly connected component. The
implementation is checked independently by enumerating every directed simple
cycle in small controls. Labels and an embedding do not enter the definition.

The observables c/Q and n_v/N are diagonal and gauge invariant. We declare
the primary projector before interpreting the results:

\[
P=1\{c\ge\lceil Q/4\rceil,\quad n_v\ge\lceil N/4\rceil\}.
\]

All other declared thresholds are retained as sensitivity controls. There is
no degree-four or connected-host filter. At N5,N6 this event requires c>=3,
n_v>=2. This is weaker than the old all-edge, four-regular target; their
probabilities and onset orders must not be conflated.

For later use, electric L1 satisfies L1>=2c+n_v. The balanced subgraph and
its capacity bound are computed separately, not substituted for charged flow.

## 2. The current is carried by hopping; graph rewiring still affects phases

Neutral rewiring removes/inserts only edges with both fields zero and leaves
the field tuple and all fermion occupations unchanged. Such edges are absent
from the directed through-flow graph. Therefore

\[
[R,c]=[R,n_v]=[R,P]=0.
\]

Electric and mass terms also commute with P. Thus C=PH(I-P)=PK(I-P),
with ||C||<=b=|t_hop|sqrt(2Q min(N,Q)). For Q=2N at frozen hopping one,
b=2N. This does NOT make formation independent of rewiring: R remains in
both diagonal blocks A=PHP and D=(I-P)H(I-P), and changes their phases and
the resulting source/memory evolution.

Writing x=P psi, z=(I-P)psi gives

\[
\dot p=2\operatorname{Im}\langle x,Cz\rangle,
\]
\[
\dot p=2\operatorname{Im}\langle x,Ce^{-iDT}\psi_0\rangle
-2\operatorname{Re}\left\langle x,
\int_0^T Ce^{-iD(T-s)}C^\dagger x(s)ds\right\rangle.
\]

The numerical source is propagated directly in the complement. Return is
computed as total current minus source; it is not an independent memory
quadrature. Central finite differences of p check the total-current sign and
magnitude. No positive, autonomous formation-rate law is assumed or derived.

## 3. Controlled streaming evolution and sharper current uncertainties

The retained basis is complete reachable electric-square budget closure,
with the original rewiring escape diagonal and every outgoing full-law
transition in B_out. The inherited signed permutation quotient is exact.
The C++ acceleration performs only integer permutations and Fock parity;
small complete Hamiltonian and outgoing-Gram comparisons agree exactly with
the legacy Python builder. Additional N6 destination comparisons are retained.

State errors use Duhamel and analytic midpoint Taylor upper sums for the
leakage vector through fifth order with a sixth-order remainder. Streaming
stores short exponential grids rather than an entire trajectory in memory.
It is checked against the inherited full-grid propagator. Optional caches
store exact constructed bases/matrices, not precomputed scientific outputs.

The observable J=i[H,P] has ||J||=||C||<=b. If ||psi-w||<=epsilon, a useful
state-dependent uncertainty is

\[
|\langle\psi,J\psi\rangle-\langle w,Jw\rangle|
\le2\epsilon\|Jw\|+b\epsilon^2.
\]

The omitted part of Jw is not silently dropped. With finite-budget projector
S, its represented parts obey

\[
\|Jw\|^2\le\|J_Sw\|^2+
(\|B_{\rm out}Pw\|+\|B_{\rm out}w\|)^2.
\]

This follows by bounding the outer component
B_out P w - P_out B_out w by the triangle inequality. It requires no guessed
event labels on omitted states.

For source current, put v=e^{-iDT}psi_0, and let v_hat have error eta. Then
the corresponding bilinear uncertainty is at most

\[
2[\epsilon\|C\hat v\|+\eta\|C^\dagger Pw\|+b\epsilon\eta].
\]

The full actions are bounded by their retained norms plus, in quadrature,
||B_out v_hat|| and ||B_out P w|| respectively. The latter are conservative
because the omitted P or complement projectors are contractions. Both sharper
and generic bounds are valid; the code takes their minimum and adds a separate
5e-10 reporting allowance. Return uncertainty is the sum of total and source
uncertainties. These are analytical bounds evaluated in ordinary floating
point, not outward-rounded certificates.

An event probability is enclosed using
max(0,sqrt(p_hat)-epsilon)^2 and min(1,(sqrt(p_hat)+epsilon)^2), including
the reporting allowance. Good norm/energy conservation alone never licenses
a formation claim.

## 4. Exact local onset without building a larger Hilbert space

At N5,N6 the minimum event L1 is eight, requiring one directed triangle and
two holes. At N7,N8 the quarter thresholds require c>=4,n_v>=2, so the minimum
L1 is ten and the circulating core must be one directed square. A simple
directed graph with four cyclic edges cannot instead be a triangle plus one
extra cyclic edge. This exhausts these minimum-order core possibilities.

At the first possible order every operation must increase L1 by one. No
diagonal or rewiring insertion contributes, and unused host edges cannot
participate. Two holes have distinct positive endpoints by the vertex Pauli
constraint. A hole on the cycle raises its positive leg from one to two;
an off-cycle hole contributes (1,0) or (0,1). All such targets and prefixes
have electric square at most 2*ell+6. In the triangle case this is twelve,
which proves the integer boundary calculation is unregulated despite its
finite target-directed support.

Fix one directed ell-cycle and enumerate all two-hole placements eta and their
allowed distinct positive endpoints. Let k_eta be the number of distinct
edges in cycle union holes. Let A_eta be the signed integer matrix element
of K^(2ell+2) into that target. Each target occurs in exactly
binom(M-k_eta,Q-k_eta) host extensions. The number of labelled directed
ell-cycles is (N)_ell/ell. Therefore

\[
\boxed{p_P(T)=\frac{(N)_\ell/\ell}
{\binom MQ[(2\ell+2)!]^2}
\sum_\eta\binom{M-k_\eta}{Q-k_\eta}|A_\eta|^2
T^{4\ell+4}+O(T^{4\ell+6}).}
\]

All interfering hopping paths are summed before squaring. Spectator filled
fermions can change an overall Fock phase but cannot change this squared
amplitude; unused edge hops would violate minimum order. Each final labelled
host has one neutral parent, so no inter-host coherence factor is inserted.
The even remainder follows from real H, real preparation and real event P.

The resulting coefficients are:

| N,Q | First probability power | Exact coefficient |
|---|---:|---:|
| 5,10 | 16 | 10264/6615 |
| 6,12 | 16 | 2228176/601965 |
| 7,14 | 20 | 8852523251/8726130000 |
| 8,16 | 20 | 4312328/2095875 |

N5/N6 agree with both an integer weighted-orbit recurrence and the recovered
integer action of the complete budget matrix. The N5 recurrence also runs
through the independent legacy Python canonicalizer. The local calculation
uses independent signed half-leg target DPs from the previous sprint and
matches the entire N5/N6 squared norms. It evaluates all four cases in well
under one second here. N7/N8 are new exact boundary data, not finite-time
formation estimates, and these formulae do not extrapolate to N>=9 without
classifying its different matter thresholds and minimal cyclic cores.

## 5. Persistence and interference in the same unfiltered state

At wait w decompose psi(w) into P and Z components, then propagate both with
the unchanged full H. Define h=P U(tau)P psi(w), e=P U(tau)Z psi(w). Then

\[
p(w+\tau)=\|h\|^2+\|e\|^2+2\operatorname{Re}\langle h,e\rangle.
\]

The sequential event probability is ||h||^2; its conditional survival is
||h||^2/p(w). Recording the initial P/Z flag deletes the last interference
term in the later probability. It is therefore a different experiment from
unmeasured evolution. The reported history overlap is not a free classical
record, a selected outcome, or survival of one identified graph.

If approximate branch errors are epsilon_h,epsilon_e, the overlap uncertainty
is at most epsilon_h||e_hat||+epsilon_e||h_hat||+epsilon_h epsilon_e.
Conditional probability intervals control numerator and denominator separately.
N5 resolves the interference term; the stated N6 bound does not resolve its
nonzero sign, despite a nonzero central estimate.

## 6. Recording costs without expanding every labelled graph branch

For complete record projectors P_g P and P_g Z, graph pinching removes physical
off-diagonal rewiring. It retains hopping only within the same event flag and
retains the physical rewiring escape diagonal Gamma*r*(M-Q)/M. In a permutation
quotient, some removed rewiring transitions may have folded into a matrix
diagonal, so keeping the numerical H diagonal alone would be incorrect.

Hopping changes n_v by one while rewiring preserves it. Thus the hopping
matrix is identifiable by this exact selection rule; adding the known physical
mass/electric/rewiring diagonal gives the pinched expectation correctly.
Uniform probabilities within each graph-label orbit give record entropy

\[
S_{\rm rec}=-\sum_{g\text{ classes},f}p_{g,f}
\log_2[p_{g,f}/m_g],
\]

where m_g is its labelled-graph multiplicity. For a pure input this equals
the measured-system entropy after the complete orthogonal record channel.
It need not increase total entropy of a unitary system-plus-apparatus.

The energy-change operator Phi(H)-H contains bounded hopping/rewiring only,
because every record projector commutes with the unbounded electric energy.
Use ||Phi(H)-H||<=2[b+Gamma Q(M-Q)/M] for its uncertainty, not a norm bound
on an unbounded energy expectation. Finite-outcome entropy continuity supplies
the entropy estimate. Apparatus preparation/reset costs are not derived.
The resulting N4 energy/entropy match the prior full physical-basis calculation.

After the explicit graph record, the single-host leakage theorem gives
same-host survival >=1-C_N^2 tau^2, C_N=Gamma sqrt(Q(M-Q))/M. It preserves
random hosts too and supplies no spatial selection mechanism.

## 7. Matter response and the remaining geometric obligation

For initial charged component X=P psi(w), source B=n0, and fixed receiver
density O=n_j, the pulse exp(-i theta B) yields the unnormalized susceptibility
2 Im <U X,O U B X>. Its conditional denominator p(w) has zero source derivative
because [B,P]=0. The source-fixing permutation quotient permits evaluation
of any fixed j!=0 as the average over those equivalent vertices. No detector
is reassigned using a final graph distance. Full outgoing errors are retained
for both response columns and the initial wait.

This is a charged-density/matter-response control, not a spatial cone test:
K5 has no receiver beyond distance one, and every N6,Q12 host has diameter
two. Circulating fields record/protect their nonzero edges, but the quarter
event does not fix the other active connections. Before a protected core can
be called an operational geometry, its metric must predict the actual matter
response including possible shortcuts through those other connections.

Quantum Hilbert space, the Born rule, Fermi statistics, compact U(1), external
Hamiltonian time and microscopic constants remain inputs. No continuum,
Lorentz symmetry, Einstein gravity or Standard Model matching is obtained
from these finite calculations.
