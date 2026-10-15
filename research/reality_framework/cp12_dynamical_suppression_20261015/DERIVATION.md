# Full-law dynamical suppression diagnostics

This packet continues `768c46c96eae404a527d9cfee8f1b09abeb61a13`.
It keeps the frozen CP-12 parameters, original neutral preparation and
unregulated Hamiltonian. New finite bases approximate that same operator.
The calculations below concern a finite quantum model; no continuum metric,
large-size phase or complete theory of reality is inferred.

## 1. Preparation, event and the meaning of the ratio

Let $M=N(N-1)/2$, with $Q$ active graph edges. A physical basis state is

$$
s=(g,b,E),\qquad
E_{ie}+E_{je}=1-n_{b_e},\qquad
\sum_{e\ni i}E_{ie}=n_i.
$$

Inactive edges have occupied edge fermions and zero electric fields. The
total number of vertex and edge fermions is $M$. With the inherited units,

$$
H=D+K+R_{\rm off},\quad
D(s)=Q/2+2n_v+\tfrac18\sum_\ell E_\ell^2
       +\tfrac1{10M}r(s)(M-Q).
$$

Here $r(s)$ counts active zero-field, occupied-edge-fermion edges. Rewiring
off-diagonal entries are $-1/(10M)$; the full escape diagonal in $D$ is
retained even if a rewiring destination lies outside the numerical basis.

The original ray is the equal coherent neutral sea over **all labelled**
fixed-$Q$ hosts, not a chosen ring or regular host. Its mean energy is $Q/2$,
variance $2Q$, and the rewiring operator annihilates it.

A directed through arc on edge $e=(i,j)$ exists when

$$
E_{ie}E_{je}<0,
$$

and points from the positive to the negative half-leg. Let
$\chi_e(s)$ be one when this arc belongs to a directed cycle, and let
$c=\sum_e\chi_e$. In particular, a hole edge with fields $(2,-1)$
can belong to the circulating core. Requiring equal field magnitudes would
discard valid charged circulation.

The common primary event in the $Q=2N$ calculations is

$$
P=\mathbf1\{c/Q\ge1/4,\ n_v/N\ge1/4\},\quad
w=e^{-iHt_w}\psi_0,\quad X=Pw,\quad p=\|X\|^2.
$$

This defines an explicitly recorded conditional diagnostic:
$\psi_{\rm formed}=X/\sqrt p$. It does not assert that the unmeasured
state autonomously selects a particular geometry. The complete binary
record is $\{P,1-P\}$, not just the successful outcome. Its mean energy change
is

$$
\Delta E_{\rm rec}=-2\operatorname{Re}\langle Pw,H(1-P)w\rangle,
$$

and its outcome entropy is the binary entropy of $p$. Record resources are
diagnostics; they are not inserted into the evolution law.

The primary event commutes with the diagonal and rewiring terms: rewiring
only replaces zero-field edges and cannot change the through arcs or
vertex filling. Thus its record work involves the off-block part of $K$.
For the unweighted $K$, fixing the rotor angles gives a fermion incidence
matrix with Frobenius square $2Q$ and rank at most $\min(N,Q)$.
Its sum of singular values bounds the many-fermion operator norm, giving
$\|K\|\le\sqrt{2Q\min(N,Q)}$ uniformly in those angles. Restricting to Gauss
states preserves this bound. The work intervals use a conservative
$2\sqrt{2Q\min(N,Q)}(2\epsilon_w+\epsilon_w^2)$ allowance.

Retain the inherited Hermitian moving-core decomposition,

$$
K_{\rm in}=\frac12\sum_e\{\chi_e,K_e\},\qquad
V_{\rm out}=K-K_{\rm in}.
$$

A physical hop $s\to t$ receives inside weight

$$
w_{\rm in}(s,t)=\tfrac12(\chi_e(s)+\chi_e(t)),\quad
w_{\rm out}=1-w_{\rm in}.
$$

The weights are $0,1/2,1$. Both operators preserve the gauge constraints and
permutation symmetry. The original neutral ray has $K_{\rm in}\psi_0=0$
and $\|V_{\rm out}\psi_0\|^2=2Q$. The present packet evolves **only $H$**.

For unnormalized $X$, define

$$
J_L(T)=\int_0^T\|L e^{-iHs}X\|^2ds,\qquad
I_L=J_L/p,\qquad \eta=J_{V_{\rm out}}/J_{K_{\rm in}}.
$$

Homogeneity eliminates the normalization from the ratio, but both $I_L$
and the formation probability must still be reported. A positive inside
action is not a demonstration of signal propagation: a local phase probe
can be completely ineffective even when this action is large.

The moving core can change during a hop. These actions include creation,
destruction and rearrangement of circulation; they do not isolate only
connections violating a fixed recorded metric. Event survival

$$
\|P e^{-iHT}Pw\|^2/p
$$

also permits different cycles at the two times. It is not the probability
that the same ring or every intermediate-time geometry survived. A fixed
recorded-core response is tested separately in the auxiliary four-site run.

## 2. Separate tests of the three proposed mechanisms

### Pauli availability

On a definite neutral active edge, $E_{ie}=E_{je}=0,n_{b_e}=1$. Orthogonal
final Fock/electric configurations give

$$
\|K_e v\|^2=\langle v,(2-n_i-n_j)v\rangle.
$$

Both filled endpoints block that edge at that instant. A single empty
endpoint permits a hop. The reported availability fraction is

$$
F_{\rm vac}=\frac{\int\langle\sum_{e\ {\rm neutral}}
 (2-n_i-n_j)\rangle ds}
 {2\int\langle\#\{e\ {\rm neutral}\}\rangle ds}.
$$

This is a bounded diagonal vacancy statistic, not a probability of completed
transfer or a coherent hopping rate. It also does not assert that the state
remains Pauli blocked under subsequent core dynamics.

### Physical diagonal detuning

An allowed hop shifts one half-leg by $\sigma=+1$ for $b_e\to a_i$,
and $\sigma=-1$ for $a_i\to b_e$. Its physical diagonal energy difference is

$$
\Delta D=2\sigma+\frac{2E_{ie}\sigma+1}{8}
 +\frac{M-Q}{10M}\,[r(t)-r(s)].
$$

The outside-hop weights used to average this quantity are $w_{\rm out}^2$.
We report the fractions of this reference action with

$$
|\Delta D|\ge4,\quad |\Delta D|\ge8,
$$

and the RMS of $\min(|\Delta D|,8)$. The cap makes the observable bounded
even for omitted, arbitrarily large electric fields. A Hilbert-norm error
does not by itself control an unbounded energy-square moment.

These are basis-diagonal differences, **not** eigenvalue gaps, dressed
couplings, transition rates or a proof that many-step dynamical detuning is
absent. For a neutral edge the first-hop cost is

$$
2+1/8-(M-Q)/(10M),
$$

and the completed transfer changes the diagonal by

$$
1/4-(M-Q)/(10M).
$$

Fields on a surrounding ring cancel from these differences. Increasing that
ring's circulation alone does not force a large cost on this neutral path.

### Interference without modifying the state

For a hopping action $L$, define in the **full labelled physical electric
and Fock basis**

$$
D_L(s)=\sum_t|\langle t|L|s\rangle|^2,\qquad
\mathcal R_L=\frac{\int\|Le^{-iHs}X\|^2ds}
 {\int\langle e^{-iHs}X,D_L e^{-iHs}X\rangle ds}.
$$

The numerator keeps every coherent cross term. The denominator is an
observable reference evaluated on the same pure trajectory; no dephasing
channel is applied. A ratio above one indicates net constructive cross
terms in this integrated action. It does not rule out cancellation in
particular source-to-receiver amplitudes.

Taking the diagonal of $L^\dagger L$ **after** permutation reduction is
incorrect for this purpose: an orbit vector already superposes physical
configurations. The independent $N=4,Q=5$ check finds a discrepancy of 5
for that shortcut. Constructing $D_L$ before reduction instead agrees with
the full labelled calculation to about $6\times10^{-15}$.

### An exact restriction on interference at a neutral exterior

Fix a recorded core edge set $C$. Let $S_C$ project onto states in which
every edge outside $C$ has zero fields and an occupied edge fermion, and
let $W_C=\sum_{e\notin C}K_e$. Then

$$
S_C W_C^\dagger W_C S_C
 =S_C\left[\sum_{e=(i,j)\notin C}n_{g_e}(2-n_i-n_j)\right]S_C.
$$

Each first exterior hop creates a hole on its own edge and a unit field at
its own endpoint. Those electric/Fock tags make different edges and
endpoints orthogonal. Within one tagged output the parent state is uniquely
reconstructible. Equivalently, the two-hop charge-transfer terms in
$K_e^2$ leave the zero-field sector and vanish between the two $S_C$
projectors; terms with different exterior edges leave an unmatched hole.

Thus arbitrary coherent amplitudes **inside this specified sector** cannot
cancel the first neutral-exterior hopping action. Zero action requires both
endpoints filled on every active exterior edge in the state's support.
It need not require every vertex filled if the exterior graph does not
touch every vertex.

If $C$ is also the complete initial through-cycle core, $W_C\psi$ has an
exterior hole while the other terms in $V_{\rm out}\psi$ do not. Therefore

$$
\|V_{\rm out}\psi\|^2\ge\|W_C\psi\|^2
 =\left\langle\sum_{e\notin C}n_{g_e}(2-n_i-n_j)\right\rangle.
$$

For the ring-plus-chord control, exterior degree is two at every vertex,
so the right-hand side is $2(N-\langle n_v\rangle)$.
The neutral_exterior.py calculation checks the identity on complex
superpositions of all $2^N$ occupation patterns for $N=5,\ldots,8$.
These are supplied validation vectors, not new formation results.

This is an instantaneous sector identity. It neither forbids later
interference after exterior fields develop nor excludes dressed spectral
detuning. Applying it to the unrestricted primary event without checking
its neutral-exterior support would be invalid.

## 3. Exact Pauli controls and a volume loophole in eta

These two families are supplied Gauss basis states used for algebraic
controls. They are not substituted for the neutral-origin preparation.

For $N\ge5$, take a spanning directed ring plus neutral chords joining

$$
i\leftrightarrow i+r\pmod N,\qquad r=\lfloor(N-1)/2\rfloor.
$$

The two edge sets are disjoint and each vertex has two ring and two chord
neighbours, giving $Q=2N$. Put one unit of circulation on the ring, with an
edge hole emitted at each outgoing tail: fields $(2,-1)$, with the positive
field at that tail. Every vertex is filled, every ring edge fermion is empty,
and every non-core active edge is neutral. Gauss holds at each vertex and
edge. Call the resulting basis vector $\phi_N$.

All neutral chord hops are initially Pauli blocked, but

$$
V_{\rm out}\phi_N=0,\quad
\|K_{\rm in}\phi_N\|^2=2N,\quad
\|V_{\rm out}H\phi_N\|^2=4N.
$$

There are two allowed core recombination hops per ring edge. Each creates
one vertex vacancy and leaves the directed ring present. Either of the two
neutral chord edges at that vacancy can then inject its edge fermion. The
resulting $4N$ two-hop targets are distinct, with unit amplitude. Diagonal
terms and graph rewiring cannot contribute to $V_{\rm out}H\phi_N$, since
they leave every vertex filled. Thus the blocked subspace is not invariant
under the full law.

Because $H,L$ and the supplied vector are real, these squared actions are
even in time. At each fixed $N$,

$$
J_{V_{\rm out}}=\frac{4N}{3}T^3+O(T^5),\quad
J_{K_{\rm in}}=2NT+O(T^3),\quad
\boxed{\eta_N(T)=\frac23T^2+O(T^4).}
$$

The leading ratio has no improvement with $N$. No uniform large-$N$
remainder bound or fixed-time limit is asserted. `pauli_algebra.py` checks
the identities in exact rational arithmetic for $N=5,\ldots,12$, including
the full rewiring and escape diagonal.

Now omit the outgoing hole at vertex 0 and leave that vertex empty. The
original ring still circulates. For this one-vacancy state $v_N$,

$$
\|V_{\rm out}v_N\|^2=2,\quad
\|K_{\rm in}v_N\|^2=2N-2,\quad
\lim_{T\to0}\eta_N(T)=\frac1{N-1}.
$$

Nevertheless the neutral chord $0\leftrightarrow r$, spanning initial
ring distance $r$, transfers a fermion from vertex $r$ into the vacancy
at 0 in two hops. The original ring fields are unchanged. The new chord
fields are $(+1,-1)$, and its edge fermion returns to its initial occupancy.
The unique two-hop coefficient has magnitude one, so

$$
\Pr(\hbox{this specified transfer target at }T)
 =\frac14 T^4+O(T^6)
$$

for each fixed $N$. The leading coefficient is independent of $N$, while
the bypassed ring distance grows. Here a decreasing global action ratio
comes from its extensive core denominator. Reporting both integrals and
nonzero core action does not close this loophole: a source-resolved
long-distance response test remains necessary.

## 4. A genuinely neutral-origin Pauli control

The five-site exploratory event $P_F$ requires all vertices filled, one
spanning charged five-cycle, total electric square 25, and every non-core
edge neutral. In the exact signed $S_5$-invariant sector its range is one
dimensional: a single orbit of size 24. Therefore

$$
P_F e^{-iHt_w}\psi_0=\alpha(t_w)\phi_F.
$$

The normalized conditional shape is exactly $\phi_F$, up to an irrelevant
global phase, whenever the formation probability is positive. Waiting-time
error belongs to $|\alpha|^2$; it does not create uncertainty in this
one-dimensional shape. This is why `rare_pauli.py` can use zero initial
shape error while still reporting a nonzero interval for the tiny original
formation probability and the complete binary record resources.

Its conditional vector satisfies

$$
\|V_{\rm out}\phi_F\|^2=0,\quad
\|K_{\rm in}\phi_F\|^2=10,\quad
\|V_{\rm out}H\phi_F\|^2=20.
$$

It exhibits the same $2T^2/3$ onset under the full law. Every vertex is
filled, so $n_i\phi_F=\phi_F$. An immediate source pulse

$$
e^{-i\theta n_i}\phi_F=e^{-i\theta}\phi_F
$$

is only a global phase. Every later density receiver has exactly zero
susceptibility to this particular pulse. This does not exclude every other
type of excitation or later probe, nor does it say the state is static.
It demonstrates why inside hopping action cannot replace an actual signal
test. Its tiny formation weight also fails the appreciable-formation gate.

## 5. Adaptive boundary reduction with the original Hamiltonian

Let $S$ be any retained set of complete signed permutation orbits. With
the exact inherited symmetry embedding, form

$$
H_S=P_SHP_S,\quad B_S=(1-P_S)HP_S.
$$

Every outgoing row of $B_S$ is generated by the unregulated physical
destination routine. There is no electric cutoff on the residual. For
$\widehat v(t)=e^{-iH_St}v$, Duhamel gives

$$
\|e^{-iHt}v-\widehat v(t)\|
 \le\int_0^t\|B_S\widehat v(s)\|ds=\epsilon_S(t).
$$

Adding selected outgoing states changes the approximation space only. It
does not renormalize the escape diagonal, remove physical hopping, select
a favourable trajectory or change the event's definition. Priority is a
deterministic numerical proxy built from the original neutral trajectory and
its charged projection. The selection itself has no physical probability
interpretation; the outgoing residual, rather than proxy coverage, controls
the error.

The six-site sequence starts from the complete electric-square budget 12,
adds 225,000 priority rows, then the 480,212 remaining original boundary
rows, then **all** outgoing rows of that enlarged basis. Independent raw
generation checks every new-to-old coupling against the stored transpose.

In a state-selected basis a rewiring destination can also be omitted.
Consequently $B_S\ne (1-P_S)KP_S$ in general. The action integrals retain
only outgoing entries changing $n_v$ as hopping, and count the remaining
rewiring entries in the **state-error residual**. Confusing these two
operators would inflate the inferred core hopping. The independent small
adaptive-basis check explicitly contains six outgoing rewiring transitions.

## 6. Squared-action quadrature and full-law intervals

Use $A=H_S-Q/2$, a harmless scalar energy shift. For a fixed retained
vector $X$ and complete action $L P_S$, set

$$
f(s)=\|L e^{-iAs}X\|^2,\quad z_k=L A^k e^{-iAs}X.
$$

At a midpoint,

$$
f''=2\|z_1\|^2-2\operatorname{Re}\langle z_0,z_2\rangle,
$$

$$
f^{(4)}=2\operatorname{Re}\langle z_0,z_4\rangle
 -8\operatorname{Re}\langle z_1,z_3\rangle+6\|z_2\|^2.
$$

One cell of width $h$ contributes

$$
\widehat J=h f+h^3f''/24+h^5f^{(4)}/1920.
$$

Let $m_k=\|A^kX\|$, constant under the retained unitary. Then

$$
|f^{(6)}|\le b_L^2\sum_{k=0}^6\binom6k m_km_{6-k},
\quad
|J_{\rm cell}-\widehat J_{\rm cell}|
\le\frac{h^7}{2^6 7!}\sup|f^{(6)}|.
$$

For inside and outside hopping, the full physical row/column bound gives

$$
b_L\le2Q.
$$

A sharper global bound uses the exact Gauss relation between filling and
holes, and applies even when the hop weights depend on the electric state.
Let $F(N,k,e)$ be the maximum number of incidences in a specified $k$-vertex
set among $e$ distinct edges of the complete host:

$$
F(N,k,e)=2\min\{e,\tbinom k2\}
 +\min\{\max(0,e-\tbinom k2),k(N-k)\}.
$$

In a sector with $k$ occupied vertices there are exactly $k$ active edge
holes. Upward hops need an occupied edge fermion and an empty vertex;
downward hops need a hole and an occupied vertex. Their physical row-count
bounds are therefore

$$
u_k=F(N,N-k,Q-k),\qquad d_k=F(N,k,k).
$$

For any Hermitian hopping with physical weights of magnitude at most one,
the block from $k$ to $k+1$ has norm at most
$a_k=\sqrt{u_k d_{k+1}}$ by the rectangular row/column bound. Applying this to
the filling-sector norms majorizes the whole operator by a finite Jacobi
matrix $J$ with off-diagonals $a_k$ and zero diagonal. It has a symmetric
spectrum, so an upper bound on its largest eigenvalue bounds $\|L\|$.

For an integer candidate $b$, the leading principal minors of $bI-J$ obey

$$
D_0=1,\quad D_1=b,\quad D_{k+2}=bD_{k+1}-a_k^2D_k.
$$

The squares $a_k^2$ and all these minors are integers. Strict positivity of
every minor proves $bI-J>0$ without floating-point eigenvalue estimation.
The hopping_bound.py script obtains bounds 8 for $N=4,Q=5$, 14 for
$N=5,Q=10$ and 17 for $N=6,Q=12$. The six-site certificate has

$$
a_k^2=(24,63,96,96,90,60),
$$

$$
(D_0,\ldots,D_7)
=(1,17,265,3434,32938,230282,950374,2339438).
$$

These bounds apply to both moving-core hopping operators on the full
unbounded-electric Hilbert space. The main five- and six-site action
calculations use 14 and 17 respectively; the exploratory rare-state
calculation retains the valid older bound 20. The independent physical
four-site check verifies the underlying hop counts and that the finite
complete-action norms lie below its bound of 8.

For a nonnegative bounded diagonal statistic $D$, apply the same formulas
with $L=\sqrt D$, $b_L^2=\sup_s D(s)$. Both retained and **all outgoing**
action blocks enter $f$; no omitted action is silently set to zero.

The outgoing norm integral for the state error uses midpoint vector Taylor
norms through order five plus an analytic sixth-order remainder, with an
independent first-order bound as a second upper estimate. Let the resulting
monotone bound including waiting error be

$$
e(s)\le\epsilon_w+\epsilon_{S,X}(s).
$$

Right-endpoint upper sums bound $E_2=\int_0^T e(s)^2 ds$. Minkowski in
time-valued Hilbert space gives the useful square-root action inequality

$$
|\sqrt{J_L^{\rm full}}-\sqrt{J_L^{\rm retained}}|
\le b_L\sqrt{E_2}.
$$

Combining it with the quadrature interval $[a,b]$ yields

$$
J_L^{\rm full}\in
\left[(\max(0,\sqrt a-b_L\sqrt{E_2}))^2,
      (\sqrt b+b_L\sqrt{E_2})^2\right].
$$

Formation probability uses the analogous contraction bound

$$
p\in[(\max(0,\sqrt{\widehat p}-\epsilon_w))^2,
       \min(1,(\sqrt{\widehat p}+\epsilon_w)^2)].
$$

Divide intervals only when the denominator's lower bound is positive.
Otherwise the ratio has no resolved finite upper bound. Numerical
intervals include a $5\times10^{-10}$ reporting allowance, but use ordinary
floating point, not outward-rounded interval arithmetic. Independent finite
spectral integrals, full labelled physical calculations and exact rational
controls test the implementation; they do not turn the floating-point
results into machine-verified proofs.

## 7. Actual matter response and its scope

For a recorded branch $X_C$, pulse $e^{-i\theta n_0}$ and a receiver
density $O_C$ fixed from the initial recorded ring,

$$
x_C=U_H(T)X_C,\quad y_C=U_H(T)n_0X_C,\qquad
\chi_C=2\operatorname{Im}\langle x_C,O_C y_C\rangle.
$$

The neighbouring receiver is the mean of the two fixed distance-one vertex
densities. The distance-two receiver is the opposite vertex of that same
ring. Multiplying $O_C$ by the final presence indicator of that **same**
ring yields a joint response with retained initial geometry. This is a
two-time recorded diagnostic, not a condition imposed at every instant.

For state errors $\epsilon_x,\epsilon_y$ and $0\le O_C\le1$,

$$
|\chi_C-\widehat\chi_C|
\le2(\epsilon_x\|O_C\widehat y_C\|
      +\epsilon_y\|O_C\widehat x_C\|
      +\epsilon_x\epsilon_y).
$$

Orthogonal classical record branches combine in Frobenius norm. The full
record's work and entropy are included. The auxiliary $N=4,Q=5$ experiment
does not belong to the $Q=2N$ sequence. Its nonzero neighbour response can
verify that this probe detects matter dynamics, while active diagonal
shortcuts still prevent a claim that the ring determines a local metric.
