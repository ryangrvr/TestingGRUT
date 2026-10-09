# Adaptive law representation and its operational boundary

**DRAFT — pending Claude Code / external scientific review. Issue #5 only.**
Standard mathematical constructions and conventional counterexamples; no
originality claim, candidate K, card, scoring, L0 price, or bank decision.
The proofs are the general evidence. Exact fixtures in [TABLES.md](TABLES.md)
and [RESULTS.json](RESULTS.json) check implementations and examples, not the
universal quantifiers. Source IDs are in [SOURCES.json](SOURCES.json).

## T1. Intervention-preserving finite-memory representation

**Domain and premises.** Let X, L, M, A and Y be nonempty finite sets: a
physical-state label, effective-rule label, internal-memory label, allowed
intervention label, and accessible-record label. No claim is made that these
sets, the memory, time order, or interface emerge. They are the supplied
operational architecture being audited. Put Z = X × L × M. Supply:

1. an initial probability μ on Z and initial readout kernel O₀(y₀|z₀);
2. a time-independent, causal joint conditional probability
   W(z′,y′|z,a), nonnegative with sum over (z′,y′) equal to one;
3. an experiment policy πₙ(a|Hₙ), where
   Hₙ=(y₀,a₀,y₁,…,aₙ₋₁,yₙ) is the accessible history.

The policy is supplied to both models and uses the same histories. It does
not gain access to hidden Z in the comparator. Initial randomization and
policy randomization have the same conditional semantics in both models.
Any finite-window rule can be put in this form by storing its required
window in M, **provided the update itself is stationary**. A finite-state
controller need not store the literal window; a sufficient internal memory
state is enough. Correlations between physical, rule, memory and record
updates are allowed: W need not factor into separate kernels.

**Theorem.** Treat z as the state of a fixed controlled Markov system, with
the single joint transition/record kernel W. For every initial μ, every
allowed policy π, every finite horizon N, and every measurable statistic
of H_N, the original adaptive description and this fixed description give
identical probabilities.

**Proof.** The probability of a finite labelled record is, in either model,

$$
\sum_{z_0,\ldots,z_N}\mu(z_0)O_0(y_0|z_0)
\prod_{n=0}^{N-1}\pi_n(a_n|H_n)
W(z_{n+1},y_{n+1}|z_n,a_n).
$$

These are the same summands on the same domain. Adaptive choice of aₙ
changes neither equality nor the induction: conditional on each common
history, the action distributions are equal. Summing over any event or
bounded statistic preserves equality. This is a single comparator for
**all** policies and initial states, not a new fitted hidden model per
experiment. Joint W also preserves causal reset interventions if resets
are actually included in A with the same update. Equality for resets not
specified in A is not asserted. ∎

**Generalization.** Replace finite sets by standard Borel spaces, sums by
integrals, and W, O₀, π by measurable probability kernels. Successive kernel
composition gives identical finite-horizon laws and a consistent path
measure. No finite dimension, finite precision, computability, locality,
energy conservation, or finite hidden-state capacity follows. The elementary
finite theorem above needs none of these analytic qualifications.

**Deterministic version.** For supplied functions
x′=f(x,l,m,a), (l′,m′)=g(x,l,m,a,x′), and y′=h(z,z′,a),
the fixed map is z′=(f(z,a),g(z,a,f(z,a))) with the same h. Stochastic
updates use their supplied joint kernel. There is no requirement to
preload future interventions or random outcomes into an initial tape.

**Freedom accounting at architecture level.** On a fixed finite Z, A, Y,
unrestricted joint W has |A||Z|(|Z||Y|−1) real simplex dimensions. An
arbitrary adaptive description using W has exactly the same freedom.
Thus the labels “law” and “memory” alone remove **zero** kernel degrees
of freedom. With a fixed deterministic state readout, unrestricted state
transition matrices instead have |A||Z|(|Z|−1) dimensions. A specified
selector or deterministic memory update can reduce these spaces, but its
restrictions remain restrictions on an ordinary fixed kernel. This is not
a charter L0 calculation, nor a claim about physical baseline dimensions.

**Novelty.** State augmentation and controlled hidden-state representations
are established. Predictive-state representations give an alternative
operational description [M03,M04]; adaptively changing network topology is
also part of a joint state [M09]. The explicit policy quantifier here makes
the comparator boundary clear; it is not new physics.

## T2. No architectural discriminator against the unrestricted lift

Let P^π_A and P^π_F be the accessible-record laws furnished by T1. Then
their total variation distance is zero for every π and finite N. For every
randomized decision test φ(H_N) in [0,1],

$$\mathbb E_A^\pi[\varphi]=\mathbb E_F^\pi[\varphi].$$

If its false-positive probability on this fixed comparator is α, its
power on the corresponding adaptive model is also α. No sample size
changes that equality. The adaptive operational class is contained in the
unrestricted fixed-state class; hence no test can reject that **whole**
fixed class while recognizing its embedded adaptive member from these
records. Equality extends to events on countably infinite records if the
path sigma-algebra is generated by finite cylinders.

This is a representation/identifiability theorem. It does **not** prove that
all adaptive proposals are unphysical, that all fixed laws are equivalent,
or that an arbitrary positive kernel has a Hamiltonian, local relativistic,
thermal, or covariant realization. Constraints such as locality, conservation,
KMS/FDT or a hidden-state bound must be imposed separately. A proposal that
violates or restricts such a specified physical comparator remains testable.
“Fixed law” in T1/T2 means the stated abstract controlled dynamical class,
not automatic membership in every conventional physical theory.

## T3. Finite memory is not finite state; infinite memory is not an escape

**Correct finite-state claim.** If all X, L, M are finite and W is
stationary, the lift uses at most |X||L||M| hidden states. A window of r
symbols from finite alphabet B can be stored with at most |B|^r states,
plus a finite initialization convention. Finite numbers of real-valued
coordinates or integer counters do not imply finite state capacity.

**Counterexample to the unqualified claim.** The deterministic sequence
yₙ=1 exactly when n=2^k, k≥0, n≥1, is memoryless given the external time n.
It has a fixed autonomous lift n′=n+1 with output 1 at powers of two, but
cannot be generated forever by an autonomous finite-state deterministic
machine with fixed state readout. Any such machine eventually repeats a
state and then is eventually periodic. This sequence is not eventually
periodic: if it had period p>0 beyond n₀, choose 2^k>max(p,n₀); periodicity
would give a one at 2^k+p, strictly between successive powers of two. ∎

The same impossibility holds for a finite time-homogeneous hidden Markov
machine producing this **deterministic output law**. Its positive-probability
state support evolves by Sₙ₊₁=Post(Sₙ). There are finitely many subsets, so
supports eventually repeat periodically. Every state in each support must
have the required certain output. Therefore that output is eventually
periodic. Stochastic emissions must also be certain on these supports;
edge-emission models can be given a finite state-output lift. ∎

This is a resource obstruction to a **finite** lift, not to a fixed law.
An arbitrary nonperiodic schedule requires an unbounded clock or other
unbounded resource; a periodic schedule only needs its finite phase.

**History lift.** Suppose a fully specified causal model gives measurable
next-step kernels Kₙ(·|h,a) on finite internal histories h. Take the state
space to be the disjoint union of all finite histories, including length
and current rule descriptions. One fixed kernel appends a successor using
K_|h|. It reproduces the model by construction. The union is standard Borel
for standard Borel record spaces. Thus infinite history dependence does
not by itself evade a fixed meta-law; it increases state/resource cost.
Growing graphs or rule descriptions can likewise use a disjoint union of
all admissible finite graphs/programs under their supplied grammar.

Computable self-modification is a fixed interpreter acting on program and
data state, provided an effective update syntax and interpreter are given.
An uncomputable kernel escapes a **computable** comparator, not the class
of all mathematical fixed kernels. Nonlocality escapes a **local** comparator,
not fixed-law representation. An unspecified rule for genuinely novel
events does not define their probability law, so is not a completed
falsifiable alternative. The theorem does not rule out an as-yet unspecified
ontology. [M06,M07]

## C02–C03. A complete finite adaptive-rule countermodel

Supply x,l,m,a in {0,1}, independent noise b with P(b=1)=1/4, and

$$x'=x\mathbin\oplus((1-l)a+lm)\mathbin\oplus b,
\quad l'=l\mathbin\oplus(ax'),\quad m'=x,\quad y'=x'.$$

Here ⊕ is binary addition modulo two; ((1−l)a+lm) is a bit. The rule
label changes and memory affects response. Encode z by 4x+2l+m. The
independent sixteen-row transition table in controls.py is the conventional
fixed comparator, with output z′ div 4. For a=1 at (0,0,0), the readout
is one with probability 3/4, a hand-checkable exact control.

All 128 deterministic binary feedback trees of depth three and all eight
initial states are checked: 1,024 policy/state comparisons and 8,192 record
cells, with exact difference zero. Twenty-four additional randomized
policy/state comparisons run to depth four. A deliberately wrong table is
detected. The proof T1 establishes the statement for every horizon.

At visible x=0 and rule label l=1, a=0 gives P(y′=1)=1/4 for m=0 and
3/4 for m=1. Identical visible reset does not reset hidden memory. That
witness is genuine observable history dependence, fully reproduced by
the fixed eight-state comparator. The control supplies its partition,
memory, discrete order, noise and readout; it derives none as fundamental.

## T4/C04. A real but restricted, measurable exclusion

For a classical hidden-state comparator with at most N states, fix
preparation protocols p_i and common subsequent probe protocol. Let
H_ij=P(response word r_j | do(p_i), do(probe)). If preparation leaves
distribution A_i over the N hidden states and the common probe has response
probabilities B_sj from state s, then H=AB, so rank(H)≤N. This uses the
same state meanings and transition/readout rules across preparations;
protocol-dependent hidden readouts or unfrozen machinery are outside that
comparator. No stationarity across different lengths is needed if the
common probe is at the same elapsed time and any clock is included in the
N-state bound. [M03]

**Countermodel and test.** A two-bit register has state (m₁,m₂), intervention
a∈{0,1}, next state (m₂,a), and output y′=m₁. Starting at 00, prepare
with each two-action word 00,01,10,11, discard preparation outputs, then
probe with actions 00 and record the next two outputs. The response word
is exactly the preparation word. Hence H=I₄, determinant 1, rank 4.
This excludes every N≤3 classical hidden-state comparator under the stated
operational definitions. A fixed four-state machine realizes it exactly.

For finite data, suppose simultaneous confidence enclosures guarantee
|H_ij−Ĥ_ij|≤δ. The spectral norm error is at most 4δ. Weyl's inequality
then gives σ_min(H)≥σ_min(Ĥ)−4δ. A positive lower bound rejects rank≤3
at the enclosure's coverage level. Without a statistically valid enclosure,
a nonzero sample determinant is not a certified rank witness.

This is established predictive-state/linear-algebra structure, not a new
selector law. The independently justified resource bound is what enables
exclusion. Quantum comparators are different: for a d-dimensional retained
quantum state and common instruments, linear response gives rank(H)≤d²,
so rank four does **not** exclude a qubit solely by this bound. Instrument-
specific quantum memory requires additional care [M05].

## C05. The examined law/state-unification matrix model has a fixed lift

Smolin's published model [M01, §2, Eq. (3)] uses antisymmetric integer
N×N matrices and the recurrence

$$X_n=2X_{n-1}-X_{n-2}+[X_{n-1},X_{n-2}],\qquad[A,B]=AB-BA.$$

Put Sₙ=(Xₙ₋₁,Xₙ). The entire recurrence is the fixed polynomial map
(A,B)↦(B,2B−A+[B,A]). For each supplied pair X₀,X₁, existence and
uniqueness at every finite step follow by recursion, without fitting.
Antisymmetry and integer entries are preserved. Matrix pairs can be
unbounded; this is finite coordinate, not finite-state evolution.

For any antisymmetric A, X₀=X₁=A is an exact fixed family. If X₀ and
X₁ commute, Xₙ=X₀+n(X₁−X₀) solves the recurrence, since all subsequent
matrices remain in their commuting linear span. Thus this recurrence
alone does not uniquely select an effective Hamiltonian, constants, or
physical observer. The supplied initial matrix data and the definition of
slow/fast effective variables matter. The noncommuting N=3 fixture checks
four subsequent steps and trace zero exactly. This audit does not rerun
the paper's large-N timescale estimates or dispute its effective-law
interpretation; it shows that this explicit version does not escape a
fixed joint dynamics. The integer matrices and antisymmetric differences
are not automatically quantum density matrices.

## T5/C08–C11. Selection exists, but its target and scope have premises

**Contraction.** On the complete real line let F_a(c)=a+(c−a)/2. Then
cₙ=a+2^-n(c₀−a) for every initial real c₀. This really selects a unique
limit independently of c₀. It selects the **supplied a**, not a derived
fine-structure constant. Replacing a changes the limit. In general a
contraction on a supplied complete domain has a unique fixed point,
but its map and domain carry the selection information.

**Basins.** F(c)=c+c(1−c²)/4 maps [-1,1] into itself: its derivative
5/4−3c²/4 is positive there, with F(±1)=±1. For 0<c<1, c<F(c)<1;
monotone bounded iteration has a positive fixed limit, necessarily 1.
For −1<c<0, the limit is −1 by odd symmetry. Zero remains zero.
Stable attractors do not imply a unique initial-condition-independent
constant. Their scales ±1 are part of the given equation and coordinate.

**Discrete selection.** For two types with supplied positive relative
fitness w, p′=wp/(1−p+wp) has the exact solution
pₙ=w^n p₀/(1−p₀+w^n p₀). If w>1 and p₀>0, the limit is one; if p₀=0,
it remains zero. At w=1 every p is fixed. Choosing which label has larger
fitness and presuming it is available supplies information. Cosmological
replication additionally needs a reproduction mechanism, mutation kernel
and measure [M08,M10]; these controls do not establish that mechanism.

**Stochastic parameter.** A two-label Markov selector with transitions
0→1 of probability q>0 and 1→0 of probability r>0 has stationary
P(1)=q/(q+r), because (1−p)q=pr. Detailed balance holds. Change r
and the concentration changes. This is a conventional fixed dynamics,
not an exact determination of any dimensionless physical constant.

## T6/C12. Historical reinforcement can be a static latent mixture

Consider the standard binary Pólya urn, with supplied positive integer
pseudocounts α,β. After n draws with k successes the next-success
probability is (α+k)/(α+β+n); a realized draw increments its own count.
For any specified word of length n with k successes, multiplying its
conditional probabilities yields

$$P(w)=\frac{(\alpha)_k(\beta)_{n-k}}{(\alpha+\beta)_n},$$

where (u)_j=u(u+1)…(u+j−1), (u)_0=1. Draw instead a **fixed** latent
parameter p from Beta(α,β), then draw independent Bernoulli(p) records.
Integrating p^k(1−p)^(n−k) against the beta density gives exactly the
same expression, using B(α+k,β+n−k)/B(α,β). This proves passive-record
equivalence for every finite word, not just a matched moment.

At α=β=1 the word probability is k!(n−k)!/(n+1)! and pair covariance
is Var(p)=1/12. Equality on every cylinder gives the same infinite record
measure. Conditional on p, the Bernoulli strong law gives empirical
frequency tending to p almost surely; integrating over p proves that
the reinforced frequency has a random beta-distributed limit, not a
uniquely selected constant. Exact code checks all 126 words
of lengths one through six. This is a conventional exchangeability
control [M19], **not** an implementation of Smolin's incompletely specified
novelty/precedence physics [M02].

**Intervention boundary.** Define an allowed intervention to force the
first record to one, and in the urn explicitly increment its success
counter. The next success probability becomes 2/3. In the static latent
model define that force as overriding the emitted record without changing
p; its next success probability remains 1/2. Both would give 2/3 after
**observing** an unforced success. Thus equality of passive record laws
does not imply equality of counterfactual response. These two precisely
defined mechanisms are separable. The fixed *count-state* urn, however,
uses the same forcing/increment intervention and remains exactly
equivalent to reinforcement under T1. The comparator intervention must
be specified before using this contrast; conditioning is not forcing.

## T7/C13. A first event is not structure-free creation

The empty set admits no normalized probability: a measure has μ(∅)=0,
whereas normalization would require μ(∅)=1. There is no function from a
nonempty domain to ∅. A stochastic process cannot start with a normalized
state on literal ∅ and then generate its first event by a defined kernel.
An initial **empty graph** is different: it is one element of a supplied
graph-state space. A birth kernel on {0,1}, with 0→1 probability q and
1 absorbing, gives P(born by n)=1−(1−q)^n. It requires that domain,
discrete order, probabilities, update and initial state. At q=0 there is
no birth. The computation gives 37/64 for q=1/4,n=3.

A calculus of distinction supplies syntax and operations; a monoidal
category supplies objects, arrows and composition; a causal set supplies
order axioms; graph rewriting supplies a grammar and rules; a quantum
cellular automaton supplies cells, algebras and locality/unitarity; a GPT
supplies states/effects/composition. Calling the initial object “zero”
removes none of those premises. This is an obstruction to **that literal
empty-domain probability construction**, not a theorem prohibiting every
cosmological origin model. [M11,M12,M13]

## T8/C14. “Best” needs an independently specified objective

Optimization over branches b requires a supplied feasible domain D and
cost J:D→R, or a supplied ordering of alternatives. On D={0,1,2},
J_A(b)=b² uniquely selects zero; J_B(b)=(b−2)² uniquely selects two.
Nothing about the unchanged branches chooses between these costs.
Even a stationary-action principle normally restricts trajectories given
an action and boundary data; it need not uniquely minimize or identify
physical constants. Free-energy objectives likewise require an ensemble,
temperature/units, model and constraints. Self-consistency c=F(c) is a
condition, not a proof of existence or uniqueness: F(c)=c admits every c,
F(c)=c+1 admits none, F_a above admits exactly its supplied a.

This is not an impossibility result for explanatory objectives. It shows
where their physical commitments and selection obligations reside.
