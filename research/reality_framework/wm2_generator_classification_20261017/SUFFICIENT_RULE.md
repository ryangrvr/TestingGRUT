# A supplied interaction rule, and its complete-evolution guarantee

**Tier: supplied generator restriction; author-derived conditional theorem,
independent review pending.** This class is a mathematical comparison, not
a retroactive repair of UC1's failed gate or an adopted new architecture.

## 1. Define the restriction before evaluating response

Retain UC1's exact physical kernel, two species, zero electric energy,
storage/receiver law, all zero-flux-gated append/pop terms and their adjoints.
The word and its occupied consecutive bonds remain generated, not a supplied
event-pair graph. Its path dimensional class remains supplied by the grammar.
The following **universal relational interaction rule RI(R,w)** is additional
microscopic data, not implied by those constraints.

* Fix an integer R>=1 independent of N. Every matter generator, for every
  species and every internal operation, is supported on a consecutive occupied
  word interval of diameter <=R. A coefficient may read only that interval's
  matter registers and occupancy/bond data. No distant spectator, root-to-all
  control, global word-dependent coefficient, or off-bond pair gate is allowed.
  New matter coefficients do not read the movable head/free flags. Each
  primitive has an explicit bounded-register implementation commuting with the
  supplied constraints; an arbitrary global P h P is not an implementation.
* Require separate species-number conservation. On a fixed occupied interval,
  the physical tree-kernel isometry represents the generator as an ordinary
  bounded local operator on that interval's matter site factors. Position
  transfers have their explicit finite-interval Wilson dressing. Internal,
  density and species-converting-level operations are included, not exempted.
  The class here excludes arbitrary additional electric or storage couplings;
  zero electric energy and the original storage law are retained.
* Each such local generator annihilates the interval matter vacuum. Pure
  vacuum constants can be removed, but geometry-dependent vacuum energies
  cannot be silently discarded. All new terms vanish in full matter vacuum.
* After taking out a sum of on-event, position-preserving internal energies,
  the total matter generator has norm <=w<=1/2 on the declared zero/one-per-
  species sectors, on every word and on every reflecting prefix restriction.
  Require an ordinary local number-preserving site-algebra extension for
  the commutator support calculation. This is a few-carrier norm condition,
  not a finite-density assertion. Extensive systems need a connected local
  interaction norm and a separate theorem.
* No new direct matter/receiver or distant matter/storage interaction is
  admitted. The only receiver contact and geometry changes are UC1's original
  movable-end primitives; their internal-state action is the identity. The
  original particle-occupation gates remain. Arbitrary diagonal internal
  pair energies are subject to the same range rule, not hidden in D.

In the original registers the rule refers only to occupied bonds and bounded
interval support. Word positions below label an already recorded path for
the proof; they are not event-label affinities inserted in the primitive law.
The rule nevertheless **supplies universal relational locality of generators**.
Deriving that selection from something deeper remains open.

This excludes the charged bare and dressed stars, the root-blind all-pairs
terms, long Wilson transfers at growing distance, and UC1's internal X or Z
products at distance>R. A finite-range multi-step transfer is permitted and
is included in the bound. Finite register count alone is insufficient.

## 2. Complete-H comparison — including the movable end and receiver

Fix a recorded word w of length m and all its actual admissible conditional
storage/receiver state. Initially all carriers lie at positions <=b. Choose
an integer ell with

    b+max(R,2) < ell < m.

Let X be this initial subspace; for a source S require S X subset X, or
enlarge b and X to include its local support. No field reset is imposed.
Let K be the space retaining the first ell recorded labels in order and
all carriers within that prefix. Keep every tail geometry, allowed receiver
state and storage mask. Flux on the tail is zero by UC1's exact tree Gauss
solution, because all compensating charges and all carriers are in the prefix.

Define a reference **on the same physical Hilbert space**. On K it has:

* the finite-prefix matter Hamiltonian H_in, with every permitted term fully
  inside the prefix, reflecting at its end;
* the full assisted tail append/pop law, reflecting at word length ell,
  with all receiver and storage configurations retained;
* the original diagonal storage/receiver Hamiltonian and on-site clocks;
* no bare geometry term and no matter term crossing the prefix boundary.

On K-perp use K-perp H K-perp; set reference cross blocks to zero. This
defines a self-adjoint reference H_ref. It is never the physical law. It
does not assert that a projected or cooled model is the experiment.

On K the exact tree-kernel coordinates factor the prefix matter from the
tail and receiver. Prefix storage bits, if superposed, are a direct-sum
spectator label: they are not changed by this reference. The prefix part
of the storage diagonal is a spectator constant in each such block. Tail
flux is zero; no tail term reads a prefix matter internal state or position.
Consequently local matter response under H_ref is precisely response under
H_in. Correlated conditional receiver states are allowed; the operator bound
below does not assume a product initial prefix/receiver state.

Take D to be the retained storage/receiver Hamiltonian and on-site clocks.
Its norm may grow with N. D does not change word, carrier positions or
prefix support. Its interaction picture preserves the relevant support:
receiver evolution is independent of matter; storage phase changes in an
append/pop read only its changed event; on-site clocks rotate only that site.
All full changing terms have norm <3 (g_max+4epsilon+w<=2.6); all reference
changing terms on K have norm <2 (g_max+w<=1.6). In particular

    ||(H-H_ref) K|| <= 5.                                     (R1)

The common D cancels in R1; no extensive ||H|| estimate is used.

For completeness, the assisted append map has a unique parent for each
destination. At fixed Q=sum q_i its column squared norm is at most
(g/N)^2 (N-Q)Q<=g^2/4: an appendable zero-storage event is among the N-Q
zero bits and n_contact+1<=Q. Adding its adjoint therefore has norm <=g.
Zero-flux gates and the reflecting prefix compression do not increase it.
Each species' directed bare append map has norm <=epsilon by the same
unique-parent count; its adjoint doubles that, for at most 4epsilon across
both species. These bounds include all allowed bosonic occupations.

## 3. First possible omitted action — a support count, not a simulation

Put c=max(R,2) and

    k = min(m-ell, ceil((ell-c-b)/R)) >=1.                     (R2)

Every reference interaction-picture word of order n<k is annihilated by
(H-H_ref) when applied to X. To see this, expand each word in its actual
primitive matrix transitions, including all tail labels and receiver modes.

After n operations the head cannot have retreated by more than n, and no
carrier can have advanced by more than nR. An omitted assisted pop first
acts at length ell, requiring at least m-ell previous pops. An omitted bare
append/pop needs a carrier on the head or its parent; with head length at
least ell this requires a carrier at position >=ell-2. An omitted matter
term crossing the cut can act only with a carrier in its interval, hence at
position >=ell-R; it annihilates the interval vacuum otherwise. These are
all omitted actions on K. Before k operations none is available.

This includes arbitrary extra append/pop detours and repeated returns.
On-site clocks and storage/receiver propagation cost zero in this count and
cannot open an omitted action. The argument does not exclude their states
from physical evolution or assume independent emission histories.

Writing tau_q(x)=sum_(n>=q) x^n/n!, the **complete** Dyson series gives

    ||(H-H_ref) U_ref(s) X|| <= 5 tau_k(2|s|).

Exact Duhamel with the full physical U on the other side then gives

    ||(U(t)-U_ref(t))X|| <= Delta_k(t),
    Delta_k(t)=min(2, (5/2) tau_(k+1)(2|t|)).                  (R3)

No truncation of actual paths is made: all physical receiver states, storage
errors, geometry changes, visits outside K and coherent returns are contained
in the physical side of Duhamel. The total-norm Dyson sum pays for the number
of paths. If the initial storage is consistent, the reference stays consistent;
physical storage inconsistency or prefix loss each has probability at most
min(1,Delta_k^2). This is a finite-time statement, not permanent locking.

## 4. Information locality, not only arrival locality

Let S and O be number-preserving physical matter observables whose fixed-
prefix representatives have disjoint interval supports separated by recorded
distance d>0. Set a=ceil(d/R). Under H_in, the interaction-picture non-on-site
part has norm <=w. Each nonzero nested primitive commutator extends support
by at most R. Thus the first a-1 orders of [S,U_in^* O U_in] vanish. On-site
rotations preserve support, and summing the entire series yields

    ||[S,U_in(t)^* O U_in(t)]||
      <= 2||S||||O|| min(1,tau_a(2w|t|)).                     (R4)

The local product-site extension ensures that disjoint representatives
commute before restriction to the fixed-number physical sector. Their
number-preserving restrictions obey the same commutator identities and norm
inequalities. This step would not be justified for an arbitrary nonlocal
physical eigenbasis encoding.

For any normalized state in X, the actual pulse response is
chi=i<[S,U^* O U]>. Apply R3 to both psi and S psi; the two bilinear terms
give |chi-chi_ref|<=4||S||||O|| Delta_k. Hence

    |chi(t)| <= min(2||S||||O||,
      2||S||||O|| tau_a(2w|t|)
      +4||S||||O|| Delta_k(t)).                              (R5)

This applies to internal observables, number observables and local gauge-
dressed motion operators within the declared physical site algebra. It is
uniform over the actual conditional receiver state and storage disorder.
It is **conditional on the record and carrier setup**, not a free distant
charged-receiver preparation or an unconditional all-state operator-norm
bound on the entire changing-geometry Hilbert space.

Since tau_q(x)<=exp(e x-q), for

    |t|<=min(a,k+1)/(8e),

R5 is bounded by

    2||S||||O|| exp(-3a/4)
      +10||S||||O|| exp(-3(k+1)/4).                            (R6)

Both d and the head/carrier buffers must grow. For example m>=ceil(N/2),
b=floor(N/16), ell=floor(N/4), fixed R and sufficiently large N give
a,k proportional to N for receivers a comparable distance from the source.
The guaranteed observation window is O(N). It is not the O(N^2) formation
window, a real-energy Schur theorem, or a universal relativistic cone.

## 5. Nonempty class and formation price

As an existence witness, keep UC1's charged hoppings j_A=1/64, j_B=1/96,
B's distinct third-level clock, and replace its internal counterexample by

    H_bond=(1/4) sum_(i!=j) B_ij Z_Ai X_Bj.                    (R7)

This is RI(1,1/2): ||H_bond||<=1/4 on the two-carrier blocks, and
2j_A+2j_B+1/4=29/96<1/2 on every word and reflecting prefix. It couples
genuinely different species and has nontrivial internal dynamics. All full
geometry, receiver and storage terms remain in the physical law. R7 is a
declared restricted-class witness, not a newly adopted WM2 successor.

In either single-species experiment H_bond vanishes. The archived UC1
adjacent gauge-pair phase protocol is therefore unchanged: its source has
variance 1/4 and at t=10^-4 both nearby responses exceed 9*10^-7 by the
same complete double-Dyson bound. Finite diagnostics replay that statement.
The restriction permits viable transport; it does not force every allowed
Hamiltonian to have nonzero coefficients or prevent deliberate cancellation.

In joint matter vacuum all matter generators vanish, all zero-flux gates
are one, and all matter-dependent bare terms vanish. This is exactly the
original WM2 formation law, including histories and reverse transitions.
No new numerical formation evidence is required for this isometry. S9's
formation conclusion is inherited only at its existing author-derived,
independent-review-pending status and with its original large resources.
Recording, injection, and any separated charged-receiver setup stay supplied.

## 6. Stability and boundaries

R3–R6 cover the entire RI class with the same R and norm caps. Small
perturbations preserving this generator class and the original geometry
and receiver rules remain covered. They do not cover all small-norm
constraint-preserving perturbations: a fixed tiny distant product has the
nondecaying double commutator O9, and on separated prepared carriers gives
an immediate nonzero response derivative. The exact O4 curve requires its
factorization hypothesis; compulsory extra terms must instead be controlled
as in O5. Already the derivative rules out a distance-uniform exponential
bound with a uniform differentiable zero-time prefactor. The locality class
is not open in that unrestricted perturbation sense. That is the assumption
price, not a numerical defect.

No theorem here derives RI, a shared speed, dimensional selection, a
finite-density matter theory, a restored spectral gap, gravity, or an
empirical prediction. The result is that a **supplied universal rule can
control the full finite-time experiment**. A selection mechanism for that
rule is a separate unsolved construction problem.
