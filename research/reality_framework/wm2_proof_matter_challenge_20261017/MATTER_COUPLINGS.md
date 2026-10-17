# A distinct matter sector: locality is not forced by the declared symmetries

**Negative architecture result, in the explicit class C of CHARTER.md.**
It is not a counterexample to S9 for its original neutral probe, nor an
impossibility theorem for every extension of WM2. It answers the bounded
question: the current symmetries and charges do not force a new matter
sector to share the generated path metric.

## 1. State the extension without altering the benchmark

Call the original probe a. Add a distinct neutral species b with vacuum
and a one-particle event space. In the b experiment a is vacuum, so the
original a-conditioned bare append/pop terms vanish exactly. Let H_G be
the unchanged storage-consistent formation Hamiltonian, with every
assisted append/pop and every physical receiver state retained.

Every new term below annihilates b vacuum. Consequently the pre-injection
preparation and formation law are exactly S9, with the same root weights,
detunings and resource sequence. The new Hamiltonian is explicitly an
extension called **WM2-MC**, not a silent edit to WM2 or CP-12. Norms here
refer to the declared vacuum/one-particle sector, not finite-density
many-body matter or a universal stress-energy theory.

## 2. The coupling class has no missing symmetry that removes shortcuts

In one-particle notation, two-event hopping is
sum_(i<j) A_ij (b_i^* b_j+b_j^* b_i), with real geometry operators A_ij
that transform covariantly with event labels and preserve the storage
charge. On-event terms are also allowed. Finite register support of a
summand does not mean finite range in a metric not yet supplied.

| Permitted example | Respects recorded path? | Reason |
|---|---|---|
| Bounded on-event potential | No transport by itself | Preserves positions, can change scattering |
| Occupied-bond hopping with bounded coefficients | Locally, with WM2 word barriers | The bond projector is a supplied coupling rule |
| Fixed bounded-length occupied-bond words | A finite propagation range on a retained prefix | Range restriction is an additional assumption |
| Normalized graph-blind all-pairs hopping | Violates uniform exponential tails | Symmetry and charge permit it |
| Root-controlled normalized star | Constant far-set arrival in fixed time | The root is an existing distinguished register |
| That star plus nonzero occupied-bond hopping | Still violates the metric | Nonzero nearby transport does not cure its shortcut |

The graph-following rows describe a sufficient class, not a class selected
by the present symmetries. In particular, setting every off-bond coefficient
to zero is not derived by listing it as a Hamiltonian convention.

## 3. Root-blind obstruction, even with a uniform norm

The fully permutation-invariant, geometry-blind hopping

    H_all=-(kappa/(N-1)) sum_(i!=j) |i><j|, kappa>0,              (M1)

has norm kappa, conserves b number, commutes with H_G and every
geometry/storage observable, and acts on two event registers per summand.
It introduces no fixed pair preference or event metric. Its distinct-site
amplitude is exactly

    <j|exp(-it H_all)|i>
       =exp(-i kappa t/(N-1)) [exp(i kappa Nt/(N-1))-1]/N.

At t=1/kappa its magnitude is at least .9/N for every N>=2. For
recorded distances proportional to N this is incompatible with any
bound C exp(vt-mu d) having C,v,mu>0 independent of N. Polynomial
suppression is not exponential spatial suppression. This counterexample
requires no reading of the distinguished root.

M1's far-set population from one localized event shrinks with N. We do
not promote this weaker example to a constant-probability bypass. The
next example supplies that stronger operational consequence.

## 4. Exact star dynamics with the same local source recipe

Let r be the conserved root of the recorded word, M=N-1, and
u_r=M^(-1/2) sum_(j!=r)|j>. Define

    H_star=-kappa(|r><u_r|+|u_r><r|), kappa=1/4.                 (M2)

Before recording the same formula is block-controlled by the existing
root register. Under a relabeling, r and every b label relabel together.
In register notation each term reads a root indicator at one event and
transfers b between two events. The coefficient is kappa/sqrt(N-1),
not an unbounded all-pairs strength. Norm H_star=kappa. It is real,
number/charge conserving and root-covariant. It commutes with H_G.
No new field, drive, reservoir or geometry transition has been added.

Inject exactly the real adjacent root pair used for a, but in b:
(|r>+|w_1>)/sqrt(2). A local phase exp(-i theta n_r) has number
variance 1/4. Put s*=pi/(2kappa), independent of N. For every label
j other than r,w_1, exact evolution gives

    p_j(s*;0)=1/(2M)+1/(2M^2),
    chi_j(s*)=-1/M^(3/2).                                        (M3)

The phase derivative is polynomially small but nonzero and independent
of recorded distance. The injection probability gives a stronger result.
Let ell=floor(N/4), d=floor(N/8), and predeclare the far receiver set

    F_w={w_d,...,w_(ell-1)}.                                     (M4)

This set lies inside the same original N/4 prefix, at root distances at
least d; for N>=16 it excludes both source labels. Its size is
ell-d>=N/16, so by (M3)

    p_star(F_w;s*)>1/32.                                        (M5)

This is independent of the conditional receiver, disorder and word
ordering, and is not favorable source or receiver postselection.

## 5. Nonzero nearby transport does not remove the counterexample

To exclude a trivial objection to M2, retain a second nonzero term

    V_loc=-j_b sum_(occupied consecutive bonds i,j)
                      (|i><j|+|j><i|), j_b=kappa/1024.
    H_MC=H_G+H_star+V_loc.                                      (M6)

Every nearest-neighbor bond matrix element has magnitude at least j_b,
independent of N. Thus the b sector includes useful local hopping, with
a different speed and an additional allowed mechanism. V_loc depends
on the changing geometry and need not commute with H_G. Its full norm
is <=2j_b, however, on the entire geometry/receiver/one-b space.
Exact Duhamel with the full operators, including all geometry histories,
gives a state difference <=2j_b s*. The probability change for F_w is
at most 4j_b s*=2pi/1024<1/128. Consequently

    p_MC(F_w;s*)>1/32-1/128=3/128.                              (M7)

No perturbation series, field reset, or deleted reverse path occurs here.
The small local coefficient is a fixed witness, not fitted to N or the
record. An open interval 0<j_b<kappa/(256pi) would suffice for the
same strict probability margin; universality fails beyond one point.

## 6. The geometry persists during this bypass

In the b experiment storage consistency is exactly invariant: b has no
extra bare geometry coupling and a remains vacuum. Put all receiver and
configuration-preserving energies in D. The changing norm is

    W<=g_max+kappa+2j_b<7/5.

Losing the first ell recorded labels requires at least
q_p=m-ell+1>=floor(N/4) end pops. Other interactions cannot reduce
this count. For N>=512, q_p>=128, eWs*<34<q_p/2, so the complete
Dyson tail bounds prefix-loss amplitude by exp(-q_p/2), and its
probability by exp(-q_p)<1/256. Combining with (M7),

    p_MC(prefix retained and b in F_w at s*)>1/64.                (M8)

This is a conditional-record theorem for every admissible field, without
using S9 formation. If S9 formation survives external review, original
trial acceptance >=1/2 gives

    Pr(accepted, consistent, prefix retained, b in F_w at s*)
       >1/128, N>=2^28.                                       (M9)

The recorded distance grows as N while s* stays fixed. For the original
a species, reaching the same receiver set from its adjacent pair needs
at least d-1 changing moves: otherwise neither the interior bonds nor
the source can have encountered the movable end. Its probability is at
most tail_(d-1)(2.1s*)^2, which vanishes exponentially for fixed s*.
Thus the two distinct sectors disagree operationally on neighborhoods,
despite persistent recorded geometry and nonzero local hopping for b.

## 7. Extra-sector resource ledger

WM2-MC supplies a new vacuum/one-particle register, the constants kappa
and j_b, the root-controlled support rule and its 1/sqrt(N-1)
normalization. These are explicit witness choices, not derived constants
or a proposed law of nature. The root-blind witness instead supplies
the 1/(N-1) normalization. In a fixed root block M2 has N-1 active
spokes; its covariant expression ranges over all root labels. Bounded
single-particle norm does not imply cheap implementation or a sparse
set of permitted event pairs.

The real pair has b interaction energy -kappa/sqrt(N-1)-j_b in M6.
Injection is external; an arbitrary absolute particle rest energy and
apparatus work remain supplied, as for a. The phase pulse changes this
pair energy by [kappa/sqrt(N-1)+j_b](1-cos(theta)). No extensive hidden
fuel, receiver reset, or time-dependent geometry drive is used. The
counterexample demonstrates what the declared coupling class permits;
it is not an efficient implementation proposal.

## 8. What has and has not been decided

M8 is a constructive negative to universality in class C; M9 is its
formation-conditioned original-trial corollary. This does not say that
every matter sector is nonlocal, that S9 is false, or that a future
selection mechanism is impossible. It does show that label covariance,
finite primitive support, neutral number conservation and the current
storage charge are insufficient. Root-blind M1 also prevents solving
the problem merely by forbidding use of the root flag.

To exclude these examples one must add and justify a restriction: for
example, an occupied-bond selection rule shared by all matter, a genuine
gauge/constraint principle enforcing it, or a dynamical argument bounding
all permitted off-bond terms. The latter must include direct terms like
M1/M2, not only virtual geometry shortcuts. None is derived by S9.
Whether such a principle can be obtained is a new architectural problem;
this bounded campaign stops at the counterexample.
