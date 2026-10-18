# A full-law finite-order test of moving exchange

Tier: author-derived finite mathematical diagnostic, not a macroscopic
operating result. The initial contact boson and short record are supplied
as frozen in CHARTER.md. They do NOT stand in for S9's actual conditional
receiver. No claim transfers from this diagnostic to its preparation class.

Use the physical N=4 B1 law with all1024 receiver modes, all word/storage
states and reverse transitions. Gauss uniquely fixes the capacity-one path
field; eliminating that deterministic redundant register deletes no physical
state. Work in the invariant two-resource, inactive-fuel-vacuum sector.
The unnormalized Fock basis has metric w(n)=prod_j n_j!; creation coefficient
one, annihilation n_j and receiver hop -n_j. Weighted Hermiticity is
h_(t,s)w(t)=h_(s,t)w(s). Color permutation is carried explicitly as a bit.

Alpha/beta are the accepted environment states defined in the charter;
their Fock norms are1 and2. Write A=(chi_++chi_-)/2 and B=(chi_+-chi_-)/2
for identity/exchange vectors. The exact primitive graph predicts:

| Component | Minimal primitive length | (H^length) coefficient, unnormalized Fock basis |
|---|---:|---:|
| A_alpha | 2 | 2J^2 |
| A_beta | 3 | -6(g/N)J^2 |
| B_beta | 7 | -3(g/N)J^2 jp^4 |
| B_alpha | 8 | (13/16)J^2 jp^4 at N4,g1 |

The six nongeometry operations for B_beta form the chain root-convert,
plus0->1,plus1->2,right-convert,plus1->0,plus2->1. The append of2 has
three allowed positions before plus1->2. Minus particles need NOT move:
Gauss fields carry the temporary separation and return to zero. This is
why tracking identifiable composites alone would miss the shorter route.

B_alpha includes both free children2 and3. Assisted append/pop contributes
12(g/N)^2 J^2 jp^4 (including the two-boson annihilation factor2).
Bare append/pop contributes4(1/8)^2 J^2 jp^4: the occupied-head factors
are1 on append and2 on pop. Omitting those paths changes the exact answer.
Every shortest path has negative hopping/creation coefficients with the
same parity of length, so its leading amplitude cannot cancel. Diagonal
and receiver terms are retained by the law; a self-loop or extra receiver
hop cannot shorten these distances. The radius8 enumeration coherently
sums all minimal words, not just the displayed witness.

## What this says about recoverability

The convergent Taylor expansions give

    A_alpha=-J^2 t^2+O(t^3),
    A_beta=-i sqrt(2)(g/N)J^2 t^3+O(t^4),
    B_beta=-i sqrt(2)3(g/N)J^2 jp^4 t^7/7!+O(t^8),
    B_alpha=O(t^8).

Consequently the normalized-component minor
W(t)=A_alpha B_beta-A_beta B_alpha has a nonzero coefficient at order9:

    W(t)=i sqrt(2) C0 t^9+O(t^10),
    C0=J^2 |(H^7)_(beta,odd)|/7! >0.

This distinguishes a common unitary exchange from an environment-correlated
successful operation. A and B are not proportional on a sufficiently short
interval. The environment includes orthogonal photon-number/tail records,
but these are traced out by the frozen diagnostic decoder. The successful
Gram determinant obeys uv-|c|^2=4 det Gram(A,B)>=4|W|^2>0.
By the recovery theorem no deterministic color-only decoder can perfectly
restore all four dimensions in that interval. The triplet subspace remains
exact. This is NOT a claim of fundamental information destruction.

## A conservative quantitative interval, not a useful operating time

For the COMPLETE N4 Hamiltonian, ||H_+||,||H_-||<=512 suffices. Bound its
terms by48 (16(Q-m)),192 (storage mismatch),6(receiver hopping),3(field),
4(mass),4(fuel),1(density),4A(degree), and offdiagonal geometry/matter/
conversion by2+5/48+2J. Their sum is below512. Nominal disorder is zero.
The bound includes unreachable high-energy states and every receiver mode.

Normalized matrix-amplitude Taylor coefficients are bounded by512^k/k!.
The minor's coefficients are thus bounded by2*1024^k/k!. On |t|<=1/1024,
the remainder after degree9 is at most6(1024|t|)^10/10!. Define exactly

    t0=min(1/1024,C0*10!/(12*1024^10)).

For 0<t<=t0, |W|>=C0*t^9/2 and uv-|c|^2>=C0^2*t^18>0.
The calculator records t0 as a rational, without floating-point decisions.
This intentionally conservative interval is extremely small and says
nothing about violating the 1% operating tolerance or useful delivery.

## What the finite computation does and does not cover

Radius8 is a COMPLETE finite-order calculation for this actual finite
Hamiltonian. It is not a complete time evolution or a replacement bath.
The graph ball is complete; accumulated weights are the minimal-order
coefficients at the four targets, not every nonminimal coefficient of
H^k xi through order8. The norm estimate pays for those remaining terms.
All1024 physical modes remain available in the neighbor generator; modes
outside the graph ball cannot affect these shortest coefficients. The norm
bound controls the subsequent analytic remainder rather than discarding it.
No continuous-window populations, large-size probabilities or actual S9
conditional receiver is calculated. Small-time nonzero noise is compatible
with excellent approximate operation later. The original task is unresolved.
