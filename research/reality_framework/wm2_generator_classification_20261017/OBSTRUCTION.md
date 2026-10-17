# IG1: the gap between observables and elementary interactions

**Tier: author-derived; independent mathematical review pending.**
The result is an obstruction to a specified Hamiltonian class, not a no-go
for every constraint, physical encoding, or microscopic theory.

## 1. Three different objects

For the finite-N model let Q_alpha denote the supplied commuting constraints,
P their joint zero projector, and H_amb the ambient register space. Define

    C = {A in B(H_amb): [A,Q_alpha]=0 for every alpha},
    A_P = {A in B(H_amb): [A,P]=0},
    A_phys = {P A P restricted to P H_amb: A in A_P}.             (O1)

C is the full constraint commutant; C is contained in A_P. In general they
are different: preserving the zero sector need not preserve every other
constraint eigenspace. As an abstract finite-dimensional algebra A_phys is
all B(P H_amb), since any operator there has a block-diagonal extension.
That abstract statement supplies neither a local register implementation
nor a cheap preparation. This packet uses explicit, bounded-support
representatives in C, not arbitrary P h P gates.

Indeed unrestricted compression of C gives the same A_phys: extending an
arbitrary physical operator by zero on P-perp commutes with every Q_alpha,
because each Q_alpha is zero on P. The ambient algebras differ, but their
unrestricted physical representations do not select interactions. Local
representatives and implementation cost, not this abstract compression,
are indispensable extra data.

A Hamiltonian specification is additional data:

    G_H subset of (A_phys)_sa,
    H = sum_a c_a g_a, with declared support and norm conditions. (O2)

G_H need not be closed under multiplication, or contain every self-adjoint
element of A_phys. A physical observable can be measurable without its
product with any distant observable being an elementary interaction.
Self-adjoint observables themselves are a real vector space; the underlying
complex observable algebra is the object closed under products.

## 2. Product-closure obstruction — precise hypotheses

Assume a sequence of accepted records has sites i_N,j_N at recorded distance
d_N tending to infinity. Assume:

* There are bounded-support physical Pauli representatives at each site,
  commuting with all Q_alpha, with independent commuting M_2 algebras.
  On a reducing physical subspace they act as M_2 tensor M_2 tensor 1.
* The subspace admits the two-qubit product preparation used below, a local
  source pulse, and a local readout. Preparation may require a paid earlier
  charged string. It is identical for every subsequent source choice.
* The admitted Hamiltonian class contains the product kappa Z_i X_j with
  kappa>0 independent of N, together with a spectator Hamiltonian independent
  of these internal qubits on this subspace. For an added mandatory term use
  the perturbation statement below instead of assuming it vanishes.

Then the class contains a Hamiltonian with fixed-time distant signaling.
Indeed products commute with every constraint by the Leibniz identity

    [Q_alpha,A_i B_j]=[Q_alpha,A_i]B_j+A_i[Q_alpha,B_j]=0.        (O3)

On the two qubits take

    H_0 = kappa Z_A X_B + H_rest,
    psi = |+y>_A |0>_B tensor phi_rest,
    S = X_A/2, O = Y_B, psi_theta=exp(-i theta S)psi.

The exact receiver expectation and response are

    <O(t)>_theta = -sin(theta) sin(2 kappa t),
    chi_0(t) = -sin(2 kappa t), Var_psi(S)=1/4.                 (O4)

At t_*=pi/(4 kappa), |chi_0|=1 regardless of d_N. No charge moves, no
flux changes, and no virtual geometry activation is needed. Noncommuting
source/readout algebras and the admissible preparation are essential
hypotheses; product closure alone is not a signaling proof.

If H=H_0+K with ||K||<=h on the invariant sector, exact Duhamel gives
||U-U_0||<=h|t|. Consequently

    |chi-chi_0| <= 4||S||||O|| h|t| = 2h|t|.                   (O5)

This includes the complete K dynamics, not a first-order expansion. Whenever
2h t_*<1, a size-independent response survives. Arbitrary compulsory
backgrounds need not satisfy this small-norm hypothesis; no claim is made
that every background or every allowed Hamiltonian signals nonlocally.

## 3. Covariance and a uniform total norm do not remove this example

For UC1's one-A, one-B sectors use the occupation-tagged covariant sum

    H_* = kappa sum_(i!=j) Z_Ai X_Bj.                          (O6)

Different ordered position pairs are orthogonal blocks. Thus ||H_*||=kappa,
not an extensive sum of norms. It is root-blind, label-covariant, physical,
number conserving, supported on two event registers per summand, and inert
in matter vacuum. Fixing the two carrier positions realizes O4.

UC1 retains the full geometry/receiver dynamics and both nonzero hoppings,
with h<=2j_A+2j_B=5/96, kappa=1/4. At t_*=pi its archived derivation gives
|chi|>7/12 and prefix survival>15/16. The current obstruction includes that
example; its conditional setup and paid O(d_N)-link receiver preparation
remain part of the scope. There is no new original-trial probability here.
Neither O3 nor O4 uses S9's formation argument.

## 4. A diagonal term can signal too

Replacing O6 by kappa sum_(i!=j) Z_Ai Z_Bj preserves every UC1 constraint
and changes no configuration-basis occupation, position, flux, or geometry.
For A initially |+y>, B initially |+x>, the same local X_A/2 pulse and
Y_B readout instead give chi(t)=sin(2 kappa t).

Therefore a proof cannot place *all* configuration-diagonal terms into an
unrestricted harmless D. A term whose coefficient reads a distant internal
register has that register in its support even if it changes only the
receiver's basis value. 'Registers changed' is not interaction support.
Only on-site or otherwise already controlled diagonal contributions may
be taken out of the information-propagation budget without further proof.

## 5. A necessary infinitesimal test, independent of decomposition

For commuting physical local algebras at separated i,j define

    Gamma_H(i,j) = sup_(||A_i||,||B_j||<=1)
                  ||[A_i,[H,B_j]]||.                         (O7)

Suppose the **state-uniform operator** locality claim is differentiable
at zero and has the form

    ||[A_i, exp(itH)B_j exp(-itH)]||
      <= C ||A_i||||B_j|| [exp(v|t|)-1] exp(-mu d(i,j)).        (O8)

Since [A_i,B_j]=0, divide by |t| and take t to zero:

    Gamma_H(i,j) <= C v exp(-mu d(i,j)).                       (O9)

For kappa Z_A X_B choose normalized A=X_A, B=Y_B. The double commutator
has norm 4kappa. Thus a uniform O8 forces kappa<=Cv exp(-mu d)/4 for
this isolated product. Constraint commutation imposes no such inequality.
Polynomially small couplings can also fail a uniform exponential claim.
The test concerns the actual H; cancellations and alternative decompositions
must be assessed there, not by treating every displayed summand as unique.

O9 is necessary under O8, not a standalone sufficient theorem for a changing
geometry, and not a statement about a single specially prepared response.
Wilming and Werner prove a stronger converse for finite-body potentials
under their fixed tensor-lattice hypotheses; we do not export that theorem
to the full WM2 configuration space without those hypotheses.

## 6. What is and is not classified

| Specification | Outcome |
|---|---|
| All bounded, constraint-preserving two-event products, independent physical internal qubits, nondecaying strength | Contains a signaling witness |
| The same class with UC1 covariance and a uniform few-carrier norm | Still contains O6 |
| A chosen Hamiltonian omitting the witness | Not ruled out by this theorem |
| Only an abelian physical algebra, no independent qubits, or no admissible local experiment | O4's hypotheses fail; no general conclusion here |
| Delocalized/error-corrected physical encoding | Requires a separate definition of operational local observables and setup cost |
| Extra finite energy penalties for the existing constraints | Do not penalize O6, which lies inside their exact kernel |

The conclusion is insufficient interaction selection, not failure of quantum
mechanics, gauge theory, geometry formation, or every possible WM2 successor.
