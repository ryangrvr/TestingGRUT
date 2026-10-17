# WM2-UC1: explicit constraints and invariant microscopic law

All results here are **author-derived, independent review pending**.
The new constraints, charges and registers are **supplied**. This is a
separately named successor matter architecture, not unchanged WM2.

## Registers and physical space

Retain WM2's full ordered-word geometry w, storage mask q, charge-allowed
bosonic receiver n_R=|q|-1, and all storage-inconsistent configurations.
Geometry coefficients receive the declared zero-flux gates below, with
every exact reverse term retained. B_ij(w) is one exactly when labels
i,j are consecutive in w. No fixed event-pair affinity or coordinates
are supplied. The path dimensional class remains grammar-supplied.

Add one independent register E_ij=-E_ji per unordered event pair with
eigenvalues -3,...,3. An index ordering only chooses an orientation;
relabeling an oppositely oriented pair also reverses the sign of E.
Define U_(i->j) to reduce E_ij by one (and increase E_ji). It is a
partial shift, zero at the finite boundary, not a unitary cyclic shift.
It obeys [E_ij,U_(i->j)]=-U_(i->j) exactly, including at the boundary.

Species A has charge z_A=1 and two internal levels; B has z_B=2 and
three internal levels. There is a vacuum and at most one carrier of
each species, including the sector with one of each. Internal levels
have equal gauge charge within a species. The charges are model U(1)
charges, not an identification with electromagnetism or known particles.

Let N_s=sum_i n_si. Two immobile compensator occupations f_A,f_B at
the conserved root r cancel their charge. Include N_A-f_A=N_B-f_B=0
as constraints. Injection creates a carrier at r together with its
compensator. It is an external charge-neutral intervention. For
R_i=1_(i=r), impose the commuting self-adjoint constraints

    G_i=sum_(j!=i) E_ij - n_Ai - 2n_Bi + R_i(f_A+2f_B)=0,
    C_ij=(1-B_ij)E_ij^2=0.                                      (C1)

The physical space P is their joint kernel with the number/compensator
constraints. C1 is exact kinematics, not a finite energy penalty or an
isolated low-energy band. G_i may read all N-1 potential incident pair
registers; the supplied constraint is not constant-size spatial locality.

### The kernel is explicit, not an arbitrary projected eigenbasis

All carriers must lie on the word. At an unused event every incident
electric field is zero by C, so G cannot accommodate positive charge
there. On a word, cut the edge after position k. Gauss gives the unique
oriented flux as minus the total charge to the right of the cut, in
the root-to-head orientation. It is 0,-1,-2 or -3. There are no cycles
and hence no independent circulation solutions. This proves that the
seven electric levels contain the **entire** physical kernel for the
declared particle sectors. Fields are fixed by word and positions;
internal states, storage and receiver remain free.

`operators.py` uses precisely this isometry. It retains all label orders,
storage masks and charge-allowed occupations in its finite two-mode
diagnostic instances. It does not remove physical electric states by
approximation, or replace the physical receiver in a formation theorem.

## An explicit invariant Hamiltonian

Use the same storage/receiver diagonal as WM2. Add B's on-event energy
sum_i n_Bi |2><2| with unit coefficient. Electric energy is chosen zero;
no massless gauge mode, gap or field dynamics is inferred from that choice.

For each original assisted append and its adjoint, multiply by the
zero-flux projector Pi_(E_uv=0) on the changed head/free-event bond.
Do the same for each bare append/pop, with the original head/parent/free
gate applied to each present species and summed. Coefficients remain
g/N and epsilon/sqrt(N), epsilon=1/4. Call these H_geom. A nonzero-flux
pop is blocked; a zero-flux append/pop leaves E unchanged. These are
explicit primitives, not implementation by the global projector P h P.

For each species introduce gauge-dressed motion

    T_s(i->j)=b_sj^* U_(i->j)^(z_s) b_si,
    H_move,s=-j_s sum_(i<j) B_ij [T_s(i->j)+T_s(i->j)^*],
    j_A=1/64, j_B=1/96.                                         (C2)

Transfer is the identity on internal levels. Both coefficients are
nonzero and independent of N, with genuinely different charge and
internal spaces. No claim that their speeds agree is made.

Finally allow, as the adversarial member of the same operator class,

    H_int=kappa sum_(i!=j) Z_Ai X_Bj, kappa=1/4,                 (C3)

where Z_Ai=|i><i|_A tensor diag(1,-1) and X_Bj is |j><j|_B tensor
the Pauli X acting only on B levels 0,1, zero on level 2. These
operators are zero in the corresponding vacuum. C3 is real, root-blind,
label-covariant, number conserving and acts on two event/internal
registers per summand, with no electric-register action. Its total
norm is kappa: the different position pairs form orthogonal blocks.

The concrete candidate is

    H_UC1=D_WM2+H_B,clock+H_geom+H_move,A+H_move,B+H_int.          (C4)

Its primitives use at most one electric register; the analysis also
states a charged-motion result for any fixed r-register primitive class.

## Invariance is an operator statement

* D and the B clock commute with every constraint.
* T_s preserves G: its matter charge changes cancel the flux shift at
  both endpoints. The finite shift boundary is zero and adds no exception.
  B_ij makes C vanish on the shifted bond; other C's are unchanged.
* An append/pop changes only B_uv on a zero-E subspace, so it commutes
  with C_uv and leaves all G's unchanged. Storage charge remains conserved.
* C3 changes neither positions nor charges nor fields nor geometry, so
  it commutes with **every** G and C, not merely with their kernel.
* Compensator occupations and species numbers are conserved throughout.

Thus [H_UC1,G_i]=[H_UC1,C_ij]=[H_UC1,N_s-f_s]=0 and P is exactly
invariant. Unphysical constraint sectors are not virtual intermediate
states of this Hamiltonian: P H (1-P)=0. This is a consequence of the
declared exact law, not an argument that a large finite penalty suffices.
Every physical geometry, receiver, storage-error and return process is
still present whenever its explicitly declared primitive is nonzero.

The changing-interaction norm on P is bounded by

    W<=g_max+4epsilon+2j_A+2j_B+kappa<5/2.                       (C5)

The assisted bound uses unique parents as in WM2. Each species' directed
bare map has norm <=epsilon; its adjoint doubles that, for 4epsilon
total. The zero-flux gate cannot increase these unique-parent bounds.
Each path hopping block has norm <=2j_s. C3's bound was given above.
Receiver, diagonal energies and the B clock can be treated exactly in
the interaction picture; an extensive diagonal norm is not substituted
into C5.

## Exact formation inheritance, with a preparation price

Augment the original seed by E=0, f_A=f_B=0 and both matter vacua.
It satisfies all constraints. H_int, H_move, clocks and bare gates vanish.
Every append/pop zero-flux gate equals one; the receiver/storage law is
exactly the original WM2 formation law. The augmented zero-E matter-vacuum
component is invariant and isometric to it for every N and physical L.

Therefore S9's **author-derived** formation result is inherited without
altering its proof status, resource sequence or observation protocol.
No new formation trajectory or independent validation is claimed. Preparing
all new registers in the physical kernel and choosing that kernel as the
state space are supplied assumptions. C1 has not been dynamically selected.
The original neutral probe remains a benchmark; it is replaced, not
silently identified with one of the new charged species.
