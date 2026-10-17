# S9 transport and composition: explicit comparison space

Same-author adversarial reconstruction, not external mathematical review.
No counterexample to the scoped S9 arrival argument was found. The
construction below removes an ambiguity about the common Hilbert space;
it does not remove transitions from the physical Hamiltonian.

## Define the two operators on one space

Fix an accepted recorded word w of m>=ceil(N/2) labels, and ell=floor(N/4).
Let G be the storage-consistent geometries whose first ell labels agree
with this record. All admissible receiver fields are retained. Let C be
the subspace of the full physical one-probe space with geometry in G and
probe on one of those ell labels. Form

    H_common = H_physical direct-sum H_dummy,
    H_dummy = span{geometry in G, field, probe index x>=ell}.

C plus H_dummy is naturally G-with-field tensor the half-line probe.
The physical operator is the original H on H_physical and only the
configuration-preserving D on H_dummy. The reference operator is the
prefix-restricted formation Hamiltonian plus independent half-line
hopping on C plus H_dummy, and only D on its orthogonal complement.

D is the actual storage/disorder energy plus receiver Hamiltonian on
physical states; on dummy states it is the same consistent-geometry
expression. It is common to both operators, is self-adjoint, and
preserves all these components. Each extension has zero changing
interaction on its unused component. Thus the physical changing norm
stays <=nu=2.1 and the reference norm stays <=g+2J<=1.6. A union of
both interactions is not used as a physical Hamiltonian.

The recorded field need not be vacuum. A pure conditional field vector
or a purification of a mixed conditional field is an admissible initial
state in C; all estimates are operator/state-uniform there.

## Check the first differing word

All words below are products of changing interactions in the exact
interaction picture of D, with arbitrary ordered insertion times.
D may mix receiver modes but cannot move the probe or change a word.

Starting with probe support at recorded indices {0,1}:

* A difference at the prefix's probe boundary needs at least ell-1
  probe hops, including the hop from ell-1 to ell.
* The formation cut first differs after m-ell+1 end pops.
* The first bare transition needs at least m-2 changes, including that
  transition. A pop can first act when the head or parent reaches the
  probe; simultaneous probe advance and head retreat cannot reduce
  their combined count. Detours cannot shorten it.

Since m>=2ell, all products of length n<q=ell-2 agree. This count is
deliberately weaker than the first-boundary count and has slack. It
applies to full amplitudes, even when histories coalesce or leave and
return later. The new diagnostic directly multiplies V_I(t_n)...V_I(t_1)
with complex field states, nonuniform detunings and exact two-mode D
evolution. Its agreement is not inferred from support or populations.
Deleting one shared probe edge is detected at order one. These are
finite diagnostic receivers, not physical-L formation calculations.

The two complete Dyson tails give

    ||U_phys(s)psi-U_ref(s)psi||<=2 tail_q(nu|s|)
                              <=2 exp(-3q/4),
    |s|<=q/(4e nu).                                               (T1)

No small norm of D, spectral gap, or exclusion of high receiver energy
appears. State-distance conversion costs at most twice (T1) for any
common probability. The reference's probe tail beyond ell costs at
most exp(-3q/2). Both costs must be retained.

## Recheck the reference's arrival probability

For the half-line shift T, H_p=-J(T+T^*), X position, the exact bounded
velocity is V=iJ(T-T^*) and i[H_p,V]=2J^2 |0><0|. For the real
adjacent pair, <V>(0)=0, p_0(0)=1/2. Since |p_0'|<=4J, on
[0,1/(16J)] one has p_0>=1/4. Integrating twice implies

    <V>(s)>=J/32 for s>=1/(16J),
    <X>(s)>=Js/64 for s>=1/(8J).

The identity X(s)=X+integral V(u)du on finite-support data yields
sqrt(<X^2>(s))<=1/sqrt(2)+2Js<=3Js for s>=1/J. Cutoffs and the
bounded velocity justify taking the position-domain limit. Applying
Paley-Zygmund at half the mean, and then weakening its event to
X>=ceil(Js/128), gives 1/147456. The ceiling is legitimate for integer X.

Subtract the two terms from (T1) and the reference exterior tail:

    p_phys[d(s)<=recorded index<ell]
       >=1/147456-4exp(-3q/4)-exp(-3q/2).                         (T2)

For q>=32, this exceeds 1/294912. It is arrival in a receiver set for
an injected particle, not a distant single-site phase susceptibility.

## Audit composition, not only individual theorems

| Stage | Required state or operation | Check |
|---|---|---|
| Formation | Original root weights; probe vacuum; physical receiver vacuum | The exact consistent sector is invariant; S9 is root-uniform |
| Time | Frozen uniform sample in [4N^2,8N^2] | Pointwise bound implies the same lower bound after sampling |
| Acceptance | Every m>=ceil(N/2), with its actual probability | No selected ordering or receiver packet |
| Record | Entire word measured | Pre-injection storage is its fixed mask, n_R=m-1 |
| Conditional field | Actual post-record field, possibly mixed | T1 and S4 are uniform over this whole field sector |
| Source | (|w_0>+|w_1>)/sqrt(2), then a local phase pulse when measuring chi | Variance of n_(w_0) is exactly 1/4 |
| Arrival | The real pair without a phase pulse | Same supplied injection; a separate observable from chi |
| Rejection | Zero signal, original probability retained | Multiply conditional lower bounds by p_acc, never renormalize it away |

S4's storage-error amplitude is <=(epsilon/nu)exp(-3q/4) in the common
window. Prefix-loss amplitude is <=exp(-3q/4). These are single-time
guarantees, not monitored survival probabilities. Subtracting both
squared errors from (T2) gives a correction <7exp(-24)<1/294912 for
q>=32. Thus the claimed joint accepted/consistent/prefix/arrival lower
bound 1/589824 follows if the formation theorem holds. The statements
compose for all conditional fields; no resetting is needed.

## Clarifications a reviewer should preserve

The response front's distance is from the **source support**. For a
receiver at recorded root distance D and the adjacent-pair source,
the elementary hop count is D-1, not D. Use tail_(D-1) when expressing
that upper bound in root coordinates. S9's q=ell-2 comparison and its
arrival lower bound already allow this offset; no asymptotic constant
or claim needs to be changed. This is a notation clarification, not a
new faster-than-front process.

The record-energy bound |Delta E_system|<=g is an aggregate system
mean-energy change, not apparatus work or an accepted-branch bound.
Injection, phase control, record storage and apparatus reset remain
supplied. The fixed-energy gap is not restored. The conditional window
grows as N, formation as N^2; their ratio shrinks. This review supplies
neither indefinite locking nor an equilibrium spatial phase.
