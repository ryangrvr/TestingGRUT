# Matter moves an increasing distance before the protected prefix is lost

Builder derivation, not independently reviewed. This is an additional
observable corollary of the **same** post-recording S4 Hamiltonian and
root-pair preparation; it changes neither frozen formation sampling nor
the numerical recovery protocol. It measures injected-particle arrival,
not the phase derivative at a far receiver. All physical geometry and
storage transitions remain in the actual evolution.

## 1. Comparison without discarding physical paths

For an accepted word w of m>=ceil(N/2) events, put ell=floor(N/4) and
q=ell-2. The first ell recorded labels define distance before propagation.
Use a mathematical reference: restrict the **formation** Hamiltonian only
at words shorter than this prefix, retain every allowed extension and
receiver state, and tensor it with a neutral probe on a half-line with
hopping -J, J=1/4. Its initial probe is (|0>+|1>)/sqrt(2). This reference
factorizes geometry/receiver and matter; the physical Hamiltonian does not.

Embed both into a common direct-sum Hilbert space, adding dummy reference
probe positions after ell only for comparison. On their common word-prefix
states use identical configuration-preserving diagonal/receiver evolution.
Extend this selfadjoint diagonal arbitrarily but identically elsewhere.
Physical changing interactions have norm <=nu=2.1; reference interactions
have norm <=g+2J<=1.6. In their interaction-picture expansions:

* a probe starting at position 0 or 1 requires at least ell-1 hops to
  reach a dummy position or a changed prefix-boundary probe edge;
* erasing the fixed prefix requires at least m-ell+1 end pops;
* a bare append/pop requires the probe at the head, its parent, or a free
  event. From the initial pair these require at least m-2 changing moves.

Consequently every Dyson word of length <q agrees on the initial state.
The diagonal/receiver evolution cannot change these word or probe counts.
The reference restriction is used only in a proved comparison; it does not
remove any physical virtual geometry state. The **two complete tail sums**
give, uniformly over the actual conditional receiver and disorder,

    ||U_physical(s)psi-U_reference(s)psi||
      <=2 tail_q(nu*|s|)
      <=2 exp(-3q/4),  |s|<=q/(4e*nu).                 (B1)

The ideal probe tail at positions >=ell is at most
tail_(ell-1)(2J|s|)^2<=exp(-3q/2). Geometry remains unrestricted outside
the protected prefix in the reference and wholly unrestricted physically.
Processes that leave and later restore the prefix are included in (B1);
this is not merely a bound on inconsistent final records.

## 2. An elementary half-line ballistic lower bound

Let T|x>=|x+1>, X|x>=x|x>, H_p=-J(T+T^dagger). Then

    V=i[H_p,X]=iJ(T-T^dagger), ||V||<=2J,
    d<V>/ds=2J^2 Pr(X=0)>=0.                           (B2)

For the real initial adjacent pair, <V>(0)=0 and Pr(X=0)=1/2.
The crude bound |d Pr(X=0)/ds|<=4J gives Pr(X=0)>=1/4 for
0<=s<=s_0=1/(16J). Thus <V>(s)>=J/32 for s>=s_0 and

    <X>(s)>=Js/64,                  s>=1/(8J).

The exact Heisenberg identity for X and ||V||<=2J imply

    <X^2>(s)<=(1/sqrt(2)+2Js)^2<=9J^2s^2,  s>=1/J.

Paley-Zygmund applied only to the actual position distribution gives

    Pr_half-line[X>=ceil(Js/128)]>=1/147456,  s>=1/J.   (B3)

This is an extremely weak constant but independent of s and N. The
half-line boundary force in (B2), not a favorable launch momentum,
produces the positive outward velocity. Finite-support initial data
justify the unbounded-X identities by cutoff and the bounded velocity.

## 3. Full-H arrival in the recorded metric

Let d(s)=ceil(Js/128), and count the particle only on recorded prefix
positions d(s),...,ell-1. This is a relational receiver set fixed before
response is measured. Equations (B1)-(B3) imply

    Pr_physical[arrival in that set]
      >=1/147456-4 exp(-3q/4)-exp(-3q/2).

For q>=32 this is >1/294912; 5 exp(-24)<5/2^24<1/294912 already suffices.
For every N>=2^28 and

    1/J<=s<=q/(4e*nu),

formation acceptance is >=1/2, hence the complete original trial has

    Pr(accepted and matter arrives at recorded distance >=d(s))
      >=1/589824.                                      (B4)

At s=q/(4e*nu), both time and d(s) grow linearly in N. The existing
S4 upper front simultaneously excludes appreciable arbitrarily faster
transport across that same prefix. This is a nonvanishing injected-matter
transport lower bound at growing distance, alongside prefix persistence;
it goes beyond a nearby response alone.

This is **not** a lower bound on phase response chi_j at a specified distant
single receiver. The detected event is an explicitly declared receiver set
and an injection-versus-vacuum experiment, with success probability paid.
Quantum mechanics and this one-dimensional growth/transport grammar remain
inputs. No universal relativistic speed, selected dimension, common causal
structure for other species, gravity or distinctive empirical prediction
follows from (B4).
