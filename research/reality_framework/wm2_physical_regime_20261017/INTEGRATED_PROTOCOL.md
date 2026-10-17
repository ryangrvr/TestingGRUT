# Paid formation, recording, injection and finite operational use

This combines the S9 analytic formation result with the unchanged S4
conditional theorem. It does not inherit S2's real-energy gap, or claim
selected dimension or original CP-12 charged-matter universality.

1. Prepare the original uniform-root neutral auxiliary seed, consistent
   storage and physical receiver vacuum, with the frozen g/disorder class.
   Evolve the complete closed Hamiltonian to a uniformly sampled time in
   [4N^2,8N^2]. For N>=2^28, accept every word m>=ceil(N/2); S9 gives
   p_acc>=1/2. No ordering is optimized. The large uncovered finite-size
   interval remains explicit.
2. Record the entire geometry word and whether it passed. Retain the actual
   conditional receiver/storage state. Rejects give zero probe response and
   retain their original probability. There is no field reset or favorable
   postselection inside the accepted event.
3. On acceptance inject S4's neutral probe (|w_0>+|w_1>)/sqrt(2), and
   pulse exp(-i theta n_(w_0)). Source number variance is exactly 1/4.
   Evolve the full probe Hamiltonian, including every geometry and storage
   transition, and sum all final records at each receiver.

## Operational guarantees, per original preparation

Let r_N=max(1,ceil(N/2)-2), a_N=max(1,floor(N/4)-2), J=epsilon=1/4,
nu=2.1. For probe duration s<=a_N/(4e nu), S4 gives storage-error
amplitude <=(epsilon/nu)exp(-3a_N/4), uniformly over the actual conditional
receiver state and every allowed independent detuning pattern. Its response
correction includes repeated complement visits and returns.

There is also an explicit extended-prefix guarantee. Before an accepted
word can lose its original first floor(N/4) events, at least a_N end pops
are required. The configuration-preserving receiver/diagonal Hamiltonian
does not change that word count; all other interactions have norm <=nu.
The interaction-picture Dyson series therefore bounds prefix-loss amplitude
by tail_(a_N)(nu*s)<=exp(-3a_N/4) in this window. This includes bare pops,
reattachments and all field states. Hence

    Pr(accepted, consistent, original N/4 prefix retained at s)
      >= (1/2)[1-(1+(epsilon/nu)^2)exp(-3a_N/2)].

The probability is a joint single-time guarantee for every allowed s,
not a continuously monitored survival probability. No repeated measurements
or Zeno locking are implied. The window grows as N, while formation takes
N^2. Permanent locking or an equilibrium phase is not established.

S4's full-H front bounds transmission beyond original core distance d by
tail_d(nu*s), with its conditional response comparison controlled as above.
At the nearby receiver, the success-weighted phase response obeys

    chi_trial(s) >= (1/2){J*s-[exp(2nu*s)-1-2nu*s]}.

At s=.005 this is >.0005139. Rejected trials are included with zero signal.
That phase-response estimate alone establishes nearby response and an upper
front. BALLISTIC_TRANSPORT.md supplies a separate injection-arrival lower
bound at growing recorded distance: for 1/J<=s<=a_N/(4e*nu), let
d(s)=ceil(Js/128) and count arrivals on prefix positions d(s),...,floor(N/4)-1.
Including formation acceptance, storage consistency and retained prefix,

    Pr(accepted, consistent, prefix retained, arrival at distance >=d(s))
      >=1/589824,       N>=2^28.

Indeed its conditional arrival bound before these last two exclusions is
1/147456-4exp(-3a_N/4)-exp(-3a_N/2). Subtract at most
[1+(epsilon/nu)^2]exp(-3a_N/2) for prefix/storage failure. For a_N>=32,
the total correction is <7exp(-24)<1/294912. Multiply by p_acc>=1/2.
At the largest allowed s, both distance and time are proportional to N.
This is arrival in a stated receiver set, not nonzero phase derivative at
each distant site, a speed common to other species, or a continuum light cone.

## Energy, record and resource account

In the relative formation convention, the seed has mean energy zero and
closed evolution conserves it. The physical receiver's positive frequencies
are offset by binding/storage energy through the original conserved charge;
no new fuel or external formation schedule is introduced.

Word recording commutes with D+dGamma(h) but erases off-diagonal assisted
coupling expectation B. For the nonselective measurement across **all**
records, Delta E_record=-<B>, |Delta E_record|<=g. This is an apparatus
energy ledger, not zero-cost recording. Conditional accepted energy need
not obey that same bound; both outcomes and their field energies remain
in the aggregate ledger. At most log_2(e*N!) geometry-label bits suffice,
apart from measurement hardware and later record processing.

Neutral probe injection is supplied and changes particle number. In S4's
relative single-probe Hamiltonian its initial pair has hopping energy -J;
any absolute probe creation/rest-energy convention and apparatus work are
additional supplied resources. The phase pulse has pair hopping-energy
change J(1-cos(theta)); it is an external intervention, not an autonomous
observer derivation.

At N=2^28 the declared receiver has L=2^62 modes and capacity N-1 quanta;
the formation window starts at 2^58 receiver-hopping time units. These
are mathematical scaling resources, not a practical experimental proposal.
The receiver adjacency is supplied receiver structure; the operational
matter path remains the separately generated, grammar-constrained word.
