# GRUT-P1 — charge-balanced preparation with explicit work

**Result: a controlled preparation regime, not an autonomous populated phase.**
P1 supplies a post-recording work protocol that loads local fuel, converts
it reversibly into opposite-charge matter, and operates the resulting
finite-density state with positive electric energy and an accounted empty
buffer. It has no root countercharge. The derivations are author-derived;
independent proof review is pending. Original-trial asymptotics additionally
depend on the archived S9 theorem, whose status is unchanged.

The charter was published at `44f8e82`; the analytical construction and
control sources at `44b2748`, before controls ran. The bounded suite passes
after one disclosed factorial-evaluation correction. D1 at `830746b` and
every frozen earlier result are unchanged. No closed formation solver was
reopened and no large dynamics were run.

## The concrete law

Two commuting hard-core species have charges +1 and -1 and positive unit
mass energies. Links have exactly three electric levels -1,0,1, with
Hamiltonian +sum E^2, fixed coefficient one. Fuel qubits have energy
3+Delta_e, |Delta_e|<=1/2. All new registers initially are vacuum.

After word formation at T_f=4N^2, accept a recorded word of length m>=N/2.
Select K=N/32 disjoint bonds in its first N/16 positions. A frequency-three
drive of duration pi/2 loads their fuel, then a reversible conversion pulse
of duration pi/2 creates + at a bond's parent, - at its child and one unit
of oriented link field. Both amplitudes may vary independently in [0.9,1.1]
without retiming. Reverse conversion and fuel return are present.

At operation start, switch both conversion directions off together and
switch on nonzero hopping j_+=1/64, j_-=1/96 and kappa=1/4 local density
interactions. Geometry and the full conditional receiver stay on throughout.
All switches have explicit work accounts. Matter-number persistence relies
on this supplied control; keeping conversion on is not proved stable.

Gauss, flux support, matter support and storage/receiver charge constraints
are preserved by explicit ambient-register operators. Opposite charges
require a new empty-head pop gate: otherwise a neutral co-located pair could
be left on a free event. This is a stated model change, not a theorem
inherited from positive-charge D1. The universal relational support rule
and finite electric capacity are also supplied.

## Earned estimates and their dependencies

N is a multiple of 64, N>=256. Let e_p<=exp(-N/4). Exact local Rabi algebra
gives each reference bond a pair probability >=49/64, uniformly over all
declared independent fuel detunings and pulse-amplitude variations. The
moving-geometry Duhamel comparison gives the full-state preparation error
e_p, including receiver/storage excursions and returns.

    Pr(P>=N/32 | accepted record)
       >=1-exp(-289N/65536)-2exp(-N/4).

No particle-number selection or favorable buffer measurement is performed.
The same controls place every reference particle inside b=N/16-1, leaving
an accounted empty end. Actual support errors are included in e_p.

If S9 is valid, its original neutral preparation has p_form>1/2 at
N>=2^28. Thus, per original trial with rejections counted zero,

    Pr(accepted AND P>=N/32)>1/4.

This statement still requires the inherited 2^62-mode receiver and formation
time 2^58 at threshold. P1 neither verifies S9 independently nor reduces
its resource threshold. The N=256 conditional formula lower bound 0.6766
is not a formation guarantee or simulation at that size.

Positive electric energy cannot use D1's matter-only local tree isometry.
Keeping finite electric links as explicit local factors gives a fixed-prefix
bond norm <2, local interaction bound J_0=4. A new common-space boundary
comparison, including extensive matter and backreaction, yields

    E_use(t)<=e_p+D(t),
    D(t)<=exp(-N/5)+3exp(-N/12), 0<=|t|<=N/65536.

Prefix-loss/storage-error probabilities are <=min(1,E_use^2). Local source
responses obey a connected propagation front plus 4||S||||O||E_use. Both
created species have nonzero nearby **bond-source** response >7.3e-8 in
magnitude at fixed t=10^-5, with nonzero source variance. These are not
single-event number phase pulses, and no growing-distance many-carrier
arrival lower bound or universal propagation speed has been derived.

## Energy is paid, and neutrality is not enough by itself

At exact resonance, three supplied work units per successful bond become
two mass units plus one positive electric unit. At detuning, residual fuel,
coherences and switching work remain in the ledger. Additional pulse/switch
system work is bounded by 6K; nonselective ideal word recording adds at most
1.1 to mean system energy. Pointer purity, memory and classical control
capability are separate supplied resources; no cyclic erasure is claimed.

Every reachable physical state has electric energy <=N-1. In the reference
operating regime its expectation is <=149K/96 for all times. A separated
single pair costs field energy proportional to its recorded separation.
This is a confining, capacity-restricted medium, not an unconfined plasma.

Zero net charge alone would not cure D1's scaling: K positive charges
followed by K negative charges cost K(2K^2+1)/3 with an unrestricted link.
P1's local preparation, fixed capacity and positive energy are substantive
new choices. The receiver's existing energy is not spent again as fuel.

## Executed checks, not claims of independent validation

| Control | Result |
|---|---|
| Independent 96-dimensional ambient two-site tensor, 12-dimensional physical block | Gauss commutators and physical entry differences exactly zero in this computation |
| 27 pulse parameter combinations | Rabi probability error <=1.12e-15; energy-ledger error <=3.56e-15 |
| Reverse pulse at resonance | Returns matter to fuel; residual matter probability <9e-32 |
| Four-site fixed path, two pumps, all physical number sectors, 272 states | Pure-state error against product solution 8.23e-15 |
| Same generated state, two nearby bond responses | Approximately -1.4842e-7 and -9.8945e-8; both exceed their magnitude lower bounds |
| Complete four-event geometry, two-mode diagnostic receiver, 21,920 states | Short-preparation error 3.286e-5 <=9.824e-4; short-operation error 9.034e-6 <=0.19093 |
| Storage-error channel retained | Nonzero probability 1.507e-11 in the operating control |
| Frozen D1 vacuum matrix comparison | Maximum entry difference 1.12e-16, from floating summation order |

The full-geometry controls use a two-mode receiver, not the physical
1,024-mode receiver at N=4; short pulse times are .02 each, not a new
primary preparation protocol. They test implementation and boundary
estimates, not macroscopic formation. All four-event physical matter-number
sectors, both fuel levels, all storage bits, fields and allowed receiver
occupations are retained. No wavefunction residual has been clipped.

The passing suite took 2.459 seconds and 117.39 MiB peak resident memory.
Ordinary double precision and agreement checks are not an outward-rounded
certificate or an external mathematical review. FIRST_RUN.txt preserves the
formula-overflow failure; parameters, times and acceptance gates did not change.

## Assessment

P1 connects recorded geometry to dynamically generated extensive matter in
one explicit **externally powered** preparation-and-use protocol. It resolves
the D1 root-countercharge/zero-electric-energy problem within a more restricted
charge-balanced candidate, at the price of fixed flux capacity, a new support
gate, classical feed-forward and conversion control. It does not establish an
autonomous matter-bearing phase, a selected spatial dimension, original CP-12
matter universality, gravity or an empirically distinctive GRUT prediction.

The next substantive choice is whether to implement this accounted controlled
medium or separately investigate an autonomous work reservoir and persistent
conversion law. That would be a new campaign, not a property to attribute
retrospectively to P1. First independently check the explicit proof obligations
in HANDOFF.md; another large solver campaign is not needed for this result.
