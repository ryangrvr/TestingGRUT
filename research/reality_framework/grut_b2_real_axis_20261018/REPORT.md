# B2-RA: actual-source scaling gate passes at author-derived status

## Decision

For the UNCHANGED nominal endpoint-source Hamiltonian on the already
frozen family L in 24N, L>=2016, the new actual-source argument gives

| Event | Uniform bound | Range |
|---|---:|---|
| Original one-source co-pair arrival A_L | >36% | Every t in I_L |
| Physical two-source joint output | >12.852225% | Every t in I_L=[L/(8J),L/(6J)] |
| Physical two-pair retention q=2 | >43.4281% | Every 2/J<=t<=L/(6J) |
| Successful decoded two-qubit/reference channel | Exact identity, success paid separately | Fixed endpoint law and same output event |

This passes the EXISTING 10% joint output, 25% retention and 1%
conditional information-error gate. Status: **author-derived, internally
challenged, exact constants checked; independent proof review pending**.
The old eight-event B2 endpoint remains failed. No law, parameter,
initial state, threshold or operating window was changed to obtain this.

## Load-bearing argument

REAL_AXIS.md factors the actual open-band zero out of the infinite-range
closed-channel kernel. The remaining factor has Wiener norm error
<=0.003482. A weighted ENTRYWISE boundary estimate gives a Neumann
constant <0.002542, including both thresholds. The exact fuel denominator
differs from -J exp(-ik) by less than 0.009192 J, so it cannot vanish
on the band. This is an actual localized inverse, not an inference from
the earlier trace-norm bound.

The actual source density in k lies between 24/25 and 26/25 times
(2/pi)sin^2 k. All unclassified source spectral weight is at most 1/25
and remains paid. Keeping the actual dispersion phase gives
|G(t)|<=2/25+1/(J|t|), hence <=29/50 whenever J|t|>=2.

A central ac spectral packet has weight >4/5. Its complete co-pair
kernel retains both traveling phases, the n+2 fuel-site offset, hard-cut
endpoint terms and a GLOBAL local remainder. Positive-word derivative
majorants prove the single velocity maximum and the amplitude-variation
estimate. In the frozen receiver/time region, the full missed central
squared norm is <3/20 (evaluated conservative upper 0.141023).
The coherent orthogonal-component inequality then gives A_L>9/25.

Finally B2-SO's existing same-space comparison transfers these source
results to the finite two-source event. It pays collisions, reflections,
two-source interference and both fuel returns. The ONE-source probability
is squared in the joint reconstruction. OPERATING_PROOF.md contains the
complete chain and exact numerical margins.

## Verification and challenge

The proof/code was frozen at
5e77eb0e9d4d950db68d5682069e95801d00011b before the arithmetic runs.
One primary and one identical-source replay passed all 110 exact
Fraction/identity checks; outputs were byte-identical. Combined child
wall time was 0.064716 s; peak RSS was 9.25 MiB. No dynamics, spectral
quadrature, time grid, large simulation or parameter search ran.

The controls enclose rational constants and square roots exactly. Display
decimals are not used for decisions. They do not validate the all-size
mathematical argument or the inherited B2-SO comparison by themselves.
CHALLENGE.md discloses the shared-context AI review, its objections and
the resulting clarifications. That review is not external validation.

## Assumption price and consequence

This is a dilute two-resource model on prepared fixed recorded geometry.
The relational generator rule, nominal couplings, fuel preparation,
quantum dynamics and external time remain supplied. Formation, moving
geometry, extensive density, dimensional selection and empirical
discrimination are not established here. Nothing here derives gravity
or supports a buoyancy interpretation of gravitation.

The specific unresolved SOURCE estimates that blocked B2-SO through
B2-BB now have an explicit author-derived proof route through the actual
kernel. This warrants independent review of that chain and, only after
its scrutiny, an explicitly proved composition with moving geometry.
It does not reopen the closed numerical branch, rescue the old endpoint,
or authorize another architecture automatically.

Read: OPERATING_PROOF.md for the result; REAL_AXIS.md and
SCATTERING_BOUNDS.md for the load-bearing estimates; REVIEW_PACKET.md
for the shortest route to falsifying them.
