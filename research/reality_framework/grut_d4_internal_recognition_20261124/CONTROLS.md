# D4 control and resource accounting

The test fault changes an ordinary label to another ordinary label. The
same fixed H keeps acting. At no point does the Hamiltonian receive the
fault time, a designation of the failed site, or a copy of the historical
answer outside the two memories written by its own dynamics. The isolated
pulse source and conventional time label describe the experiment. The
historical label is used to score ensemble probabilities.

| Case | State immediately after the fault on a clean written triple | What the Hamiltonian can do |
|---|---|---|
| Undamaged | Parities 00, baths vacuum | There is no repair-bath emission or change of S. The writer alone can revert to blank. |
| S changes | Parities 11, historical majority retained | The two drains receive two parity excitations in successful repair. |
| M1 changes | Parities 10, historical majority retained | The first drain can take one parity excitation. |
| M2 changes | Parities 01, historical majority retained | The second drain can take one parity excitation. |
| Two sites change together | One outlier relative to the opposite majority | Parity detects disagreement. The repair heads toward the opposite codeword. The historical joint score is exactly zero after nominal writing. |
| All three change together | Opposite unanimous word, parities 00 | No discrepancy is detected. The original historical color cannot return under H. |

The last two conclusions hold for the full writer-on H. It preserves
the majority coordinate during every fault-free interval. The coordinate
extends to blank components by retaining the S value. A clean post-write
state has no blank amplitude at nominal timing.

A two-site flip from a clean codeword is physically the same as a
single-site flip of the remaining site from the opposite codeword.
Thus the one-fault theorem also gives a lower bound above 99.9865% for
the **incorrect** unanimous codeword after a simultaneous pair and a
registered wait. The controller repairs the diagnosed minority even when
the historical majority has been lost. A triple flip leaves both drains
vacuum with the opposite source color. Its opposite ready codeword has
probability `cos^2(kappa t)` after elapsed t.

The crucial ambiguous pair is history 0 with an S flip and history 1
with M1 and M2 flipped together. Both give exactly
`|1,0,0>|1_W>|vac_B1,vac_B2>`, up to an irrelevant phase. There is no
physical difference for a detector, energy measurement or later repair
to use. The guarantee for separated faults depends on the bounds on the
residual parity at each later pulse. The proof includes the possibility
that a new fault arrives before the prior discrepancy has completely
disappeared.

## What remains in the finite reservoirs

For each i, the exact conserved charge between faults is `Q_i=p_i+N_i`.
From a clean ready word a single error raises the charge by its bit mask.
When repair succeeds, p_i is zero and N_i carries that charge. No bath
reset, erasure of an older syndrome or fresh external cold spin occurs
between the registered pulses. The supremum of possible total charge
after an ordered sequence is `q=sum(|mask|)`, with each bath at most the
sum of its incoming parity bits. The proof also bounds the probability
of the correct codeword with **both** expected counts present.

Cold preparation uses 4608 vacuum spins in the two chains, one ready
writer spin and two blank memory qutrits. The full evolution with fixed
test pulses preserves global von Neumann entropy. It starts at one bit
for a fair historical mixture and pure other degrees of freedom. When
faults are randomized, their source carries additional preparation or
entropy, depending on how that source is modeled. The finite hardware
stores the error charges as excitations when repair succeeds. A terminal
count record of one or many faults disperses spatially. Preparation of
vacuum, storage capacity and recurrent return of excitations are all
physical inputs or consequences. [INFORMATION_BOUND.md](INFORMATION_BOUND.md)

## Energy and timing

At nominal writing, the first unflagged S fault injects mean work
`2 Delta=8g`. Either first memory fault injects `Delta=4g`. No energy is
borrowed at the repair interface. Each subsequent pulse has exact work
`Tr[rho^-(F^† H F-H)]` with rho^- the pre-pulse state. All inter-pulse
energy changes are zero in the full system. The fault apparatus is
external. Its pulses are label blind. The asymmetric first-fault work
comes from the two chosen comparisons with S and does not discriminate
between the two historical colors. The state collision has equal full
energy distributions.

The finite XX hopping spectrum in charge sector Q_i lies inside
[-2g Q_i,2g Q_i]. The writer exchange has norm kappa and Q_W=1
throughout. Thus for total net pulse work W and a proof that the joint
historical-codeword and expected-count event has probability at least
p, we have

$$
(\Delta-2g)pq-\kappa\le W\le(\Delta+2g)q+\kappa .
$$

The same event implies a mean terminal total bath occupancy between
pq and q. For four S faults, q=8 and the conservative bound guarantees
more than 7.9355 stored quanta on average. The corresponding
conservative net work interval is approximately [15.8711g,48.00002g].
These are resource bounds for the full system. No claim of a measured
per-pulse work value in the large regime is made.

All repair interactions are always on. There is no gate scheduler or
clock that tells the model when to detect the discrepancy. The laboratory
parameter t, the writer's acquisition window and external fault pulses
are still supplied. The model has no energy-storing implementation of
the pulse source. It addresses whether recognition and repair after a
fault can proceed under fixed microscopic dynamics.

## Recurrence and transient changes

With vacuum correction chains and no faults, both charges stay zero and
there is zero false repair-bath activation. The writer can coherently
return both memories to blank. After the maximum registered duration
576/g this has probability about 7.7246e-5 at nominal writing. At a delay
pi/(2 kappa), reversion is certain in that undamaged control.

When both correction chains have one spin, a single M1 fault has an
exact recurrence at `2 pi/sqrt(g^2+kappa^2)`, where the discrepancy
returns with probability one. If that bath is instead initialized with
its interface spin excited at the beginning of a diagnostic, it can
create a discrepancy in an intact triple with probability
`g^2/(g^2+kappa^2)` at `pi/(2 sqrt(g^2+kappa^2))`. This resource stress
test shows why the vacuum condition matters. It does not enter the
vacuum false-activation score.

Within a repair, lowering one parity can move the outlier between
sites. For an isolated S fault in the ready comparison, each parity
has the same single-source occupation probability p(t). Each initially
correct memory is temporarily wrong with probability p(t)(1-p(t)).
This can reach one quarter. The uniform bounds score one terminal
read after adequate time, not an uninterrupted record or continuous
absence of false changes. [PROOF.md](PROOF.md)
