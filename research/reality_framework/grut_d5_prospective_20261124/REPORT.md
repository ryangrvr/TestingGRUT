# GRUT-D5 prospective developmental identity

25 November 2026. A single finite, time-independent Hamiltonian starts
with a simple head and a wholly blank tape, acquires a record, and uses
it to build any of four ordered words. If an already written cell is
flipped when only a prefix exists, the system rebuilds the damaged prefix
and fills the unbuilt suffix. The whole word first appears after that
rebuild on the ordered internal path. The result uses a compact supplied
local rule, a routed head, and two prepared blank archive banks.

| Criterion | Verdict |
|---|---|
| Generative development | Exact. One rule bit selects uniform or alternating growth, with seed color selecting the phase. Every displayed cell starts blank. |
| Acquisition of developmental records | Exact on the ordered path. The assay writes the seed color into a blank memory before any deposit. Archived cells later record actual prior site values. |
| Causal contribution of memory | Exact within this Hamiltonian. Holding the damaged substrate fixed, scrambling the acquired value halves whole-word success. Inversion or disabling reduces it to zero. The pilot-null passes the intact case and is insensitive to the memory interventions. |
| Resumption following disturbance | Exact at the terminal time for internal early and late fault choices. A finite bound also covers one external unflagged flip at either early stage or after full construction. |
| Scaling and resource efficiency | Two data bits select N output symbols. The hardware, time, and archive capacity grow linearly in N. No claim of optimal resource cost. |
| Autonomous dynamics | The internal episode evolves under one fixed local Hamiltonian. The clock route, stage gates, and fault selectors are supplied. The external probes are separate disturbances. |
| Emergence versus supplied programming | Repeated local synthesis with a fixed rule and visit order. The structure is generated from that rule, seed and geometry, with no sitewise target copy. Collective selection of the rule or route remains unestablished. |
| Physical novelty | Standard reversible dynamics and engineered spin transfer. This is a constructive example. |

## What develops

The assay acquires the initial color h as M. The generative bit G
chooses whether the next-site color alternates. A parity bit is
updated by each local visit. The head's location, the acquired M,
the rule G and the displayed tape are physically separate. For
N=8, the first two possible disturbances occur with 2 and 4 filled
cells. At either early stage the damaged cell is recovered and
the remainder is built, which tests more than returning to a
completed prior configuration. The later fault tests recovery
of a completed tape. Replaced values remain in local qutrits.
[MODEL.md](MODEL.md)

## Exact and finite-time results

The head path is unitarily equivalent to a spin-L/2 rotation.
The ordered endpoint of a no-probe run arrives exactly at
T=pi/omega, and all internal intact cases score one. Both early
cases and the late case archive the actual flipped value. The
common energy distribution is analytically equal across all
internal selector arms, including masked and disabled records.
The matching null also shares that spectrum and preparation law.

The full quantum state has clock tails during laboratory
probes. We apply a physical flip to a tape cell without reading
the clock. The proof bounds failure of **terminal arrival, whole
word and flipped-value archive together** by less than
**7.32e-17** at N=8 and below that at larger sizes. This
conservative analytic bound holds over intervals of clock phase
around the three sampled times. The exact sample calculation
keeps all coherent clock contributions and records the remaining
tiny failure probabilities in log form. A distinct finite late
interval supports a persistence bound on the output. Recurrence
precludes indefinite persistence. [PROOF.md](PROOF.md) ·
[RESULTS.json](RESULTS.json)

Strict "never previously existed" here concerns the ordered
internal path. In the unobserved real-time state there is
small amplitude at other clock coordinates. We bound this
amplitude at the early probes and do not equate a
clock-conditioned path with a monitored monotone history.

## What the controls isolate

For the same damaged tape, the original whole-word scores are
1 with either input intact, 1/2 for either fair scramble, and
zero for inversion and disabling. A null copying the same M and
using h directly gives the same intact joint evolution, then
retains full success when M is manipulated. Observation of shape
and record correlation alone therefore cannot show that the
record participates in development. The intervention tests the
write rule's causal input. [CONTROLS.md](CONTROLS.md)

If four distinct intended future patterns are equally likely
and their complete accessible pre-action states are identical,
no present mechanism can succeed at matching more than 1/4 of
them. This follows by pulling back the mutually exclusive
whole-word effects to an input state-discrimination problem.
D5's seed and acquired record supply information about its
target. An independent future choice supplies no such information
in advance. [LIMITS.md](LIMITS.md)

## Resource and reproduction checks

At N=8,16,32,64 we tested 128 intact and fault trajectories,
224 internal control cases, and 96 external single-probe cases.
The internal controls retain dimensions and initial entropy.
Two label bits select the output and a blank memory acquires
one bit correlated with h. There are 256N+1 head nodes and 2N
blank archive qutrits. The transfer time grows linearly in N
under the fixed maximum link scale. The exact internal energy
law, zero mean work of the prescribed real tape flip, and the
resource costs of its external actuator are distinguished.
[RESOURCES.md](RESOURCES.md)

The design was recorded before production. A separate internal
reviewer independently built the reachable physical Hamiltonian
and directly propagated through 27 external pulses, including
small-clock cases where the intended-stage gate fails. The
coherent endpoint amplitudes agree within 7.87e-15. The reviewer
also directly propagated 18 small internal control cases. The
main solver matches the independent stage and joint scores
within 3.34e-15. The physical-model comparison is separate
from the analytic large-clock certificate. External
RAI review remains open. [CRITIC.md](CRITIC.md) · [README.md](README.md)
