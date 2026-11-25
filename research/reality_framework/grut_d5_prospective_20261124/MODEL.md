# GRUT-D5 finite prospective development

25 November 2026. A moving seed head grows and reconstructs a finite tape.
The two rule choices produce uniform or alternating output. The seed's
initial physical color is acquired as a record and used in the writes.

## Registers and preparation

N>=8 is a multiple of four. Each of N tape sites has a displayed cell A
and archive cells B,C, all three-level with symbols 0,1,e (blank).
All 3N cells start blank. One head starts at position 0 on a path
0 through L=256N. Its registers are:

| Register | Alphabet and preparation |
|---|---|
| Physical seed color h | 0 or 1, initially fair |
| Rule G | 0,1,e, initially fair on 0 and 1 independently of h |
| Acquired memory M | 0,1,e, initially e |
| Phase P | bit, initially 0 |
| Independent coin c | fair bit |
| Backup slots E_G,E_M | three levels each, initially e |
| Fault selector F | 0,1,2,3, with 0 meaning no diagnostic fault |
| Control selector K | intact, G inverted, G scrambled, G disabled, M inverted, M scrambled, M disabled |

The intended tape for initial h,g has entry `h xor [g and (i mod 2)]`
at site i, numbered from zero. The write rule at the moving head reads
`M xor [G and P]`, or requests blank if either G or M is e. Define T_r
as exchanging e and r for r in {0,1}, and T_e as identity. X exchanges
0 and 1 and fixes e. The initial h and g label scores and no extra
reference register is coupled. The original seed color remains present
in the head after the assay.

## Microscopic exchange

Head edge j connects nearby positions j and j+1. Set hbar=J=1 and
omega=4J/(L+1). The time-independent Hamiltonian in the conserved
one-head sector is

$$
H=\sum_{j=0}^{L-1}\frac{\omega}{2}\sqrt{(j+1)(L-j)}
\left(|j+1\rangle\langle j|\otimes U_j+
      |j\rangle\langle j+1|\otimes U_j^\dagger\right).
$$

Physically each route site has a vacuum or a head carrying its
registers. An exchange moves the head from one adjacent site into a
vacancy, with a matrix acting on its internal state and at most one
tape site. Each tape site is next to its successive visit edges on
separate rails. Idle edges route the head and supply waiting regions.
The route is a supplied bounded-degree interaction graph. Each term has
support on at most three sites and norm at most J. The head position
is the internal stage variable.

## Local updates and route

- Assay uses T_h on M and leaves h unchanged.
- Grow uses T_color on A, then toggles P.
- Rebuild uses T_color on fresh B, swaps B and A, then toggles P.
- Maintain toggles P first, uses T_color on fresh C, then swaps C and A.
  The different order keeps the phase aligned during the reverse walk.
- A fault edge flips the indicated A if its selector matches. At the
  central edge apply K after the fault. Inversion exchanges the two
  labels of the chosen G or M register, scrambling does that only for
  coin 1, and disabling swaps with the corresponding backup slot.
  All updates are reversible permutations. The writes read neither
  the fault selector nor the control selector.

Set k=N/4 and m=N/2. The only declared nonidentity links are:

| j, edge from j to j+1 | Update |
|---|---|
| 0 | Assay |
| 1+i, i=0,...,k-1 | Grow i |
| L/4 | F=1 flips A_0 |
| L/4+1+(i-k), i=k,...,m-1 | Grow i |
| L/2 | F=2 flips A_1, then apply K |
| L/2+1+m+i, i=0,...,N-1 | Rebuild i |
| 3L/4 | F=3 flips A_0 |
| 3L/4+1+i, i=0,...,N-1 | Maintain N-1-i |

All other links are identity updates. The m idle edges after the
central edge route the head back to the beginning. The entire path
and its couplings stay present through the experiment. The central rule
and record controls preserve dimensions, with backups and an unused
fair coin present in the intact arm. The flags select diagnostics in
the initial preparation. The complete target first occurs on the
ordered intact history at the last rebuild write, after the two early
diagnostic sites. The maintenance pass follows completion.

## Experiments and scores

T=pi/omega is the terminal time. Score the event of head at L and A
equal to the target from the original labels on all trials, including
misses. Authenticating faults also requires the flipped color in the
corresponding archive B_0, B_1 or C_0. The pre-control record and the
numbers of filled cells at both early stages are separate scores.

External disturbance trials set F=0 and K=intact. They apply one laboratory
X to A_0, A_1 or A_0 when sin^2(omega t/2) equals 1/8, 3/8 or 5/8,
respectively. They make no stage measurement. The active work is constant
on each corresponding idle band. Clock tails can occupy other stages
at the laboratory pulse. The exact ordered history and the laboratory
time probabilities are distinct statements.

The matched null retains all registers, the assay, diagnostics and
couplings but reads the original seed h in place of M for color
deposition. Thus it produces an acquired record and the same intact
patterns while leaving the record causally unused by the writes.
