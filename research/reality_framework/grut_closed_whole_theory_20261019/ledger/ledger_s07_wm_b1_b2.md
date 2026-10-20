# Campaign ledger entry — s07_wm_b1_b2: WM / B1 / B2 matter and moving channels

Sector key: `s07_wm_b1_b2`. Repository `/home/user/TestingGRUT`, branch
`codex/grut-closure-fork-20261019`. The handoff commit `96a7e7e` is an ancestor
of the checked-out HEAD `da8f925`; `git diff --stat 96a7e7e HEAD` over every
source directory used below is empty, so the sources read here are byte-identical
to the handoff state. All paths are relative to the repository root;
`RF` = `research/reality_framework`. Nothing tracked was modified; the only
execution was the closeout's own read-only rational-arithmetic script, run from
the scratchpad (section 4.9).

---

## 1. Question asked in this sector

Closed-whole parent question as it was asked here: in a **finite closed quantum
model** whose "bath" for the matter subsystem is the rest of the same model — the
recorded geometry word, storage bits, tail registers and the retained
`64N^2`-mode conditional bosonic receiver, never reset — can one **supplied,
time-independent Hamiltonian** (i) start from an *accepted prepared record*,
(ii) convert fuel into matter, (iii) *deliver* it to a readout region,
(iv) *retain* it throughout a stated window, and (v) make an encoded input
*recoverable*; and does anything in that force a law beyond ordinary quantum
mechanics? (model: `RF/grut_b2_feasibility_closeout_20261018/SPECIALIST_REVIEW.md:15-56`;
`RF/grut_b2_consolidation_20261018/MOVING_GEOMETRY.md:7-41`; "no unmodeled energy
sink or controller is used": `RF/grut_b1_neutral_escape_20261018/REPORT.md:13-14`.)

The B2 closeout sub-question: does any **currently available quantitative route**
establish the moving-geometry thresholds — output `>=1/10`, retention `>=1/4`,
conditional information trace distance `<=1/100` — within the authorized compute
envelope? (`SPECIALIST_REVIEW.md:85-88`; `STOP_DECISION.md:6-13`.)

## 2. Handoff row, as read

`RF/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md:120`:

> | WM/B1/B2 matter and moving channels | Keep preparation acceptance,
> recoverability, delivery and retention distinct. The B2 closeout is a
> feasibility stop within a specified computational envelope, not a physical
> impossibility proof. Results from incompatible state spaces are not
> composable. | [closeout disclosure](../grut_b2_feasibility_closeout_20261018/REVIEW_STATUS.md),
> original WM/B reports via the continuation/index documents |

Checked against primary sources in section 13.

## 3. Assumptions / model class

- **Supplied law and constants.** Hard-core plus/minus species, three-level
  electric links `E in {-1,0,1}` with energy `+sum E^2` and Gauss law, neutral
  fuel slots, a `64N^2`-mode open-chain bosonic receiver with nearest-mode
  element `-1`, `g in [0.9,1.1]`, `delta_i = bar_delta + epsilon_i`,
  `|bar_delta|<=1/2`, `|epsilon_i|<=1/4`; `jp=1/64`, `jm=1/96`, `A=13/46080`,
  `J=lambda=1/3840`, `Omega=2-2A`; "These are the actual source/B1 choices,
  including supplied locality and resonance." (`SPECIALIST_REVIEW.md:17-56`;
  `RF/grut_b2_consolidation_20261018/THEOREM.md:9-33`.)
- **Preparation is supplied, not autonomous.** "Recording and loading are
  supplied." (`MOVING_GEOMETRY.md:41`); "This loading and the record are
  supplied, not autonomous formation claims." (`SPECIALIST_REVIEW.md:70`);
  P1 is "a controlled preparation regime, not an autonomous populated phase"
  (`RF/grut_p1_preparation_20261017/REPORT.md:3`).
- **Readout, decoder and thresholds are supplied.** Output window
  `[L/(8J),L/(6J)]`, retention window `[2/J,L/(6J)]`, thresholds `1/10`,
  `1/4`, `1/100` (`SPECIALIST_REVIEW.md:75-88`).
- **Status of every derivation.** Author-derived, independent mathematical
  review pending; internal "Challenger" shares ancestry and is "not independent
  mathematical review" (`REVIEW_STATUS.md:3-7`; `RF/grut_b2_consolidation_20261018/CHALLENGE.md:1-6`;
  `RF/grut_b1_neutral_escape_20261018/CHALLENGE_AND_INTEGRATION.md:3-5`).
- **Computational envelope for the closeout.** "at most 25 seconds per child,
  60 combined child seconds, 1 GiB memory and 25 MiB new files"
  (`FEASIBILITY.md:84-86`).
- **Acceptance threshold inherited from S9.** `N>=2^28, L=64N^2, 4N^2<=t<=8N^2
  => Pr(m>=ceil(N/2))>1/2` (`RF/wm2_physical_regime_20261017/REPORT.md:7-8`),
  "a theorem derivation about a supplied-resource quantum path model, not
  external mathematical validation or a distinctive prediction about nature"
  (same file `:10-11`).

## 4. Executed evidence, with sources and exact numbers

### 4.1 WM2 line (one-line verdicts from `head` skims)

| Packet | One-line verdict (as printed) | Source |
|---|---|---|
| WM2-S1 defect assembly | "Its closed seeded dynamics fail appreciable formation." "The full WM2 success theorem has therefore not been achieved." | `RF/wm2_defect_assembly_20261016/REPORT.md:6-10` |
| WM2-S2 fuel lock | resolves S1 "at a supplied exact resonance"; "rejects a robust untuned-phase claim" | `RF/wm2_fuel_lock_20261016/REPORT.md:3-10` |
| WM2-S3 detuning receiver | "a finite-receiver stability interval, with a spectral limitation"; "does not establish a generic spatial phase ... or a distinctive empirical prediction" | `RF/wm2_detuning_receiver_20261016/REPORT.md:1-11` |
| WM2-S4 prepared disorder | "prepared bulk protection, with disordered formation still open"; "does not establish large-size disordered formation" | `RF/wm2_prepared_disorder_20261016/REPORT.md:1-8` |
| WM2-S5 history formation | five-event formation lower bounds "above 77% and 80%"; "six-event target remains unresolved under the declared compute budget" | `RF/wm2_history_formation_20261016/REPORT.md:1-9` |
| WM2-S6 conditioned receiver | "formation frontier unchanged"; "Six-event formation at T=144 remains unresolved." | `RF/wm2_conditioned_receiver_20261016/REPORT.md:1-8` |
| WM2-S7 positive residual | "five-event recovery passes; the six-event target remains computationally unresolved"; "the stop is not a negative formation result" | `RF/wm2_positive_residual_20261017/REPORT.md:1-8` |
| WM2-S8 Schmidt receiver | "both original S8 recovery components stop before the target"; "numerical-method outcomes, not evidence against physical formation" | `RF/wm2_schmidt_receiver_20261017/REPORT.md:1-7` |
| WM2-S9 physical regime | "an analytic operational regime is derived; the bounded finite recovery fails"; threshold `N>=2^28` | `RF/wm2_physical_regime_20261017/REPORT.md:1-11` |
| WM2 proof/matter challenge | "shared geometry is not yet enforced"; S9 self-check "found no failed load-bearing estimate" but "This is the same Builder checking its work" | `RF/wm2_proof_matter_challenge_20261017/REPORT.md:1-12` |
| WM2-UC1 charge constraint | "Primary architecture verdict: negative in the frozen class U."; "universal locality still fails" | `RF/wm2_charge_constraint_20261017/REPORT.md:1-8` |
| WM2-IG1 generator classification | "physical observables do not select the interaction law"; locality rule "is an input, not an emergence result" | `RF/wm2_generator_classification_20261017/REPORT.md:1-10` |
| WM2-D1 many carrier | "a conditional many-carrier regime"; "not an unchanged-law result, independent validation, or a completed theory of reality"; `P_N=N/16, N_A=N_B=N/32, 0<=|t|<=N/4096`, error `<=exp(-N/5)+3exp(-N/17)` | `RF/wm2_many_carrier_20261017/REPORT.md:1-12` |

### 4.2 P1 — preparation acceptance with explicit supplied work

- Protocol: after word formation at `T_f=4N^2`, accept a recorded word of length
  `m>=N/2`; select `K=N/32` bonds in the first `N/16` positions; drive of
  duration `pi/2` then conversion pulse `pi/2`; amplitudes in `[0.9,1.1]`
  (`RF/grut_p1_preparation_20261017/REPORT.md:24-29`).
- Earned bound: `Pr(P>=N/32 | accepted record) >= 1-exp(-289N/65536)-2exp(-N/4)`,
  `N` multiple of 64, `N>=256`, per-bond pair probability `>=49/64`
  (`REPORT.md:46-52`). With S9: `Pr(accepted AND P>=N/32)>1/4` per original
  trial, still requiring "the inherited 2^62-mode receiver and formation time
  2^58 at threshold" (`REPORT.md:59-67`). Use-window comparison
  `E_use(t)<=e_p+D(t)`, `D(t)<=exp(-N/5)+3exp(-N/12)`, `0<=|t|<=N/65536`
  (`REPORT.md:73-75`).
- Controls (double precision, not certificates): 27 pulse combinations, Rabi
  error `<=1.12e-15`; four-event two-mode-receiver control with 21,920 states,
  short-preparation error `3.286e-5<=9.824e-4`; suite `2.459 seconds`,
  `117.39 MiB` (`REPORT.md:105-114, 123`). "The full-geometry controls use a
  two-mode receiver, not the physical 1,024-mode receiver at N=4"
  (`REPORT.md:116-117`).
- Failure retained: first run overflowed (`OverflowError: int too large to
  convert to float`) and was corrected by log-gamma evaluation with a floor
  `exp(-700)`; "No model coefficient, time, size, error gate, parameter class
  or resource budget was changed." (`RF/grut_p1_preparation_20261017/FIRST_RUN.txt:1-20`.)
- Verdict: "externally powered"; "does not establish an autonomous
  matter-bearing phase, a selected spatial dimension, original CP-12 matter
  universality, gravity or an empirically distinctive GRUT prediction"
  (`REPORT.md:130-136`).

### 4.3 A1 — autonomous retention gate refuted (one permitted point)

- Gate: `Pr(P(t)>=N/64)>=1/4` at **every** `t in [pi/2, pi/2+N/65536]`,
  `N>=2^18` multiples of 64 (`RF/grut_a1_fixed_conversion_20261017/REPORT.md:14-16`).
- Result at `Omega=13/4`, `lambda=1`: `Pr(P(pi/2)>=K)>0.946-2D_A(pi/2)`,
  `Pr(P(pi)>=K)<0.215+2D_A(pi)<1/4`, `D_A<=exp(-N/5)+exp(-N/12)`
  (`REPORT.md:36-44`). "One point refutes uniformity on the whole class; it
  does not settle other parameter choices or observation protocols."
  (`REPORT.md:46-47`.) "No energy or physical system disappears when the
  declared matter count falls." (`REPORT.md:53`.)
- Finite controls: one cell `0.999999995` at `pi/2`, `0.00065247` at `pi`;
  eight events, all 773 states, `~1` at `pi/2`, `0.00345813` at `pi`, energy
  conserved within `5.4e-15`; suite `0.9872 s, 121.41 MiB` (`REPORT.md:78-88`).
- Separate surviving spectral result: lowest mobile co-located-pair band with
  hopping `t_0=1/3840`, dispersion error `<8.1e-6`; adjacent one-pair initial
  weight in that band `<=1/529`; other bands `-1/23040`, `1/4608`
  (`REPORT.md:96-107`).
- Methods disclaimer: "Rabi conversion, local perturbation control and
  strong-coupling bound-pair bands—are conventional quantum physics ... not an
  earned new law of nature" (`REPORT.md:118-120`).

### 4.4 B1 — a conditional one-composite delivery window

- For `j_+=1/64, j_-=1/96, kappa=1/4, J=1/3840`: matter present for every
  `t in [7680,15360]` with probability `>0.570975`; at `t=15360`, intact prefix
  and co-pair at recorded distance `>=4` with probability `>3/5`; decode on
  arrival an "Exact conditional identity channel; missing trials are erasures"
  (`RF/grut_b1_neutral_escape_20261018/REPORT.md:16-24`).
- Fixed-path theorem uses `>=64` events; the moving extension uses
  `m>=N/2`, `N>=2^28`, is "an author-derived argument requiring independent
  review, not a simulation of that enormous system"; `2^62` receiver modes
  "are not reduced by B1" (`REPORT.md:26-33`).
- Executed: exact matrices on 2, 8, 16 events; two resources on 8 events with
  all 773 states; 16-event control `0.996551` matter, `0.879595` arrival at
  `t=15360`; two-event control "returns below threshold"; worker `0.312144 s`,
  `71.59 MiB` (`REPORT.md:37-49`). Two-event matter probability `0.182092` at
  `theta=4` is "an explicit counterexample to inferring retention from the
  local conversion cell alone" (`FAILURE_MODES.md:29-32`). Two-fuel 773-state
  matrix: `p(two pairs)=0.751352`, mean pairs `1.735023` at `theta=4`, "not a
  tensor product of isolated one-resource channels" (`FAILURE_MODES.md:43-45`).
- Boundary: "Every finite system has recurrence. No permanent locking,
  irreversible matter creation, asymptotic exponential lifetime, or
  parametrically large ratio tau_loss/tau_transport follows."
  (`FAILURE_MODES.md:36-39`.) "The same leading predictions follow from a
  conventional tight-binding channel. Empirical distinction remains absent."
  (`CHALLENGE_AND_INTEGRATION.md:19`.) Supplied: "New neutral conversion,
  near-band fuel energy, weak contact, endpoint-degree counterterm and color
  symmetry" (`REPORT.md:55-57`).

### 4.5 B2 — original eight-event endpoint failure and its diagnosis

- Sector account at `theta=6`: no matter, both fuels excited `48.89%`; one
  matter-pair quantum `43.14%`; two matter-pair quanta, fuels empty `7.97%`;
  prescribed accepted outgoing event `1.93%`; "70.46% was EACH fuel's
  occupation, not the joint both-full probability." "No listed probabilities
  are multiplied." (`RF/grut_b2_diagnosis_20261018/REPORT.md:12-20`.)
- Cause: from `theta=2` to `6` each endpoint receives `0.690129` resource back
  across its boundary bond while fuel increases by `0.703486`; remainder
  `0.013356` (`REPORT.md:25-30`). Spectral weight in 38 low modes; high-band
  weight `1.95e-11`, bounded `<8.65e-7`; interference adds `0.3764` to
  double-fuel return and suppresses accepted output by `0.09356`
  (`REPORT.md:33-40`). Conditional fixture `|gamma|` about `0.99999583`,
  "nearly unitary, not a large demonstrated loss", but these returns "do not
  rescue its 1.93% accepted output" (`REPORT.md:51-56`).
- Verdict: "retain escape as a conditional transport mechanism"; "B2 still
  fails its original endpoint task" (`REPORT.md:3-8`); "No distinct empirical
  prediction follows." (`ARCHITECTURAL_VERDICT.md:39`.)

### 4.6 B2-RA fixed-path pass and why it does not compose with the moving law

- B2-RA (commit `38fa6a3`): for `L in 24N, L>=2016`, `J=1/3840`: joint co-pair
  output `>12.852225%` on every `t in [L/(8J),L/(6J)]`; two-pair retention
  `>43.4281%` on every `t in [2/J,L/(6J)]`; conditional successful-output
  channel exact. "The geometry is prepared and FIXED, not dynamically retained
  by this theorem." Status: author-derived, review pending
  (`RF/grut_b2_consolidation_20261018/REPORT.md:5-8, 16-31`;
  `RF/grut_b2_real_axis_20261018/REPORT.md:5-19`). The 39 exact checks used a
  67-dimensional four-event scalar matrix, "an algebra control, not an operating
  simulation"; wall `0.109865473 seconds`, `10,496 KiB` (`REPORT.md:70-78`).
- Non-composition: the inherited comparison needs an EMPTY RESOURCE BUFFER
  `d_buffer=ell-2-b>=0`; for the exact B2 reference `d_buffer=-1`, "The
  support hypothesis fails BEFORE evolution"; moving the weight origin gives
  `sqrt(W(0)) exp(-d_buffer/2)>=sqrt(2e)` so "The exponential suppression
  therefore cancels" (`MOVING_GEOMETRY.md:72-95`). Exact local witnesses:
  `<target|(H_full^2-H_ref^2)|initial>=J jp=1/245760`, minus analogue
  `J jm=1/368640`, diagonal differs by `+A`; assisted append element
  `-g/N sqrt(n_contact+1)` even with fuel on the old head
  (`MOVING_GEOMETRY.md:112-145`). Endpoint order identity not inherited
  (`:147-159`). Horizon: `T_L=L/(6J)=640L`, `1,290,240` at `L=2016`, "84 times"
  B1's `15360` (`:161-169`; `REPORT.md:54-55`).
- The missing estimate is Eq. (3): `||[exp(-itH_full)-exp(-itH_ref)]V psi||<=epsilon_geo`,
  `0<=t<=T_L`; if `epsilon_geo<=1/1000` were achieved: `p_out>20449/160000=0.12780625`,
  `p_ret>108241/250000=0.432964`, trace-distance error `<=2/715<1/100`;
  "It is NOT established in this packet." (`MOVING_GEOMETRY.md:171-203`.)
  "A new proof may work; the inherited buffer theorem does not supply it."
  (`REPORT.md:93-94`.)

### 4.7 B2-RB, B2-IS, B2-AS, MC-E1

- **B2-RB** (right boundary): exact `B*B=(g^2/N^2)[(N-m)(c+1)+cG_h]+(1/(16N))[(N-m)n_h^2+G_h n_p^2]`;
  leakage `gamma_w^2=g^2(N-m)/N^2[1+Tr(sigma_w c)]`; "For m<N ... no nonzero
  invariant embedding contained wholly in this prescribed P exists"; the
  absolute Duhamel certificate over `T=640m` is `>=288 sqrt(N-m)`, "a
  limitation of the positive bound, NOT a lower bound on actual error"
  (`RF/grut_b2_right_boundary_20261018/REPORT.md:13-30, 47-54`).
- **B2-IS** (instrument stability): "outcome 3: a specific unresolved
  actual-channel estimate"; success determined by Gram entries `u,v,c`;
  worst success `min(u,v)`, worst retention `min(u_T,v_T)`; crossing operator
  `T_X^2=Pi_X`, `B_c=-jp T_X`; the signed kernel's "actual long-time weight is
  UNBOUNDED here" (`RF/grut_b2_instrument_stability_20261018/REPORT.md:5-9, 13-23, 36-48`).
- **B2-AS** (actual source): `p_out>=[[sqrt(A_L)-sigma]_+^2-delta_L]_+^2`;
  hypotheses `A_L>=17/50`, `|G|<=2/3` would give `>11.458225%` and `>30.75%`;
  "The two source hypotheses themselves remain unproved."
  (`RF/grut_b2_actual_source_20261018/REPORT.md:5-8, 21-29`.)
- **MC-E1** (moving channel encoding): optimal deterministic color-only
  recovery `F_max=(9u+v+6|c|)/(4(3u+v))`; perfect four-dimensional recovery
  requires `u=v=|c|>0`; at `u=v=p` error exactly `(1-|c|/p)/2`
  (`RF/grut_moving_channel_encoding_20261018/REPORT.md:10-20`). Four-event
  diagnostic of the actual B1 law, 1024-mode receiver, radius-8 ball of 4059
  states, coefficients `1/7372800`, `-1/9830400`, `-1/329853488332800`,
  `13/3958241859993600`; `uv-|c|^2>=C0^2 t^18>0` for `0<t<=t0`,
  `C0=1/24513993491703083827200000`,
  `t0=1/102761833874829111436196589800363649819557756928000`
  (`REPORT.md:33-57`). MC-E1-T qutrit protocol has exact conditional identity
  but "is a LOWER-DIMENSION MESSAGE PROTOCOL, not the original arbitrary
  two-qubit task"; "Neither u>=.1 nor u_T>=.25 has been proved on the
  macroscopic windows. There is no operating pass." (`REPORT.md:61-72`.)
  64 checks, `1.167924309 seconds`, `13184KiB`; Challenger "caught an omitted
  electric-field term before freeze" (`REPORT.md:88-94`).

### 4.8 B2 feasibility closeout — the exact numbers

All from `RF/grut_b2_feasibility_closeout_20261018/`; parent pin
`71b5465547c01ffac8ca4666485bb513ce8217bf` (`STOP_DECISION.md:3`;
`MANIFEST.json:15`), which exists in this repository.

| Quantity | Exact value as printed | Source |
|---|---|---|
| Route A quality target | `e^2<=1/25000`; amplitude display `0.00632455532...` | `FEASIBILITY.md:10-12`; `FEASIBILITY.json:43-47` |
| Positive access integral limit | `<=sqrt(128/3125)`, about `0.202386` | `FEASIBILITY.md:14-16`; `FEASIBILITY.json:30-34` |
| Only available occupation bound | `0<=P_X<=1` gives `e<=T/32=20L` | `FEASIBILITY.md:22-23` |
| Crossing certificate, formal minimum `L=2016` | `40,320` (vs target about `0.00632456`) | `FEASIBILITY.md:23-24`; `STOP_DECISION.md:32-35` |
| Crossing certificate, whole word | `2,684,354,880` | `FEASIBILITY.md:24-25`; `STOP_DECISION.md:38-39` |
| Uniform population majorant that would suffice (whole word) | `1/180144028044495360000000` (about `5.55e-24`), "sufficient for this particular method, not a necessary condition" | `STOP_DECISION.md:39-43`; `FEASIBILITY.json:61, 66-69` |
| Same, formal minimum | `1/40642560000000` (about `2.46e-14`) | `FEASIBILITY.json:60, 62-65` |
| Off-diagonal current norm majorant | `||K||<=21/10+5/48+1/1920=1411/640` (geometry `<=g+4epsilon<=21/10`; hopping `<=4(jp+jm)=5/48`; conversion `<=2J=1/1920`); projection commutator needs "no extra factor two" | `FEASIBILITY.md:35-46`; `FEASIBILITY.json:3-6` |
| Absolute-current probability bound | `0<=p_O(t)<=min(1,(1411/640)t)`, "with no positive lower bound" | `FEASIBILITY.md:49-50` |
| Integrated current majorant at `T=640L` | `1411L`; `2,844,576` at `L=2016`; `189,381,236,784` at whole word | `FEASIBILITY.md:52-53`; `FEASIBILITY.json:16-23` |
| Inherited whole-word eligibility | `N>=2^28` (`268435456`), `L>=N/2`, first multiple of 24: `L=134,217,744`, `T=85,899,356,160`; "an eligibility intersection, not a selected favorable record or a new preparation probability" | `FEASIBILITY.md:66-70`; `FEASIBILITY.json:15, 55, 70-73` |
| Formal fixed-path minimum | `L=2016`, `T=1290240` | `FEASIBILITY.md:66`; `FEASIBILITY.json:10-14` |
| Receiver size | `64N^2>=2^62` modes (`4611686018427387904`); dense one-boson complex128 `16*2^62=2^66` bytes (`73786976294838206464`, "64 EiB") vs `2^30` authorized; ratio `2^36` (`68719476736`); actual receiver holds `L-1` bosons | `FEASIBILITY.md:72-75`; `STOP_DECISION.md:51-53`; `FEASIBILITY.json:2, 8-9, 48` |
| Margins | output `114089/4000000=0.02852225`; retention `184281/1000000=0.184281`; success-map budget `514089/404000000`, about `0.0012725` | `FEASIBILITY.md:57-60`; `FEASIBILITY.json:24-27, 50-59` |
| Envelope | 25 s/child, 60 combined child seconds, 1 GiB, 25 MiB new files | `FEASIBILITY.md:84-86` |
| Execution | "only read-only exact rational feasibility arithmetic and deterministic packaging/integrity checks"; `physical_trajectories: 0`; `research_route_selected: null`; `physical_task_verdict: unresolved`; `decision: specialist_review_stop` | `REVIEW_STATUS.md:16-19`; `FEASIBILITY.json:7, 28, 49` |

Status preserved at the stop: "B2-RA: author-derived fixed-path pass,
independent review pending. B2-RB: exact endpoint-sector embedding excluded,
unchanged. B2-IS: exact instrument reduction, not a moving operating pass.
Moving task: neither established nor refuted. No new law, preparation, buffer,
decoder, threshold or interpretation." (`STOP_DECISION.md:59-63`.) Checked
against `RF/grut_b2_real_axis_20261018/REPORT.md:1-19`,
`RF/grut_b2_right_boundary_20261018/REPORT.md:1-9, 29-30`,
`RF/grut_b2_instrument_stability_20261018/REPORT.md:1-9`: consistent.

### 4.9 Verification performed for this ledger entry (read-only, light)

- Ran `python3 -I RF/grut_b2_feasibility_closeout_20261018/feasibility.py`
  from the scratchpad; its output, re-serialized with sorted keys, is
  **byte-identical** to `FEASIBILITY.json`. Independent rational checks:
  `21/10+4(1/64+1/96)+2/3840=1411/640`; `((2^27+23)//24)*24=134217744` (divisible
  by 24); `T0=85899356160`; `T0/32=2684354880`; `(1411/640)T0=189381236784`;
  `(717/2000-1/1000)^2=20449/160000`; `(659/1000-1/1000)^2=108241/250000`;
  `2/715<1/100` true.
- Every `git_blob` in `SOURCE_INDEX.json:6-71` matches `git hash-object` of the
  named live file (11/11 MATCH).
- Commit pins `71b5465, 38fa6a3, d6fc811, 2d72076, 22916ef, b0787a7, 3d26397,
  43bf3af, e6d2b6a, 44f8e82, 44b2748, 830746b, 07609ae` all resolve to commits.
- No trajectory, solver or scan was run; nothing tracked was modified.

## 5. The four notions, kept distinct

| Notion | What the sources establish | What they do not | Sources |
|---|---|---|---|
| **Preparation acceptance** | Accepted record `m>=N/2` at `T_f=4N^2` is conditional on S9's `N>=2^28` threshold; P1 adds `Pr(P>=N/32|accepted)>=1-exp(-289N/65536)-2exp(-N/4)`; the whole-word `L=134,217,744` is an eligibility intersection only | No new preparation probability; recording and loading supplied; S9 not independently verified | `wm2_physical_regime REPORT.md:7-11`; `P1 REPORT.md:24, 51-52, 59-67`; `SPECIALIST_REVIEW.md:60-70`; `FEASIBILITY.md:66-70` |
| **Recoverability** | Fixed-path: exact conditional identity channel once success is paid; moving: `F_max=(9u+v+6|c|)/(4(3u+v))`, perfect iff `u=v=|c|>0`; qutrit (MC-E1-T) exact even with moving geometry; small diagnostic has strict obstruction `uv-|c|^2>=C0^2 t^18` | Actual macroscopic `u,v,c` unknown; "an exchange is not automatically information loss, but unprovided access cannot decode it"; `|gamma|~0.99999583` fixture is not an outgoing transport result | `consolidation REPORT.md:24`; `MC-E1 REPORT.md:14-17, 53-57, 61-72`; `SPECIALIST_REVIEW.md:87-88`; `diagnosis REPORT.md:51-56` |
| **Delivery (output)** | Fixed path `>12.852225%` on `[L/(8J),L/(6J)]`; original 8-event B2 `1.93%`; B1 one-composite arrival `>3/5` at `t=15360` | Moving `u>=1/10` unproved; information access alone "leaves common production, arrival and retention unproved" | `consolidation REPORT.md:22`; `diagnosis REPORT.md:17`; `B1 REPORT.md:23`; `FEASIBILITY.md:30-31`; `MC-E1 REPORT.md:71-72` |
| **Retention** | Fixed path `>43.4281%` on `[2/J,L/(6J)]`; A1 every-time gate fails (`<0.215+2D_A<1/4` at `pi`); B1 matter present `>0.570975` on `[7680,15360]`; margin `0.184281` | Moving `u_T>=1/4` unproved; "do not substitute output or a favorable endpoint"; "Every finite system has recurrence" | `consolidation REPORT.md:23`; `A1 REPORT.md:36-44`; `B1 REPORT.md:22`; `FEASIBILITY.md:57-58`; `STOP_DECISION.md:22`; `B1 FAILURE_MODES.md:36-39` |

Not multiplied across notions: "No listed probabilities are multiplied."
(`diagnosis REPORT.md:20`); "It does not multiply independent
readability/arrival events." (`MOVING_GEOMETRY.md:199`); "One-source 34% is
never a joint 34% claim" (`ARCHIVE_STATUS.md:27-28`).

## 6. The three confirmations requested

**(a) Feasibility stop, not impossibility proof — confirmed.** "Stop is a
method/resource decision; the physical verdict remains unresolved."
(`FEASIBILITY.md:89-90`.) "These large numbers are failed ERROR CERTIFICATES,
not lower bounds on actual error or probabilities of failure."
(`STOP_DECISION.md:47-48`.) "This is not a theorem ruling out a better
analytic/compressed method." (`REVIEW_STATUS.md:14`.) "Dense enumeration is
out of budget, not proof that all analytic/compressed methods fail."
(`SPECIALIST_REVIEW.md:115-116`.) "Moving operation remains unestablished and
unrefuted." (`REVIEW_STATUS.md:26`.) The stop is also a *method* stop: the
available norm majorants are vacuous at the required times independent of
compute — "More compute alone is not prescribed." (`STOP_DECISION.md:55`.)

**(b) Envelope numbers — confirmed exactly** (section 4.8): `1411/640`
majorant; `N=2^28`; `L=134,217,744`; `T=85,899,356,160`; `2^62` modes;
`2^66` bytes vs `2^30`; 25 s / 60 s / 1 GiB / 25 MiB.

**(c) Incompatible state spaces not composable — supported in specific forms,
not as a general principle.** The sentence itself appears only in the handoff
row (`NEXT_CHAT_HANDOFF.md:120`); grep over all B/WM/handoff markdown finds no
primary-source occurrence. What the sources actually say:
(i) the fixed-path reference and the moving law can be compared only on the
SAME embedded space, and the inherited hypothesis fails there (`d_buffer=-1`)
(`MOVING_GEOMETRY.md:43-53, 85-95, 205-210`);
(ii) no invariant endpoint sector exists inside the fixed-word subspace for
`m<N` (`B2-RB REPORT.md:29-30`);
(iii) small diagnostics are not the physical space — "The 4059 states are a
graph ball, not the full physical state space" (`MC-E1 REPORT.md:81-82`);
two-mode receivers are "not the physical 1,024-mode receiver at N=4"
(`P1 REPORT.md:116-117`; `A1 REPORT.md:90-91`);
(iv) one-resource channels do not tensorize to two (`B1 FAILURE_MODES.md:43-45`;
`CHALLENGE_AND_INTEGRATION.md:17`);
(v) A1 and B1 differ in "conversion operator, prepared resource count, energy,
contact strength and degree-dependent diagonal law" so neither rescues the
other (`FAILURE_MODES.md:17-25`).
The sources explicitly leave open that a future same-space comparison could
compose them: "A new proof may work" (`consolidation REPORT.md:93-94`);
"They are not a theorem that no future proof can do so." (`MOVING_GEOMETRY.md:202-203`.)

## 7. Answer

Within the stated supplied-Hamiltonian class, the sector established:
a conditional, externally powered preparation (P1); a size-uniform refutation
of autonomous every-time retention at one permitted point (A1); a conditional
one-composite delivery/retention window at author-derived status (B1); an
author-derived fixed-path two-source pass (`>12.852225%`, `>43.4281%`, exact
conditional channel) for `L in 24N, L>=2016` (B2-RA); exact reasons the
fixed-path pass does not transfer to the moving law (buffer failure,
`J jp=1/245760` witness, lost endpoint identity, horizon `640L`); an exact
instrument reduction to Gram entries `u,v,c,u_T,v_T` (B2-IS); an optimal
recovery formula and a strict small-diagnostic obstruction to perfect
four-dimensional recovery (MC-E1); and a documented stop because every
available certificate (`20L`, `1411L`, dense `2^66` bytes) is vacuous or out
of envelope at the required `N>=2^28` scale. The moving-geometry operating
question is **unresolved: neither established nor refuted**. The original
eight-event B2 endpoint remains **failed** (`1.93%`).

## 8. Closure verdict

**Nothing additional beyond known physics was forced in this class, and the
computations that show it are the sector's own.** Every result is ordinary
finite closed-system quantum mechanics applied to a supplied Hamiltonian with
supplied preparation and readout: "conventional quantum physics ... not an
earned new law of nature" (`A1 REPORT.md:118-120`); "The same leading
predictions follow from a conventional tight-binding channel. Empirical
distinction remains absent." (`B1 CHALLENGE_AND_INTEGRATION.md:19`); "No
distinct empirical prediction follows." (`diagnosis ARCHITECTURAL_VERDICT.md:39`);
"No new law ... gravity interpretation or theory-of-reality claim follows."
(`MC-E1 REPORT.md:102-103`); "There is no new physical result"
(`REVIEW_STATUS.md:23`).

What closure ("bath = rest of the model, never reset") did do here is
structural and was computed, not postulated:
1. With no external sink, retention cannot be bought from the environment.
   A1's matter returns to fuel with total energy conserved
   (`A1 REPORT.md:49-53, 71-72`), and B1's only surviving route is internal
   spatial transport in a finite window with recurrence
   (`B1 REPORT.md:63-66`; `FAILURE_MODES.md:36-39`). This is finite-system
   recurrence, i.e. known physics.
2. Because the complement (receiver, tail, storage) stays in the law, every
   operating estimate must control it; the fixed-path/buffer proof cannot
   (`d_buffer=-1`, `MOVING_GEOMETRY.md:87-88`), and dense computation cannot
   (`2^62` modes, `FEASIBILITY.md:72-75`). Closure therefore raised the
   **proof and resource cost**; it did not supply a new dynamical principle.
3. The information question reduces exactly to Gram entries of the closed
   evolution (`u,v,c`), and the only SWAP primitive is a hop across the
   current bond (`B2-IS REPORT.md:36-44`). Recovery limits are standard
   channel theory (`F_max` formula, `MC-E1 REPORT.md:14`).

No GRUT-specific law, selection rule, prediction or gravity statement was
earned in this sector.

## 9. Classification

`mixed`. Components:
- **Recovered known physics:** Rabi conversion, local perturbation bounds,
  bound-pair bands (`A1 REPORT.md:118-120`); tight-binding transport
  (`B1 CHALLENGE_AND_INTEGRATION.md:19`); finite recurrence
  (`FAILURE_MODES.md:36-39`); optimal recovery fidelity of a two-vector
  Gram instrument (`MC-E1 REPORT.md:14-17`).
- **Conditional constraints (model class stated):** B2-RA bounds for
  `L in 24N, L>=2016` (`consolidation REPORT.md:16-26`); P1 bound for `N`
  multiple of 64, `N>=256` (`P1 REPORT.md:46-52`); B1 window for `>=64`
  events (`B1 REPORT.md:16-26`); A1 refutation at `Omega=13/4, lambda=1`
  (`A1 REPORT.md:30-47`); D1 regime (`wm2_many_carrier REPORT.md:10-12`).
- **Supplied model choices:** conversion law, fuel energy, contact, degree
  counterterm, color symmetry, recording, loading, decoder, thresholds
  (`B1 REPORT.md:55-57`; `SPECIALIST_REVIEW.md:54-56, 70`; `P1 REPORT.md:37-42`).
- **Rejected mechanisms:** autonomous A1 retention (`A1 REPORT.md:3-8`);
  inherited empty-buffer composition for moving B2 (`MOVING_GEOMETRY.md:85-95`);
  gap-only two-level alternative (`FAILURE_MODES.md:3-15`); dense enumeration
  within envelope (`FEASIBILITY.md:72-75`); positive-access and absolute-current
  routes as currently available (`STOP_DECISION.md:17-23`).
- **Conjecture / unresolved:** moving-geometry operation at `N>=2^28`
  (`STOP_DECISION.md:62`); many-carrier theory (`B1 REPORT.md:66`).
- **Unexecuted tasks:** independent review of every author-derived proof;
  Eq. (3) comparison (`MOVING_GEOMETRY.md:178-183`); bounds on `u,v,c,u_T,v_T`
  (`SPECIALIST_REVIEW.md:92-101`); specialist review (packet prepared, "no
  sending action was taken", `REVIEW_STATUS.md:24-25`); six-event formation
  targets in S5/S6/S7/S8.

## 10. Remaining freedom

- A phase-sensitive estimate of the signed kernel `R U_-(t-s) T_X U_+(s) xi_w`
  or a certified compressed/analytic moving-window representation
  (`FEASIBILITY.md:28-29, 76-78`; `SPECIALIST_REVIEW.md:103-109`).
- The same-space comparison Eq. (3) with `epsilon_geo<=1/1000`
  (`MOVING_GEOMETRY.md:171-199`).
- Other parameter points/protocols for autonomous retention — A1 refutes
  uniformity at one point only (`A1 REPORT.md:46-47`); a successor must
  "exhibit a mechanism that defeats this local return estimate"
  (`A1 REPORT.md:126-129`).
- MC-E1-T: a restricted qutrit message protocol with exact conditional
  identity, whose `u>=.1`, `u_T>=.25` remain to be shown (`MC-E1 REPORT.md:61-72`).
- Reopening conditions: "an explicit quantitatively adequate route, a
  substantive proof correction, or independent mathematical assistance"
  (`STOP_DECISION.md:69-70`).
- Not free: a following-head program, new recording time, favorable receiver,
  time averages, new thresholds (`MOVING_GEOMETRY.md:32-35`;
  `SPECIALIST_REVIEW.md:64-65, 117-118`; `A1 REPORT.md:16`).

## 11. Failures, withdrawals, corrections and unexecuted tasks retained

- WM2-S1 closed-seed formation fails (`wm2_defect_assembly REPORT.md:6-8`);
  S2 robust untuned phase rejected (`wm2_fuel_lock REPORT.md:7-9`); S5/S6/S7/S8
  six-event targets unresolved or stopped (`wm2_history_formation REPORT.md:8-9`;
  `wm2_conditioned_receiver REPORT.md:8`; `wm2_positive_residual REPORT.md:3-8`;
  `wm2_schmidt_receiver REPORT.md:3-6`); S9 bounded finite recovery fails
  (`wm2_physical_regime REPORT.md:3`); UC1 universal locality fails
  (`wm2_charge_constraint REPORT.md:1-8`); shared geometry not enforced
  (`wm2_proof_matter_challenge REPORT.md:3-7`).
- P1 first run `OverflowError`, corrected (`FIRST_RUN.txt:1-20`).
- A1 autonomous retention gate refuted; stopped (`A1 REPORT.md:3-8, 126`).
- B1 two-event control `0.182092` below threshold (`FAILURE_MODES.md:29-32`).
- Original eight-event B2 endpoint failed, `1.93%`, unchanged by all later
  work (`diagnosis REPORT.md:12-19`; `ARCHIVE_STATUS.md:13, 33`).
- Moving composition: inherited theorem hypotheses fail at `t=0`
  (`MOVING_GEOMETRY.md:87-88`); B2-RB exact embedding excluded
  (`B2-RB REPORT.md:29-30`); B2-IS and B2-AS unresolved
  (`B2-IS REPORT.md:5-9`; `B2-AS REPORT.md:5-8`).
- MC-E1: strict obstruction to perfect four-dimensional recovery in the small
  diagnostic; Challenger-caught omitted electric-field term before freeze
  (`MC-E1 REPORT.md:53-57, 91-93`).
- Closeout: no route selected, `physical_trajectories: 0`, specialist packet
  not sent (`FEASIBILITY.json:48-49`; `REVIEW_STATUS.md:24-25`).
- Independent mathematical review: absent for every result in this sector
  (`REVIEW_STATUS.md:3-7`; `CHALLENGE.md:1-6`; `ARCHIVE_STATUS.md:34`).

## 12. What must not be promoted

- `>12.852225%` / `>43.4281%` are **fixed-path, L in 24N, L>=2016,
  author-derived** bounds; not moving-geometry results and not reviewed
  (`consolidation REPORT.md:16-31`).
- The MOVING_GEOMETRY margins `0.12780625`, `0.432964`, `2/715` are
  **conditional on an unachieved** `epsilon_geo<=1/1000`
  (`MOVING_GEOMETRY.md:185-193`; `REVIEW_DOSSIER` not read).
- `40,320`, `2,684,354,880`, `2,844,576`, `189,381,236,784`, `5.55e-24` are
  **failed certificates**, not error estimates or failure probabilities
  (`STOP_DECISION.md:45-48`).
- `2^66` bytes is a dense-storage count, "not a memory lower bound for every
  algorithm" (`FEASIBILITY.md:74-75`).
- `|gamma|~0.99999583`, `0.751352`, `0.996551`, `0.879595` are finite-control
  observations, not all-size theorems (`diagnosis REPORT.md:51-56`;
  `B1 REPORT.md:40-42`).
- `70.46%` is a per-fuel marginal, not a joint probability
  (`diagnosis REPORT.md:19`).
- S9's `N>=2^28` and P1's `Pr(accepted AND P>=N/32)>1/4` are conditional on
  the unreviewed S9 theorem (`P1 REPORT.md:59-67`).
- MC-E1-T's exact identity is for a qutrit message, not the original two-qubit
  task (`MC-E1 REPORT.md:69-72`).
- The B1 "3/10 per original seed trial" is conditional on S9 surviving review
  (`B1 REPORT.md:30-32`).
- No empirical prediction, gravity statement or theory-of-reality claim from
  any packet in this sector.

## 13. Handoff-row discrepancies (`NEXT_CHAT_HANDOFF.md:120`)

1. **"Results from incompatible state spaces are not composable"** is a
   paraphrase that occurs in no primary source; the sources make narrower
   claims (section 6c) and explicitly leave composition open to a future
   same-space proof (`MOVING_GEOMETRY.md:202-203`; `consolidation REPORT.md:93-94`).
   Read as a general principle it over-generalizes.
2. **"via the continuation/index documents"** resolves to no path. No file
   named continuation/index exists under `RF/`; `GRUT_RECORD_RECONCILIATION_INDEX_01.md`
   and `NO_GO_LEDGER.md` contain no B2/WM2/B1 entries (grep). The WM/B reports
   were located by directory listing instead.
3. The row frames the stop as "within a specified computational envelope"
   only; the sources frame it as both a missing phase-sensitive **method** and
   a resource limit: "Stop is a method/resource decision" (`FEASIBILITY.md:89`);
   "More compute alone is not prescribed." (`STOP_DECISION.md:55`). Not a
   contradiction, but an incomplete characterization.
4. The row omits that the fixed-path pass itself (B2-RA) is author-derived with
   independent review pending (`STOP_DECISION.md:59`; `consolidation REPORT.md:5-6`),
   and that the original eight-event B2 endpoint remains failed
   (`ARCHIVE_STATUS.md:33`). A reader of the row alone could take "feasibility
   stop" to be the only open item.
5. The row's first sentence (keep the four notions distinct) is consistent
   with the sources (`MC-E1 REPORT.md:99-101`; `SPECIALIST_REVIEW.md:70-73`;
   `diagnosis REPORT.md:20`). No discrepancy.

## 14. Unread or missing

- The "continuation/index documents" named in the handoff row: not found.
- `REVIEW_PACKET.zip`: not opened; its sources are pinned in
  `SOURCE_INDEX.json` and all 11 blob ids match the live files.
- Read only in part: `consolidation/THEOREM.md:1-40`; `B2-RB REPORT.md:1-70`;
  `B2-IS REPORT.md:1-70`; `B2-AS REPORT.md:1-30`; `B2-RA REPORT.md:1-30`;
  `diagnosis ARCHITECTURAL_VERDICT.md:1-40`; `B1 FAILURE_MODES.md:1-45`;
  `B1 CHALLENGE_AND_INTEGRATION.md:1-40`; every `wm2_*/REPORT.md` to line 12
  only.
- Not read: `consolidation/REVIEW_DOSSIER.md`, `EXECUTION_NOTES.md`;
  `MC-E1 RECOVERY_THEOREM.md`, `FINITE_DIAGNOSTIC.md`, `PREDICTIVE_INTERFACE.md`;
  `B1 CENTRAL_RESULT.md`, `GRUT_CANDIDATE.md`, `PREDICTIVE_INTERFACE.md`;
  `P1 MODEL.md`, `HANDOFF.md`; `A1 RETENTION.md`, `NEUTRAL_MATTER.md`,
  `CHALLENGE_AND_HANDOFF.md`; B2 sub-packets SA/SO/SS/SB/BB/kernel note/
  collision channel/scaling analysis/source bridge (known only through
  `ARCHIVE_STATUS.md:9-22`); the original eight-event B2 packet directory
  (`d6fc811`) was not located or read — its numbers come via the diagnosis
  report and archive status.
- No independent (external) review of any result exists in the repository per
  the sources themselves.

## 15. Sources

- `RF/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md` (1-60, 111-124, 172-205); `START_NEXT_CHAT.txt` (1-9)
- `RF/grut_b2_feasibility_closeout_20261018/REVIEW_STATUS.md` (1-26), `FEASIBILITY.md` (1-90), `STOP_DECISION.md` (1-70), `SPECIALIST_REVIEW.md` (1-124), `feasibility.py` (1-50), `FEASIBILITY.json` (1-74), `MANIFEST.json` (1-16), `SOURCE_INDEX.json` (1-73)
- `RF/grut_b2_consolidation_20261018/REPORT.md` (1-95), `MOVING_GEOMETRY.md` (1-214), `ARCHIVE_STATUS.md` (1-41), `CHALLENGE.md` (1-43), `THEOREM.md` (1-40)
- `RF/grut_b2_diagnosis_20261018/REPORT.md` (1-75), `ARCHITECTURAL_VERDICT.md` (1-40)
- `RF/grut_b1_neutral_escape_20261018/REPORT.md` (1-71), `FAILURE_MODES.md` (1-45), `CHALLENGE_AND_INTEGRATION.md` (1-40)
- `RF/grut_moving_channel_encoding_20261018/REPORT.md` (1-108)
- `RF/grut_p1_preparation_20261017/REPORT.md` (1-142), `FIRST_RUN.txt` (1-20)
- `RF/grut_a1_fixed_conversion_20261017/REPORT.md` (1-132)
- `RF/grut_b2_right_boundary_20261018/REPORT.md` (1-70); `RF/grut_b2_real_axis_20261018/REPORT.md` (1-30); `RF/grut_b2_instrument_stability_20261018/REPORT.md` (1-70); `RF/grut_b2_actual_source_20261018/REPORT.md` (1-30)
- `RF/wm2_*/REPORT.md` (1-12 each, thirteen directories)
- Top-level greps: `GRUT_RECORD_RECONCILIATION_INDEX_01.md`, `NO_GO_LEDGER.md`, `AGENT_COORDINATION.md` (no sector entries found)
