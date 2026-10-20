# Campaign ledger entry — s07_wm_b1_b2: WM / B1 / B2 matter and moving channels

**Adversarially verified and corrected entry (verifier pass, 2026-10-20).**

Sector key: `s07_wm_b1_b2`. Repository `/home/user/TestingGRUT`, branch
`codex/grut-closure-fork-20261019`. The handoff commit `96a7e7e` is an ancestor of
the checked-out HEAD (`1489773`; the earlier reader's HEAD `da8f925` is also an
ancestor). `git diff --stat 96a7e7e HEAD` over every sector directory used below
(closeout, consolidation, diagnosis, B1, MC-E1, P1, A1, all `grut_b2_*`, all
`wm2_*`) is empty, so every source read here is byte-identical to the handoff
state. All paths are relative to the repository root; `RF` =
`research/reality_framework`. Nothing tracked was modified. The only execution
was the closeout's own read-only rational-arithmetic script, rerun from the
scratchpad with `python3 -I` (section 4.9).

**Verifier's verdict on the earlier entry:** every quoted number matched its
source. The entry was NOT fully correct: (a) it reported the original eight-event
B2 packet as "not located", but it exists at
`RF/grut_b2_collision_channel_20261018/` (commit `d6fc811`) and was read here;
(b) several citations pointed at the wrong line ranges; (c) it omitted several
failures, pre-freeze corrections and unexecuted tasks that the sources record;
(d) two phrases in the closure verdict were the reader's characterization
presented as source content. All are corrected below and listed in section 16.

---

## 1. Question asked in this sector

Closed-whole parent question as asked here: in a **finite closed quantum model**
whose "bath" for the matter subsystem is the rest of the same model — the recorded
geometry word, storage bits, tail registers and the retained `64N^2`-mode
conditional bosonic receiver, never reset — can one **supplied, time-independent
Hamiltonian** (i) start from an *accepted prepared record*, (ii) convert fuel into
matter, (iii) *deliver* it to a readout region, (iv) *retain* it throughout a
stated window, and (v) keep an encoded input *recoverable*; and does any of that
force a law beyond ordinary quantum mechanics?
(Model: `RF/grut_b2_feasibility_closeout_20261018/SPECIALIST_REVIEW.md:15-56`,
receiver "n_R=Q-1" at `:18-19`; "no reset" at
`RF/grut_b2_consolidation_20261018/MOVING_GEOMETRY.md:38-41`; "no unmodeled
energy sink or controller is used" at
`RF/grut_b1_neutral_escape_20261018/REPORT.md:13-14`.)

B2 closeout sub-question: does any **currently available quantitative route**
establish the moving-geometry thresholds — output `>=1/10`, retention `>=1/4`,
conditional successful information trace distance `<=1/100` — within the
authorized compute envelope? (`SPECIALIST_REVIEW.md:85-88`;
`STOP_DECISION.md:6-13`.)

## 2. Handoff row, as read

`RF/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md:120`:

> Keep preparation acceptance, recoverability, delivery and retention distinct.
> The B2 closeout is a feasibility stop within a specified computational envelope,
> not a physical impossibility proof. Results from incompatible state spaces are
> not composable. | [closeout disclosure](../../grut_b2_feasibility_closeout_20261018/REVIEW_STATUS.md),
> original WM/B reports via the continuation/index documents

The handoff is not treated as a primary source; discrepancies are in section 13.

## 3. Assumptions / model class (all supplied, none derived)

1. **Supplied law and constants.** Hard-core plus (states vac,0,1) and minus
   (vac,occupied) species; three-level electric links `-1,0,1` with `+sum E^2`
   and Gauss law `sum_j E_ij - n_+,i + n_-,i = 0`; neutral fuel slots; the
   `64N^2`-mode open-chain receiver with nearest-mode element `-1`; `g in
   [0.9,1.1]`; `delta_i = bar_delta + epsilon_i`, `|bar_delta|<=1/2`,
   `|epsilon_i|<=1/4`; `jp=1/64, jm=1/96, A=13/46080, J=lambda=1/3840,
   Omega=2-2A`. (`SPECIALIST_REVIEW.md:17-56`;
   `RF/grut_b2_consolidation_20261018/THEOREM.md:9-33`.)
2. **Recording, fuel loading and the accepted record are supplied**, "not
   autonomous formation claims" (`SPECIALIST_REVIEW.md:60-70`;
   `MOVING_GEOMETRY.md:41` "Recording and loading are supplied.";
   `RF/grut_p1_preparation_20261017/REPORT.md:3`).
3. **Readout windows, decoder and thresholds are supplied**: output window
   `[L/(8J),L/(6J)]`, retention window `[2/J,L/(6J)]`; output `>=1/10`,
   retention `>=1/4`, conditional error `<=1/100` (`SPECIALIST_REVIEW.md:75-88`).
4. **Every derivation is author-derived, independent review pending**: "All
   prior derivations are author-derived, independent review pending. No external
   reviewer has been contacted." (`SPECIALIST_REVIEW.md:4-5`). The internal
   Challenger "shared context, repository and mathematical ancestry. This is
   internal challenge, not independent mathematical review."
   (`REVIEW_STATUS.md:3-7`; also `RF/grut_b2_consolidation_20261018/CHALLENGE.md:1-6`;
   `RF/grut_b1_neutral_escape_20261018/CHALLENGE_AND_INTEGRATION.md:3-5`).
5. **Closeout computational envelope**: "at most 25 seconds per child, 60
   combined child seconds, 1 GiB memory and 25 MiB new files"
   (`FEASIBILITY.md:84-86`).
6. **Preparation acceptance inherits the unreviewed S9 threshold**: `N>=2^28,
   L=64N^2, 4N^2<=t<=8N^2 => Pr(m>=ceil(N/2))>1/2`
   (`RF/wm2_physical_regime_20261017/REPORT.md:7-8`; "author-derived; independent
   mathematical review is still pending",
   `RF/wm2_proof_matter_challenge_20261017/REPORT.md:11-12`).
7. **Length families**: fixed-path B2-RA family `L in 24N, L>=2016`; inherited
   whole-word moving family `N>=2^28`, `L>=N/2` (`FEASIBILITY.md:66-69`;
   `SPECIALIST_REVIEW.md:60-63`).
8. **P1's persistence is a supplied control**: "At operation start, switch both
   conversion directions off together ... Matter-number persistence relies on
   this supplied control; keeping conversion on is not proved stable."
   (`P1 REPORT.md:31-35`). A1 is the test of leaving conversion on.

## 4. Executed evidence, with sources and exact numbers

### 4.1 WM2 line (one-line verdicts from `head` skims, all 13 packets)

| Packet | One-line verdict (as printed) | Source |
|---|---|---|
| WM2-S1 defect assembly | "Its closed seeded dynamics fail appreciable formation." "The full WM2 success theorem has therefore not been achieved." | `RF/wm2_defect_assembly_20261016/REPORT.md:6-10` |
| WM2-S2 fuel lock | resolves S1 "at a supplied exact resonance"; "rejects a robust untuned-phase claim"; "A separately tested branching proposal reopens a constant-length protected shortcut." — "a constructive conditional result and two architecture failures" | `RF/wm2_fuel_lock_20261016/REPORT.md:3-12` |
| WM2-S3 detuning receiver | "a finite-receiver stability interval, with a spectral limitation"; "does **not** retain S2's isolated real-energy consistency band"; "does not establish a generic spatial phase ... or a distinctive empirical prediction" | `RF/wm2_detuning_receiver_20261016/REPORT.md:1-11` |
| WM2-S4 prepared disorder | "prepared bulk protection, with disordered formation still open"; "does not establish large-size disordered formation" | `RF/wm2_prepared_disorder_20261016/REPORT.md:1-8` |
| WM2-S5 history formation | five-event lower bounds "above 77% and 80%"; "The six-event target remains unresolved under the declared compute budget." | `RF/wm2_history_formation_20261016/REPORT.md:1-9` |
| WM2-S6 conditioned receiver | "formation frontier unchanged"; "Six-event formation at T=144 remains unresolved." | `RF/wm2_conditioned_receiver_20261016/REPORT.md:1-8` |
| WM2-S7 positive residual | "five-event recovery passes; the six-event target remains computationally unresolved"; "the stop is not a negative formation result" | `RF/wm2_positive_residual_20261017/REPORT.md:1-8` |
| WM2-S8 Schmidt receiver | "both original S8 recovery components stop before the target"; "numerical-method outcomes, not evidence against physical formation" | `RF/wm2_schmidt_receiver_20261017/REPORT.md:1-7` |
| WM2-S9 physical regime | "an analytic operational regime is derived; the bounded finite recovery fails"; `N>=2^28, L=64N^2, 4N^2<=t<=8N^2 => Pr(m>=ceil(N/2))>1/2` | `RF/wm2_physical_regime_20261017/REPORT.md:1-11` |
| WM2 proof/matter challenge | "shared geometry is not yet enforced"; S9 self-check "found no failed load-bearing estimate" but "This is the same Builder checking its work." | `RF/wm2_proof_matter_challenge_20261017/REPORT.md:1-12` |
| WM2-UC1 charge constraint | "Primary architecture verdict: negative in the frozen class U."; "universal locality still fails" | `RF/wm2_charge_constraint_20261017/REPORT.md:1-8` |
| WM2-IG1 generator classification | "physical observables do not select the interaction law"; constraints "do **not** force information locality"; relational rule "is an input, not an emergence result" | `RF/wm2_generator_classification_20261017/REPORT.md:1-10` |
| WM2-D1 many carrier | "a conditional many-carrier regime"; "not an unchanged-law result, independent validation, or a completed theory of reality"; `P_N=N/16, N_A=N_B=N/32, 0<=|t|<=N/4096`, error `<=exp(-N/5)+3exp(-N/17)` | `RF/wm2_many_carrier_20261017/REPORT.md:1-12` |

### 4.2 P1 — preparation acceptance with explicit supplied work

- "a controlled preparation regime, not an autonomous populated phase"; "The
  derivations are author-derived; independent proof review is pending."
  (`RF/grut_p1_preparation_20261017/REPORT.md:3-9`).
- Protocol: `T_f=4N^2`; accept `m>=N/2`; `K=N/32` bonds in the first `N/16`
  positions; two `pi/2` pulses; amplitudes independently in `[0.9,1.1]`;
  per-bond pair probability `>=49/64` (`REPORT.md:24-29,46-48`).
- Earned bound, `N` a multiple of 64, `N>=256`:
  `Pr(P>=N/32 | accepted record) >= 1-exp(-289N/65536)-2exp(-N/4)`
  (`REPORT.md:46-53`).
- Conditional on S9: `Pr(accepted AND P>=N/32)>1/4` at `N>=2^28`, "still
  requires the inherited 2^62-mode receiver and formation time 2^58"; "The N=256
  conditional formula lower bound 0.6766 is not a formation guarantee"
  (`REPORT.md:59-67`).
- Operation comparison: `E_use(t)<=e_p+D(t)`, `D(t)<=exp(-N/5)+3exp(-N/12)`,
  `0<=|t|<=N/65536` (`REPORT.md:74-75`).
- Controls: 21,920-state four-event geometry with a **two-mode** diagnostic
  receiver, "Short-preparation error 3.286e-5 <=9.824e-4"; suite "2.459 seconds
  and 117.39 MiB"; "The full-geometry controls use a two-mode receiver, not the
  physical 1,024-mode receiver at N=4" (`REPORT.md:112,116-117,123`).
- **Failure and correction retained**: first run `OverflowError: int too large
  to convert to float` at `N=1024` in `checks.py` `tail`; corrected by log-gamma
  for `k>170` with floor `exp(-700)`; "No model coefficient, time, size, error
  gate, parameter class or resource budget was changed." (`FIRST_RUN.txt:1-20`;
  `REPORT.md:12-13,125-126`).
- No empirical claim: "does not establish an autonomous matter-bearing phase, a
  selected spatial dimension, original CP-12 matter universality, gravity or an
  empirically distinctive GRUT prediction" (`REPORT.md:134-136`).

### 4.3 A1 — autonomous every-time retention gate refuted (one permitted point)

- "**Decision: refuted by a permitted counterexample, at author-derived status
  pending independent proof review.**" (`RF/grut_a1_fixed_conversion_20261017/REPORT.md:3-8`).
- Gate: `Pr(P(t)>=N/64)>=1/4` at **every** `t in [pi/2, pi/2+N/65536]`,
  `N>=2^18`, multiples of 64; "No time average can substitute" (`REPORT.md:14-16`).
- Point: `Omega=13/4`, `lambda=1` (`REPORT.md:29-30`). Result:
  `Pr(P(pi/2)>=K)>0.946-2D_A(pi/2)`, `Pr(P(pi)>=K)<0.215+2D_A(pi)<1/4`,
  `D_A<=exp(-N/5)+exp(-N/12)`; "One point refutes uniformity on the whole class;
  it does not settle other parameter choices or observation protocols."
  (`REPORT.md:39-47`).
- Closed-system accounting: "failure is matter returning to fuel ... No energy
  or physical system disappears when the declared matter count falls."
  (`REPORT.md:49-53`); "Total energy thereafter is conserved." (`REPORT.md:72`).
- Controls: one cell `0.999999995` at `pi/2`, `0.00065247` at `pi`; 773-state
  two-cell `about 1` at `pi/2`, `0.00345813` at `pi`, "energy conserved within
  5.4e-15"; suite `0.9872 s, 121.41 MiB` (`REPORT.md:78,84-85`); "The small
  receiver is not the physical 1,024-mode receiver at N=4" (`REPORT.md:90-91`).
- Surviving module: lowest mobile co-pair band hopping `t_0=1/3840`, dispersion
  error `<8.1e-6`; adjacent one-pair initial weight in that band `at most 1/529`
  (`REPORT.md:96-101`).
- Status of the methods: "Rabi conversion, local perturbation control and
  strong-coupling bound-pair bands—are conventional quantum physics ... not an
  earned new law of nature." (`REPORT.md:118-120`). "Stop A1 under this gate."
  (`REPORT.md:126`).

### 4.4 B1 — a conditional one-composite delivery window

- "supported within the declared conditional one-resource class, at
  author-derived status pending independent proof review"
  (`RF/grut_b1_neutral_escape_20261018/REPORT.md:3-5`).
- Guarantees (`j_+=1/64, j_-=1/96, kappa=1/4, J=1/3840`): matter present
  `>0.570975` for every `t in [7680,15360]`; at `t=15360` intact prefix and
  co-pair at distance `>=4` with probability `>3/5`; exact conditional identity
  channel on arrival; theorem uses `>=64` events (`REPORT.md:16-26`).
- Moving extension: `m>=N/2`, `N>=2^28`, "comparison error is exponentially
  small"; "If S9 also survives review, accepted formation plus arrival exceeds
  3/10 per original seed trial"; "Its 2^62 receiver modes at threshold are not
  reduced by B1." (`REPORT.md:26-33`).
- Controls: 16-event `0.996551` / `0.879595` at `t=15360`; two-event control
  below threshold; primary `0.312144 s`, `71.59 MiB`; total `1.063435 s`
  (`REPORT.md:40-49`). Two-event control `0.182092` at `theta=4` "violates the
  proposed operating threshold if one removes the long-path hypothesis"
  (`FAILURE_MODES.md:29-32`). Two-fuel 773-state: `p(two pairs)=0.751352`,
  `mean pairs=1.735023` at `theta=4`; "not a tensor product of isolated
  one-resource channels" (`FAILURE_MODES.md:43-45`); `0.997750` at `theta=2` vs
  `0.751352` at `theta=4` (`CHALLENGE_AND_INTEGRATION.md:76-77`).
- Recurrence boundary: "Every finite system has recurrence. No permanent locking
  ... The proven statement is one finite use interval." (`FAILURE_MODES.md:36-39`).
- Supplied structure: conversion, near-band fuel energy, weak contact,
  endpoint-degree counterterm, color symmetry (`REPORT.md:55-57`); narrow fuel
  interval `|delta Omega|,|delta lambda|<=J/1000` (`FAILURE_MODES.md:59-60`);
  "The uncorrected architecture has not passed this frozen test."
  (`FAILURE_MODES.md:63-64`).
- Known-physics equivalence: "The same leading predictions follow from a
  conventional tight-binding channel. Empirical distinction remains absent."
  (`CHALLENGE_AND_INTEGRATION.md:19`). "a many-carrier theory is unproved"
  (`REPORT.md:66`).

### 4.5 B2 — original eight-event endpoint failure and its diagnosis

Original packet (located by this verifier: `git show --stat d6fc811` lists only
`RF/grut_b2_collision_channel_20261018/`; the directory is in the live tree):

- "**Decision:** the specified eight-event fixed-record B1 experiment fails its
  outgoing-matter gate at its frozen final time. Its successful-event channel
  does preserve both ordered qubits exactly." (`RF/grut_b2_collision_channel_20261018/REPORT.md:3-6`).
- Table at `theta=Jt=6` (`t=23040`): near encounter at `theta=3`
  `0.66721185 to 0.66721452` (>=1/4 pass); joint outgoing `0.01932660 to
  0.01932751` (>=1/10 **fail**); change from cancelling the co-pair contact term
  `0.00550774 to 0.00550942` (>=1/100 **fail**); conditional two-qubit
  information error "Exactly zero in endpoint-path model" (<=1/100 pass)
  (`REPORT.md:14-22`).
- Two pairs `0.07967519`; each fuel occupation `about 0.70461148`; mean
  remaining pair count `0.59077704`; earlier `0.99774997` two-pair at `theta=2`
  and outgoing `0.32619522` at `theta=4`, "Neither replaces the frozen theta=6
  readout." (`REPORT.md:29-41`). Worker `0.539908 s`, `79.58 MiB`
  (`REPORT.md:78`). "The prepared review request has not been sent; independent
  mathematical review is pending." (`REPORT.md:88-89`).

Diagnosis packet:

- Verdict: "retain escape as a conditional transport mechanism ... B2 still
  fails its original endpoint task." (`RF/grut_b2_diagnosis_20261018/REPORT.md:3-8`).
- Sectors at `theta=6`: no matter, both fuels `48.89%`; one pair `43.14%`; two
  pairs `7.97%`; accepted output `1.93%`; "70.46% was EACH fuel's occupation,
  not the joint both-full probability ... No listed probabilities are
  multiplied." (`REPORT.md:12-21`).
- Refill: `0.690129` back across each boundary bond, fuel `+0.703486`, residual
  `0.013356` (`REPORT.md:25-27`). 38 low modes; high-band weight `1.95e-11`,
  bound `<8.65e-7`; interference `+0.3764` / `-0.09356`; "Pinching and time
  averages are NOT replacement protocols." (`REPORT.md:32-40`).
- Fixture `|gamma|` about `0.99999583`; 773 states; 193-point diagnostic;
  `1.169 s`, `104.11 MiB` (`REPORT.md:51-52,60-68`).
- "No distinct empirical prediction follows." (`ARCHITECTURAL_VERDICT.md:39`);
  "Independent review of B1's formation/use composition ... is still pending."
  (`ARCHITECTURAL_VERDICT.md:36-37`).

### 4.6 B2-RA fixed-path pass and why it does not compose with the moving law

- B2-RA (`RF/grut_b2_real_axis_20261018/REPORT.md:5-19`): for `L in 24N,
  L>=2016`, one-source arrival `A_L>36%`, joint output `>12.852225%` on
  `I_L=[L/(8J),L/(6J)]`, two-pair retention `>43.4281%` on `[2/J,L/(6J)]`,
  exact identity channel; "author-derived, internally challenged, exact
  constants checked; independent proof review pending"; "The old eight-event
  B2 endpoint remains failed." Its run: 110 exact checks, `0.064716 s`,
  `9.25 MiB` (`REPORT.md:52-56`). Pre-freeze: "Two written-bridge gaps were
  corrected before source freeze" (`CHALLENGE.md:39-43`). It relies on "B2-SO's
  existing same-space comparison" (`REPORT.md:44-47`).
- Consolidation (`RF/grut_b2_consolidation_20261018/REPORT.md:5-8,16-31`):
  `p_out>514089/4000000=0.12852225`, `p_ret>434281/1000000=0.434281`
  (`THEOREM.md:95-99`); "The geometry is prepared and FIXED, not dynamically
  retained by this theorem." Consolidation's own run: 39 checks, 67-dim local
  matrix, `0.109865473 s`, `10,496 KiB` (`REPORT.md:72-77`).
- Why it does not compose (`MOVING_GEOMETRY.md`): buffer hypothesis
  `d_buffer=ell-2-b>=0` fails at `t=0` with `d_buffer=-1` (`:76-88`); shifting
  the weight origin gives `sqrt(W(0)) exp(-d_buffer/2)>=sqrt(2e)` (`:92-95`);
  exact witnesses `<target|(H_full^2-H_ref^2)|initial>=J jp=1/245760`, minus
  analogue `J jm=1/368640`, diagonal differs by `+A` (`:112-131`); assisted
  append `-g/N sqrt(n_contact+1)` even with fuel at the head (`:135-142`);
  endpoint no-passing identity lost after extension (`:147-159`); horizon
  `T_L=640L=1,290,240` at `L=2016`, 84x B1's `15360` (`:161-169`).
- Available margin, NOT achieved: if `epsilon_geo<=1/1000` in Eq.(3) then
  `p_out,full>(717/2000-1/1000)^2=20449/160000=0.12780625`,
  `p_ret,full>(659/1000-1/1000)^2=108241/250000=0.432964`, error
  `<=2/715<1/100`; "It is NOT established in this packet." (`:171-199`).
  "They are not a theorem that no future proof can do so." (`:202-203`);
  "A new proof may work; the inherited buffer theorem does not supply it."
  (`REPORT.md:93-94`).
- Archive: "Old finite B2 task: FAILED, unchanged."; "Independent mathematical
  validation: NOT performed" (`ARCHIVE_STATUS.md:32-37`); sub-packets SA/SO/SS/
  SB/AS/kernel/BB/RA listed with dispositions (`ARCHIVE_STATUS.md:9-22`).

### 4.7 B2-RB, B2-IS, B2-AS, MC-E1

- **B2-RB** (`RF/grut_b2_right_boundary_20261018/REPORT.md`): exact
  `B*B=(g^2/N^2)[(N-m)(c+1)+c G_h]+(1/(16N))[(N-m)n_h^2+G_h n_p^2]`;
  `gamma_w^2=g^2(N-m)/N^2 [1+Tr(sigma_w c)]`; "For m<N ... no nonzero invariant
  embedding contained wholly in this prescribed P exists"; at `m=N` a nonzero
  fourth-order coefficient `J jp jm/sqrt(N)` (`:13-34`). Absolute Duhamel
  certificate `>=288 sqrt(N-m)` over `T=640m`, "NOT a lower bound on actual
  error" (`:47-54`). "The parent epsilon<=1/1000 target remains unachieved";
  `717/202000` sufficient ceiling (`:64-68`). 31 controls, 84 local states,
  160 boundary image states; "internal challenge corrected two specifications
  before freezing" (`:72-81`).
- **B2-IS** (`RF/grut_b2_instrument_stability_20261018/REPORT.md`): "specific
  unresolved actual-channel estimate" (`:5-9`); operation determined by Gram
  entries `u,v,c` and `u_T,v_T`: `min(u,v)`, a quadratic, a cubic, `min(u_T,v_T)`
  (`:13-22`); "A common change of accepted environment preserves the instrument
  exactly when its two-vector Gram matrix is preserved. Large departure from the
  old fixed-word sector can therefore be operationally harmless." (`:29-32`);
  the only SWAP primitive is `T_X` with `T_X^2=Pi_X`, `B_c=-jp T_X`; signed
  kernel `R U_-(t-s) T_X U_+(s) xi_w` "actual long-time weight is UNBOUNDED
  here" (`:36-49`). Unachieved `eta_R=eta_T=1/1000` would give output
  `>0.12752225`, retention `>0.433281`, error `<=4000/510089` (`:72-75`). 49
  checks; 571 configurations, 2,538 transitions, 44 exchange hops;
  `0.150232178 s`, `11,648 KiB` (`:79-88`).
- **B2-AS** (`RF/grut_b2_actual_source_20261018/REPORT.md`): "**UNRESOLVED.**"
  (`:5-8`); hypotheses `A_L>=17/50` and `|G|<=2/3` would give `>11.458225%` and
  `>30.75%`; "The two source hypotheses themselves remain unproved." (`:24-29`);
  off-real iteration certifies contraction only at `Im z/J>72` while `.001`
  smoothing at `L=2016` would use `Im z/J<=1/336000` (`:69-73`); 45 checks,
  `0.738991 s`, `10.875 MiB` (`:84-90`); "Prepare focused mathematical
  assistance, then stop ... No reviewer has been contacted" (`:94-98`).
- **MC-E1** (`RF/grut_moving_channel_encoding_20261018/REPORT.md`):
  `F_max=(9u+v+6|c|)/(4(3u+v))`; perfect iff `u=v=|c|>0`; error `(1-|c|/p)/2`
  at `u=v=p` (`:10-18`); four-event diagnostic with the physical 1024-mode
  receiver, radius-8 ball of 4059 states, coefficients `1/7372800`,
  `-1/9830400`, `-1/329853488332800`, `13/3958241859993600` (`:33-46`);
  `uv-|c|^2>=C0^2 t^18`, `C0=1/24513993491703083827200000`,
  `t0=1/102761833874829111436196589800363649819557756928000` (`:48-57`);
  MC-E1-T qutrit protocol "a LOWER-DIMENSION MESSAGE PROTOCOL, not the original
  arbitrary two-qubit task ... Neither u>=.1 nor u_T>=.25 has been proved ...
  There is no operating pass." (`:61-72`); "The 4059 states are a graph ball,
  not the full physical state space." (`:81-82`); 64 checks, `1.167924309 s`,
  `13184 KiB`; Challenger "caught an omitted electric-field term before freeze"
  (`:88-94`; "a missing +3 electric-field term", `INTERNAL_CHALLENGE.md:22-23`);
  "The closed 07609ae comparison branch stays closed." (`:102`).

### 4.8 B2 feasibility closeout — the exact numbers

All from `RF/grut_b2_feasibility_closeout_20261018/`:

- Decision: "neither available route passes quantitative feasibility within the
  existing resources. Stop this moving-geometry branch for specialist review.
  Actual operation remains unresolved, not refuted." (`STOP_DECISION.md:6-8`);
  "no primary A or B execution selected" (`:25-28`). JSON:
  `decision: specialist_review_stop` (`FEASIBILITY.json:7`),
  `physical_task_verdict: unresolved` (`:28`), `physical_trajectories: 0`
  (`:29`), `research_route_selected: null` (`:49`).
- Route A: `e^2<=1/25000` (`FEASIBILITY.md:10-12`; amplitude display
  `0.00632455532...`, `FEASIBILITY.json:43`); `e<=2jp integral sqrt(P_X)`,
  integral `<=sqrt(128/3125)`, about `0.202386` (`FEASIBILITY.md:14-16`);
  only available bound `0<=P_X<=1` gives `e<=T/32=20L`: `40,320` at `L=2016`,
  `2,684,354,880` whole word (`:22-25`); sufficient uniform majorant
  `1/180144028044495360000000` (about `5.55e-24`) whole word
  (`STOP_DECISION.md:39-43`), `1/40642560000000` (about `2.46e-14`) formal
  (`FEASIBILITY.json:60-65`).
- Route B: `||K||<=21/10+5/48+1/1920=1411/640` (`FEASIBILITY.md:35-43`);
  `0<=p_O(t)<=min(1,(1411/640)t)` (`:49-50`); at `T=640L` majorant `1411L`:
  `2,844,576` at `L=2016`, `189,381,236,784` whole word (`:52-53`); "These
  large numbers are failed ERROR CERTIFICATES, not lower bounds on actual error
  or probabilities of failure." (`STOP_DECISION.md:45-48`).
- Margins: output `114089/4000000=0.02852225`; retention
  `184281/1000000=0.184281`; success-map norm budget `514089/404000000`, about
  `0.0012725` (`FEASIBILITY.md:57-60`).
- Size: `N>=2^28=268435456`, `L>=N/2`, first multiple of 24 `L=134,217,744`,
  `T=85,899,356,160` (`FEASIBILITY.md:66-69`; `FEASIBILITY.json:15,55,70-73`);
  receiver `64N^2>=2^62=4611686018427387904` modes (`FEASIBILITY.json:48`);
  dense one-boson complex128 `16*2^62=2^66=73786976294838206464` bytes (64 EiB)
  vs `2^30=1073741824` authorized, ratio `2^36=68719476736`
  (`FEASIBILITY.md:72-75`; `FEASIBILITY.json:2,8-9`; `STOP_DECISION.md:51-53`);
  "it is not a memory lower bound for every algorithm" (`FEASIBILITY.md:75`).
- Formal fixed-path minimum `L=2016`, `T=1290240` (`FEASIBILITY.json:10-14`).
- Nature of the stop: "Stop is a method/resource decision; the physical verdict
  remains unresolved." (`FEASIBILITY.md:89-90`); "More compute alone is not
  prescribed." (`STOP_DECISION.md:55`); "There is no new physical result ... The
  ZIP is prepared for the owner to review or share; no sending action was
  taken." (`REVIEW_STATUS.md:23-25`).

### 4.9 Verification performed for this entry (read-only, light)

- `python3 -I feasibility.py` from the scratchpad: output byte-identical to
  `FEASIBILITY.json` (diff empty). All eight `MANIFEST.json` sha256 values match
  the live files (`MANIFEST.json:3-12`). All 11 `SOURCE_INDEX.json` `git_blob`
  values match `git hash-object` of the live repository paths
  (`SOURCE_INDEX.json:6-71`). Pins `71b5465, 38fa6a3, d6fc811, 2d72076,
  22916ef, b0787a7, 3d26397, 43bf3af, e6d2b6a, 44f8e82, 44b2748, 830746b,
  07609ae, da8f925` all resolve to commits.
- Independent rational checks: `21/10+4(1/64+1/96)+2/3840=1411/640`;
  `((2^27+23)//24)*24=134217744`; `(1411/640)*85899356160=189381236784`;
  `(717/2000-1/1000)^2=20449/160000`; `(659/1000-1/1000)^2=108241/250000`;
  `2/715<1/100`.

## 5. The four notions, kept distinct (never multiplied)

| Notion | Where it is established | Status | Source |
|---|---|---|---|
| Preparation acceptance | S9 threshold `N>=2^28`; P1 `Pr(P>=N/32|accepted)` | S9 unreviewed; P1 author-derived; P1 switches conversion off (supplied control) | `wm2_physical_regime REPORT.md:7-8`; `P1 REPORT.md:31-35,46-67` |
| Delivery (arrival/output) | B1 `>3/5` at `t=15360` (one composite, `>=64` events); B2-RA `>12.852225%` fixed path | author-derived, fixed geometry; moving-geometry delivery unresolved | `B1 REPORT.md:22-26`; `consolidation REPORT.md:20-27`; `STOP_DECISION.md:20,62` |
| Retention | A1 autonomous every-time gate refuted; B1 `>0.570975` on `[7680,15360]`; B2-RA `>43.4281%` on `[2/J,L/(6J)]` | A1 negative at one point; B1/B2-RA fixed-path only; moving unresolved | `A1 REPORT.md:39-47`; `B1 REPORT.md:22`; `consolidation REPORT.md:23`; `STOP_DECISION.md:22` |
| Recoverability | exact identity channel on fixed endpoint path; MC-E1-T exact for a qutrit; arbitrary two-qubit conditional on Gram `u,v,c` | exact only conditionally; strict small-diagnostic obstruction to perfect four-dim recovery | `THEOREM.md:78-88`; `collision_channel REPORT.md:22`; `MC-E1 REPORT.md:14-18,55-57,61-72` |

"No listed probabilities are multiplied." (`diagnosis REPORT.md:21`); "It does
not multiply independent readability/arrival events." (`MOVING_GEOMETRY.md:199`);
"One-source 34% is never a joint 34% claim" (`ARCHIVE_STATUS.md:27-28`).

## 6. The three confirmations requested

**(a) Distinct notions** — confirmed, section 5.

**(b) Feasibility stop, not impossibility** — confirmed. The sources state it
in the strongest terms: "Actual operation remains unresolved, not refuted."
(`STOP_DECISION.md:8`); "These are limitations of this route, not physical
violations." (`FEASIBILITY.md:55`); "This rules out DENSE enumeration within
this envelope; it is not a memory lower bound for every algorithm."
(`FEASIBILITY.md:74-75`); "Positive-access/current norm bounds are vacuous at
the required times, not evidence of actual large error." (`SPECIALIST_REVIEW.md:114-115`);
"Moving operation remains unestablished and unrefuted." (`REVIEW_STATUS.md:26`).
Exact envelope numbers: `1411/640` majorant, `N=2^28`, `L=134,217,744`,
`T=85,899,356,160`, `2^62` modes, `2^66` bytes vs `2^30` (section 4.8).

**(c) "Incompatible state spaces are not composable"** — supported only in the
specific forms the sources state: the same-space comparison hypothesis fails at
`t=0` with `d_buffer=-1` (`MOVING_GEOMETRY.md:85-95`); no invariant endpoint
embedding for `m<N` (`B2-RB REPORT.md:29-30`); small diagnostics are not the
physical state space (`MC-E1 REPORT.md:81-82`; `P1 REPORT.md:116-117`;
`A1 REPORT.md:90-91`); one-resource channels do not tensorize
(`B1 FAILURE_MODES.md:43-45`; `CHALLENGE_AND_INTEGRATION.md:17`). As a general
principle it is NOT what the sources say: B2-IS states that "Large departure
from the old fixed-word sector can therefore be operationally harmless" when the
Gram matrix is preserved (`B2-IS REPORT.md:29-32`), and the consolidation leaves
a future same-space proof open (`REPORT.md:93-94`; `MOVING_GEOMETRY.md:202-203`).

## 7. Answer

Within the supplied-Hamiltonian class this sector established: an externally
powered, conditional preparation with conversion switched off at operation start
(P1); a size-uniform refutation of autonomous every-time retention at one
permitted parameter point (A1); a conditional one-composite delivery/retention
window at author-derived status (B1); an author-derived fixed-path two-source pass
(`>12.852225%` output, `>43.4281%` retention, exact conditional channel) for
`L in 24N, L>=2016` (B2-RA); exact reasons that pass does not transfer to the
moving law (buffer hypothesis fails at `t=0` with `d_buffer=-1`; witnesses
`J jp=1/245760`, `J jm=1/368640`; lost endpoint order identity; horizon `640L`);
an exact instrument reduction to Gram entries `u,v,c,u_T,v_T` (B2-IS); an
optimal color-only recovery formula plus a strict small-diagnostic obstruction to
perfect four-dimensional recovery (MC-E1); and a documented stop because every
available certificate (`20L`, `1411L`, dense `2^66` bytes) is vacuous or out of
envelope at `N>=2^28`. The moving-geometry operating question is unresolved:
neither established nor refuted. The original eight-event B2 endpoint remains
failed at `0.01932660 to 0.01932751` (and its contact-cancel control at
`0.00550774 to 0.00550942`), while its conditional information channel was
exact. The four notions stay distinct and are never multiplied. The closeout is
a feasibility stop within a stated envelope, explicitly not an impossibility
proof. "Incompatible state spaces not composable" holds only in the specific
forms the sources state; the sources leave a future same-space proof open.

## 8. Closure verdict

Nothing additional beyond known physics was forced in this class, and the
sector's own computations show it. Every result is finite closed-system quantum
mechanics applied to a supplied Hamiltonian with supplied preparation and
readout: "conventional quantum physics ... not an earned new law of nature"
(`A1 REPORT.md:118-120`); "The same leading predictions follow from a
conventional tight-binding channel. Empirical distinction remains absent."
(`B1 CHALLENGE_AND_INTEGRATION.md:19`); "No distinct empirical prediction
follows." (`diagnosis ARCHITECTURAL_VERDICT.md:39`); "There is no new physical
result" (`REVIEW_STATUS.md:23`).

What closure (bath = rest of the model, never reset) did do, as recorded in the
sources — the grouping into three points is this ledger's synthesis, not a
source statement:

1. With no external sink, retention cannot be bought from the environment: A1's
   matter returns to fuel, "No energy or physical system disappears"
   (`A1 REPORT.md:51-53`), "Total energy thereafter is conserved."
   (`A1 REPORT.md:72`); B1's surviving route is internal spatial transport in a
   finite window subject to recurrence (`B1 FAILURE_MODES.md:36-39`) — known
   finite-system physics.
2. Because the complement stays in the law, every operating estimate must
   control it: this is why the fixed-path/buffer proof fails to compose
   (`d_buffer=-1`, `MOVING_GEOMETRY.md:87-88`) and why dense computation is out
   of envelope (`2^62` modes, `FEASIBILITY.md:72-75`). Closure raised the proof
   and resource cost; no new dynamical principle was supplied by any packet.
3. The information question reduces exactly to Gram entries of the closed
   evolution with a single SWAP primitive (`B2-IS REPORT.md:36-44`) and to an
   author-derived optimal-fidelity formula (`MC-E1 REPORT.md:10-19`, "external
   review pending", `:3`).

No GRUT-specific law, selection rule, prediction or gravity statement was earned
here (`P1 REPORT.md:134-136`; `B1 REPORT.md:59-61`; `B2-RA REPORT.md:69-72`;
`MC-E1 REPORT.md:102-103`; `STOP_DECISION.md:63`).

## 9. Classification: mixed

- **Recovered known physics**: Rabi conversion, local perturbation bounds,
  bound-pair bands (`A1 REPORT.md:118-120`); tight-binding transport
  (`B1 CHALLENGE_AND_INTEGRATION.md:19`); finite recurrence
  (`B1 FAILURE_MODES.md:36-39`); two-vector Gram recovery fidelity
  (`MC-E1 REPORT.md:10-19`, author-derived).
- **Conditional constraints**: B2-RA bounds for `L in 24N, L>=2016`
  (`consolidation REPORT.md:16-26`); P1 bound for `N` multiple of 64, `N>=256`
  (`P1 REPORT.md:46-53`); B1 window for `>=64` events (`B1 REPORT.md:16-26`);
  D1 regime (`wm2_many_carrier REPORT.md:10-12`); S5 five-event bounds
  (`wm2_history_formation REPORT.md:5-7`).
- **Supplied model choices**: conversion law, fuel energy, contact, degree
  counterterm, color symmetry, recording, loading, decoder, thresholds,
  conversion-off control in P1 (`B1 REPORT.md:55-57`;
  `SPECIALIST_REVIEW.md:54-56,70`; `P1 REPORT.md:31-42`).
- **Rejected mechanisms**: autonomous A1 retention at `Omega=13/4, lambda=1`
  (`A1 REPORT.md:3-8,29-30`); inherited empty-buffer composition for moving B2
  (`MOVING_GEOMETRY.md:85-95`); gap-only two-level alternative
  (`B1 FAILURE_MODES.md:3-15`); dense enumeration within envelope
  (`FEASIBILITY.md:72-75`); positive-access and absolute-current routes as
  currently available (`STOP_DECISION.md:17-23`); S2's robust untuned phase and
  branching shortcut (`wm2_fuel_lock REPORT.md:7-12`); universal locality in
  class U (`wm2_charge_constraint REPORT.md:3-8`).
- **Conjecture / unresolved**: moving-geometry operation at `N>=2^28`
  (`STOP_DECISION.md:62`); many-carrier theory (`B1 REPORT.md:66`); B2-AS source
  hypotheses (`B2-AS REPORT.md:24-29`).
- **Unexecuted**: independent review of every author-derived proof; Eq.(3)
  comparison (`MOVING_GEOMETRY.md:178-183`); bounds on `u,v,c,u_T,v_T`
  (`SPECIALIST_REVIEW.md:92-101`); specialist packet prepared but "no sending
  action was taken" (`REVIEW_STATUS.md:24-25`); earlier B2-AS assistance request
  with "No reviewer has been contacted" (`B2-AS REPORT.md:94-98`); original B2
  "prepared review request has not been sent" (`collision_channel REPORT.md:88-89`);
  six-event formation targets in S5/S6/S7/S8.

## 10. Remaining freedom

1. A phase-sensitive estimate of the signed kernel `R U_-(t-s) T_X U_+(s) xi_w`,
   or a certified compressed/analytic moving-window representation
   (`FEASIBILITY.md:28-29,76-78`; `SPECIALIST_REVIEW.md:103-109`).
2. The same-space comparison Eq.(3) with `epsilon_geo<=1/1000`, which would
   yield `p_out>0.12780625`, `p_ret>0.432964`, error `<=2/715`
   (`MOVING_GEOMETRY.md:171-199`); or, unnormalized-operation errors
   `eta_R=eta_T=1/1000` (`B2-IS REPORT.md:72-75`).
3. Other parameter points or observation protocols for autonomous retention;
   A1 refutes uniformity at one point only; a successor must defeat the local
   return estimate with a specified mechanism (`A1 REPORT.md:46-47,126-129`).
4. MC-E1-T restricted qutrit protocol; its `u>=.1` and `u_T>=.25` remain to be
   shown (`MC-E1 REPORT.md:61-72`).
5. Reopening conditions for the moving branch: "an explicit quantitatively
   adequate route, a substantive proof correction, or independent mathematical
   assistance" (`STOP_DECISION.md:69-70`).
6. Not free: a following-head program, new recording time, favorable receiver,
   time averages, new thresholds (`MOVING_GEOMETRY.md:32-35`;
   `SPECIALIST_REVIEW.md:64-65,117-118`; `A1 REPORT.md:16`;
   `diagnosis REPORT.md:38`).

## 11. Failures, withdrawals, corrections and unexecuted tasks retained

1. WM2-S1 closed-seed formation fails (`wm2_defect_assembly REPORT.md:6-10`);
   S2 robust untuned phase rejected AND branching proposal reopens a protected
   shortcut — "two architecture failures" (`wm2_fuel_lock REPORT.md:7-12`); S3
   does not retain S2's consistency band (`wm2_detuning_receiver REPORT.md:6-7`);
   S4 large-size disordered formation open (`wm2_prepared_disorder REPORT.md:8`);
   S5/S6/S7/S8 six-event targets unresolved or stopped
   (`wm2_history_formation REPORT.md:8-9`; `wm2_conditioned_receiver REPORT.md:8`;
   `wm2_positive_residual REPORT.md:3-8`; `wm2_schmidt_receiver REPORT.md:3-6`);
   S9 bounded finite recovery fails (`wm2_physical_regime REPORT.md:3`); UC1
   universal locality fails (`wm2_charge_constraint REPORT.md:1-8`); IG1
   constraints do not force information locality
   (`wm2_generator_classification REPORT.md:6-7`); shared geometry not enforced
   (`wm2_proof_matter_challenge REPORT.md:3-7`).
2. P1 first run `OverflowError` at `N=1024`; corrected by log-gamma for `k>170`
   with floor `exp(-700)`; no parameter changed (`P1 FIRST_RUN.txt:1-20`).
3. A1 autonomous every-time retention gate refuted; A1 stopped
   (`A1 REPORT.md:3-8,126`).
4. B1 two-event control `0.182092` at `theta=4` below threshold
   (`B1 FAILURE_MODES.md:29-32`); "The uncorrected architecture has not passed
   this frozen test." (`:63-64`); two-pair probability drops `0.997750 ->
   0.751352` from `theta=2` to `4` (`CHALLENGE_AND_INTEGRATION.md:76-77`).
5. Original eight-event B2: joint outgoing `0.01932660 to 0.01932751` fails
   `>=1/10`; contact-cancel change `0.00550774 to 0.00550942` fails `>=1/100`
   (`collision_channel REPORT.md:20-21`); "The prepared review request has not
   been sent" (`:88-89`); unchanged by all later work (`ARCHIVE_STATUS.md:33`;
   `diagnosis REPORT.md:17`; `B2-RA REPORT.md:18`).
6. Moving composition: inherited hypotheses fail before evolution, `d_buffer=-1`
   (`MOVING_GEOMETRY.md:87-88`); B2-RB exact embedding excluded
   (`B2-RB REPORT.md:29-30`); B2-IS and B2-AS unresolved (`B2-IS REPORT.md:5-9`;
   `B2-AS REPORT.md:5-8`); B2-AS off-real iteration estimate fails
   (`Im z/J>72` vs `<=1/336000`, `B2-AS REPORT.md:69-73`).
7. Pre-freeze corrections disclosed: B2-RA "Two written-bridge gaps were
   corrected before source freeze" (`B2-RA CHALLENGE.md:39-43`); B2-RB "internal
   challenge corrected two specifications before freezing"
   (`B2-RB REPORT.md:80-81`); MC-E1 "a missing +3 electric-field term" corrected
   before freeze (`MC-E1 INTERNAL_CHALLENGE.md:22-23`; `REPORT.md:91-93`).
8. MC-E1 strict obstruction to perfect four-dimensional recovery in the small
   diagnostic; no operating pass (`MC-E1 REPORT.md:53-57,71-72`).
9. Closeout: no Route A/B selected, `physical_trajectories: 0`, specialist
   packet prepared but not sent (`FEASIBILITY.json:29,49`;
   `REVIEW_STATUS.md:24-25`).
10. Independent mathematical review absent for every result in this sector
    (`REVIEW_STATUS.md:3-7`; `consolidation CHALLENGE.md:1-6`;
    `ARCHIVE_STATUS.md:34`; `diagnosis ARCHITECTURAL_VERDICT.md:36-37`).

## 12. What must not be promoted

- `>12.852225%` and `>43.4281%` are fixed-path, `L in 24N, L>=2016`,
  author-derived, unreviewed bounds; not moving-geometry results
  (`consolidation REPORT.md:16-31`; `THEOREM.md:147-151`).
- `0.12780625`, `0.432964`, `2/715` are conditional on an unachieved
  `epsilon_geo<=1/1000` (`MOVING_GEOMETRY.md:185-193`); `>0.12752225`,
  `>0.433281`, `<=4000/510089` on unachieved `eta=1/1000`
  (`B2-IS REPORT.md:72-75`); `>11.458225%`, `>30.75%` on unproved source
  hypotheses (`B2-AS REPORT.md:24-27`).
- `40,320`, `2,684,354,880`, `2,844,576`, `189,381,236,784` are failed error
  certificates, not error estimates or failure probabilities
  (`STOP_DECISION.md:45-48`); `5.55e-24` is a sufficient majorant requirement
  for one method, "not a necessary condition" (`STOP_DECISION.md:39-43`).
- `2^66` bytes is a dense-storage count, not a memory lower bound for every
  algorithm (`FEASIBILITY.md:74-75`).
- `|gamma|` about `0.99999583`, `0.751352`, `0.996551`, `0.879595`,
  `0.99774997`, `0.32619522` are finite-control observations at specific times,
  not all-size theorems or replacement readouts (`diagnosis REPORT.md:51-56`;
  `B1 REPORT.md:40-42`; `collision_channel REPORT.md:35-37`).
- `70.46%` is a per-fuel marginal, not a joint both-full probability
  (`diagnosis REPORT.md:19`).
- S9's `N>=2^28` acceptance and P1's `Pr(accepted AND P>=N/32)>1/4` are
  conditional on the unreviewed S9 theorem (`P1 REPORT.md:59-67`).
- MC-E1-T exact identity is for a qutrit message, not the original two-qubit
  task (`MC-E1 REPORT.md:69-72`).
- B1's `>3/10` per original seed trial is conditional on S9 surviving review
  (`B1 REPORT.md:30-32`).
- No empirical prediction, gravity statement or theory-of-reality claim from any
  packet in this sector (`P1 REPORT.md:134-136`; `B1 REPORT.md:59-61`;
  `ARCHITECTURAL_VERDICT.md:39`; `MC-E1 REPORT.md:102-103`;
  `REVIEW_STATUS.md:23`; `B2-RA REPORT.md:69-72`).

## 13. Handoff-row discrepancies (`NEXT_CHAT_HANDOFF.md:120`)

1. **"Results from incompatible state spaces are not composable"** is a
   paraphrase occurring in no primary source. `grep -rn -i "incompatible state
   space|not composable" --include=*.md` finds it only in the handoff row and in
   the derived ledger files of `RF/grut_closed_whole_theory_20261019/`
   (`CAMPAIGN_LEDGER.md:922,960`; `ledger/ledger_s07_wm_b1_b2.md`), which were
   written from the earlier reader's entry after `96a7e7e` and are not primary.
   The sources make narrower claims (section 6c) and B2-IS explicitly says a
   large departure from the fixed-word sector "can therefore be operationally
   harmless" when the Gram matrix is preserved (`B2-IS REPORT.md:29-32`). The
   sources leave composition open to a future same-space proof
   (`consolidation REPORT.md:93-94`; `MOVING_GEOMETRY.md:202-203`). As a general
   principle the row over-generalizes.
2. **"original WM/B reports via the continuation/index documents"** resolves to
   no path: `find` for `*continuation*` returns nothing; the only `*index*`
   files are `GRUT_RECORD_RECONCILIATION_INDEX_01.md` (root) and the closeout's
   `SOURCE_INDEX.json`. `GRUT_RECORD_RECONCILIATION_INDEX_01.md` and
   `NO_GO_LEDGER.md` contain 0 matches for `B2|WM2|B1|B2-RA`. The WM/B reports
   were located by directory listing.
3. The row frames the stop as "within a specified computational envelope" only;
   the sources frame it as a missing phase-sensitive method AND a resource
   limit: "Stop is a method/resource decision" (`FEASIBILITY.md:89`); "More
   compute alone is not prescribed." (`STOP_DECISION.md:55`). Incomplete rather
   than contradictory.
4. The row omits that the fixed-path pass B2-RA is itself author-derived with
   independent review pending (`STOP_DECISION.md:59`; `consolidation
   REPORT.md:5-6`), that the original eight-event B2 endpoint remains failed
   (`ARCHIVE_STATUS.md:33`; `collision_channel REPORT.md:20`), and that A1's
   autonomous retention gate is refuted (`A1 REPORT.md:3`). A reader of the row
   alone could take the feasibility stop to be the only open item.
5. The row's first sentence (keep the four notions distinct) is consistent with
   the sources (`MC-E1 REPORT.md:99-101`; `SPECIALIST_REVIEW.md:70-73`;
   `diagnosis REPORT.md:19-21`). No discrepancy.
6. The row's source link points only at `REVIEW_STATUS.md`, which mentions the
   `1411/640` majorant and `N=2^28` (`:10-11`) but carries none of the other
   envelope numbers; those are in `FEASIBILITY.md`, `FEASIBILITY.json` and
   `STOP_DECISION.md`. Omission, not contradiction.

## 14. Unread or missing

- The "continuation/index documents" named in the handoff row: not found at any
  path (section 13.2).
- `REVIEW_PACKET.zip`: not opened (binary); its 11 sources are pinned in
  `SOURCE_INDEX.json` and every `git_blob` matches the live file.
- Read in full by this verifier: all four closeout `.md` files, `feasibility.py`,
  `FEASIBILITY.json`, `MANIFEST.json`, `SOURCE_INDEX.json`; consolidation
  `REPORT.md`, `MOVING_GEOMETRY.md`, `ARCHIVE_STATUS.md`, `CHALLENGE.md`,
  `THEOREM.md` (151 lines); diagnosis `REPORT.md`, `ARCHITECTURAL_VERDICT.md`;
  B1 `REPORT.md`, `FAILURE_MODES.md` (85), `CHALLENGE_AND_INTEGRATION.md` (80);
  MC-E1 `REPORT.md`; P1 `REPORT.md`, `FIRST_RUN.txt`; A1 `REPORT.md`; B2-RB
  `REPORT.md` (87), B2-RA `REPORT.md` (83), B2-IS `REPORT.md` (96), B2-AS
  `REPORT.md` (105); `collision_channel REPORT.md` (89).
- Read only in part: every `wm2_*/REPORT.md` to line 40; B2-RA `CHALLENGE.md:35-47`;
  MC-E1 `INTERNAL_CHALLENGE.md:8-26`; handoff `NEXT_CHAT_HANDOFF.md:1-60,100-130,160-210`.
- Not read: consolidation `REVIEW_DOSSIER.md`, `EXECUTION_NOTES.md`; MC-E1
  `RECOVERY_THEOREM.md`, `FINITE_DIAGNOSTIC.md`, `PREDICTIVE_INTERFACE.md`; B1
  `CENTRAL_RESULT.md`, `GRUT_CANDIDATE.md`, `PREDICTIVE_INTERFACE.md`; P1
  `MODEL.md`, `HANDOFF.md`; A1 `RETENTION.md`, `NEUTRAL_MATTER.md`,
  `CHALLENGE_AND_HANDOFF.md`; `collision_channel` `ONE_COMPOSITE_CHANNEL.md`,
  `TWO_COMPOSITE_CHANNEL.md`, `REVIEW_REQUEST.md`; B2 sub-packets
  `grut_b2_source_output`, `source_spectrum`, `scaling_analysis`,
  `source_bridge`, `bulk_boundary`, `kernel_note` (known only through
  `ARCHIVE_STATUS.md:9-22`). A grep for `withdraw|retract|corrected|correction`
  across all sector `.md` files found no withdrawal; the corrections found are
  listed in section 11.7.
- No independent (external) mathematical review of any result in this sector
  exists in the repository, per the sources themselves.

## 15. Sources

See the `sources` list in the structured return; every path:lines cited above.
Git: `96a7e7e` and `da8f925` are ancestors of HEAD `1489773`;
`git diff --stat 96a7e7e HEAD` over all sector directories is empty; commit pins
`71b5465, 38fa6a3, d6fc811, 2d72076, 22916ef, b0787a7, 3d26397, 43bf3af,
e6d2b6a, 44f8e82, 44b2748, 830746b, 07609ae` resolve.

## 16. Verifier corrections to the earlier entry

| # | Earlier claim | Problem | Deciding check | Fix |
|---|---|---|---|---|
| 1 | "The original eight-event B2 packet directory (commit d6fc811) was not located or read" | False: `git show --stat d6fc811` lists `RF/grut_b2_collision_channel_20261018/`, present in the live tree | `collision_channel REPORT.md:1-89` | Read; primary numbers added (4.5); contact-cancel fail and unsent review request retained |
| 2 | "39 checks, 67-dim local matrix, 0.109865473 s, 10,496 KiB" cited at `consolidation REPORT.md:5-8,16-31; B2-RA REPORT.md:5-19` | Numbers are at `consolidation REPORT.md:72-77`; B2-RA's own run is 110 checks, `0.064716 s`, `9.25 MiB` | `consolidation REPORT.md:72-77`; `B2-RA REPORT.md:54-56` | Ranges corrected; both runs listed |
| 3 | `physical_trajectories: 0 ... (FEASIBILITY.json:7,28,48-49)` | `physical_trajectories` is at `:29`; `:48` is `receiver_modes_at_minimum_N` | `FEASIBILITY.json:28-29,48-49` | Cited `:29,49` |
| 4 | "Every derivation is author-derived ... (REVIEW_STATUS.md:3-7)" | That text is `SPECIALIST_REVIEW.md:4-5`; `REVIEW_STATUS.md:3-7` says internal challenge is not independent review | `SPECIALIST_REVIEW.md:4-5` | Both cited |
| 5 | "probabilities are not multiplied (diagnosis REPORT.md:20)" | Sentence is at `:21` | `diagnosis REPORT.md:21` | Cited `:19-21` |
| 6 | "matter returns to fuel with total energy conserved (A1 REPORT.md:49-53)" | "Total energy thereafter is conserved" is at `:72` | `A1 REPORT.md:72` | Added `:72` |
| 7 | "standard channel-fidelity formulas (MC-E1 REPORT.md:14)" | Source calls it a proved optimal fidelity, author-derived, review pending; "standard" is the reader's word | `MC-E1 REPORT.md:3,10-19` | Reworded |
| 8 | "5.55e-24 ... failed error certificates" | `5.55e-24` is a sufficient population-majorant requirement, "not a necessary condition" | `STOP_DECISION.md:39-43` | Separated |
| 9 | Closure-verdict points (1)-(3) presented as source content | They are the ledger's synthesis | — | Labelled as synthesis |
| 10 | Omitted: contact-cancel control fail; B2-RA/B2-RB/MC-E1 pre-freeze corrections; B2-AS unsent assistance request and `Im z/J` failure; S2 second architecture failure; B1 uncorrected-architecture failure and `0.997750->0.751352`; P1 conversion-off supplied control; B2-RB `717/202000`, `J jp jm/sqrt(N)`; B2-IS `eta=1/1000` margins and run numbers; S3/S4/IG1 verdicts absent from the JSON evidence | Failures/corrections present in sources but omitted | Sections 4.1, 4.4-4.8, 11 | Added |
| 11 | "HEAD da8f925" | HEAD is now `1489773`; both are ancestors; sector diff still empty | `git merge-base --is-ancestor` | Updated |
| 12 | Partial-read list ("THEOREM.md:1-40", "B2-RB 1-70", etc.) | These files have 151/87/96/105/83 lines and were read in full here | `wc -l` | Updated in 14 |
