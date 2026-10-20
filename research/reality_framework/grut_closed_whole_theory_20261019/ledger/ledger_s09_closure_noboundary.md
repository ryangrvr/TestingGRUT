# Campaign ledger entry: s09_closure_noboundary (adversarially verified and corrected)

Sector: Closure / no-boundary (C1, C2, B three-prescription control, Obadia coefficient, curvature sign).
Repository `/home/user/TestingGRUT`, branch `codex/grut-closure-fork-20261019`. All paths are relative to the repository root. Line numbers are from the files at HEAD `1489773` ("WIP checkpoint 7"). `git diff --stat 96a7e7e HEAD` over `research/reality_framework/grut_closure_fork_20261019`, `research/reality_framework/grut_closure_control_20261019` and `research/reality_framework/grut_closed_whole_handoff_20261019` prints nothing. `git diff --stat 0b6df6f HEAD` over `GRUT_III/01_PRINCIPLES.md` and the v11/v14 calc files also prints nothing. Every cited line is therefore the same at the handoff commit `96a7e7e`. The reader entry gave HEAD as `da8f925`. That commit exists and is an ancestor of HEAD, but HEAD has moved on since. No cited content changed.

Sector commits (all ancestors of HEAD; checked with `git show --stat`):
- `0b6df6f` "Assess closure fork with classical channel gaps and cosmology counterchecks" created the fork packet (15 files).
- `f28b49d` "Clarify observer closure condition and correct B control sources" added `CLOSURE_SELECTION_RESPONSE.md` (177 lines) and amended `MANIFEST.json`, `README.md` and `SOURCES.md`.
- `1997d20` "Construct GRUT-M1 memory model and autonomous quantum writer" is the only commit that touches `grut_closure_control_20261019/`. It added `STATUS.md` (10 lines) and `protocol.json` (47 lines).

Coordination-log entries: `AGENT_COORDINATION.md:5794-5798`, `:5800-5804`, `:5806-5808`.

## 0. Verifier's corrections to the reader entry

I checked the reader entry's citations, numbers and verdict line by line against the sources. Its core verdict holds: C2 selects nothing, the coefficient is prior art, the curvature sign alone is not a numerical falsifier, and B is unexecuted. The corrections and additions are:

1. **Number format.** The answer prints `-7.67e-53` at N=60. No source prints that rounding. `results.json:416` prints `-7.667648073722e-53`. The source prose prints two other forms, and both are quoted below: `-4.54e-5` and `-2.06e-9` (`REPORT.md:107`), and `-4.5402e-5` and `-2.06115e-9` (`DERIVATION.md:150`).
2. **The "executed additions" are rechecks of archived results.**
   - The second-order residual coefficient `s^2/4 - s^3 coth(s/2)/12` is already in archived v11 revision 3: "W_eps / W_dS(H_eff) = 1 + eps^2 [x^2/4 - (x^3/12) coth(x/2)] + O(eps^3)   (reviewers' coefficient confirmed)" (`GRUT_III/calc/v11_quasi_ds_kms_defect.py:15-16`). The packet itself says "This rechecks the conformal part of GRUT III v11" (`DERIVATION.md:101`).
   - The archived v14 docstring already records, qualitatively, the O(eps) non-thermality of an endpoint-anchored (one-sided) detector kernel. The one-sided rate "differs at O(eps) by a frequency-dependent term from Hdot (it is not a pure temperature shift in either sector)" (`GRUT_III/calc/v14_massless_sector_first_order.py:10-12`).
   - What this packet adds is the 60-digit exact power-law verification and the explicit time-domain endpoint coefficient `(s/2)[2 - s coth(s/2)]` (`DERIVATION.md:109-118`).
3. **Omitted corrections from files the entry itself cites.**
   - v11 revision 3 carries a "CORRECTION OF REVISION 2": the earlier headline "c ~ 1.014 and a '1%-of-epsilon frequency dependence'" was "finite-epsilon contamination" (`v11_quasi_ds_kms_defect.py:6-13`).
   - v14 exists because "v13's zero-mode coefficient turned out to depend on the stochastic coarse-graining scale sigma" (`v14_massless_sector_first_order.py:4-5`).
4. **Omitted scope limit (new in this verification).** v11's docstring states "The present universe is NOT slow-roll: eps_0 = (3/2) Omega_m = 0.4725 (Omega_m = 0.315), eta_0 = -3 Omega_L = -2.06". It also says "the present-day kernel is not described by the first-order thermal form at H_eff" and "the adiabatic expansion is not usable at eps_0 ~ 0.5, eta_0 ~ -2" (`v11_quasi_ds_kms_defect.py:27-40`). So `T_eff=H(1-epsilon)/(2 pi)` must not be applied to the present epoch.
5. **Omitted rejected mechanism.** The entry leaves out two statements: "The known no-boundary-to-Euclidean-vacuum connection is a conditional matter state construction, not a theorem that every closed universe is thermal" and "Closure alone fixes neither the observer's restriction nor its thermal spectrum" (`REPORT.md:85-91`).
6. **Omitted caveat on sequestering.** The `k>0` sequestering result also needs separate treatment of "Stiff singular endpoints and their cutoff" (`DERIVATION.md:210-211`).
7. **The coarse-epsilon failure, stated more sharply.** The refined criterion (`probe.py:105`) judges only the last three halving ratios per separation against 4 with tolerance `.1`. It passes only because that window leaves out the earlier ratios:
   - At separation 3, `3.6864287316097104` and `3.8479396999925513` (`results.json:368-369`) lie 0.3136 and 0.1521 from 4. Both are outside `.1`.
   - At separation 2, `4.121409752180324` (`results.json:310`) also lies outside `.1`.
   - Under the original all-ratio bound `3.7<r<4.3` (`probe.py:105`, else-branch), the six-epsilon grid would still fail at separation 3.
   The source's own framing, "tests asymptotic order rather than loosening the initial bound" (`refinement_protocol.json:4`), is kept verbatim alongside these facts.
8. **The check names claim more than the checks test (new).**
   - `endpoint_correction_is_linear_and_nonzero` compares only the smallest-epsilon slope with the analytic coefficient, within 2% (`probe.py:102`, `:107`).
   - `proper_time_translation_is_not_stationary_at_first_order` tests only that `|mean_time_translation_change| > 0.01*epsilon` at epsilon `0.000625` (`probe.py:108`). The translation shifts the mean time by one unit (`probe.py:98`).
   - Neither check tests scaling across the epsilon grid. The linear scaling shows only in the recorded rows. For example, at s=0.5 the translation change is `0.0007739172028123652` at epsilon `0.02` and `2.5877658050247224e-05` at `0.000625` (`results.json:155`, `:190`).
   - Citation fix: the entry cites `probe.py:98, :108` for the endpoint term, but those are the translation lines. The endpoint lines are `probe.py:97`, `:101-102`, `:107`.
9. **Classification fixes.**
   - (i) The entry labels "some explicit closure condition could select a state/cycle" as a conjecture. The sources do not assert it. They pose it as an open question whose "central result is open" (`NEXT_QUESTION.md:4`), with "No such selection rule has been established" (`REPORT.md:157-158`). It is reclassified as an open research question / unexecuted task.
   - (ii) The entry lists Ivo–Li–Maldacena trace matching as "recovered known physics". It was read, not reproduced: "No reproduction of the subregion solution here" (`SOURCES.md:31`). It is an imported literature fact.
   - (iii) Gibbons–Hawking is imported, not recovered (`curvature_transfer_protocol.json:4`).
   - (iv) The archived v11/v14 finite-separation statements are GRUT III's own calculations, and their priority is open (`SOURCES.md:22-25`; `NEXT_QUESTION.md:46`). They are conditional constraints in a stated model class: neither established known physics nor new.
10. **"(a)'s inverse-Gaussian failure" is not a result from this repository.** It is a literature statement about the FLT prescription (`CLOSURE_SELECTION_RESPONSE.md:120`) and an unexecuted pre-registered expectation, `B<0` (`grut_closure_control_20261019/protocol.json:22`).
11. **The B-control description was incomplete.** The control would compute more than widths. It would also compute **saddle selections** from "oriented upward-flow intersections" (`protocol.json:28-29`, `:38`). A pass needs "All three prescription-specific saddle selections and canonical widths reproduced" (`:41`). "If only prescribed saddle widths are obtained, label Step 1 partial" (`:38`). The large-q check compares with "the published asymptotic value with its explicit O(1/q1) correction; do not demand equality at finite size" (`:37`). The single-branch policy is at `:5`.
12. **Provenance qualification (new).** Two statements cannot be checked against git. `registered_before_control_execution: true` (`protocol.json:3`) is the protocol's own declaration, and "The original registration is preserved unchanged" (`STATUS.md:5-6`) is the status file's. `protocol.json` and `STATUS.md` first appear together in `1997d20`, so git holds no earlier version to compare with. No execution artifact exists either, so nothing contradicts the declarations.
13. **Observation sign (new; verifier note, not a source claim).** By the packet's own convention, `Omega_K=-1/(a^2 H^2)` for k=+1 (`REPORT.md:100`), so closed sections give `Omega_K<0`. Both quoted DESI DR2 central values, `.0023` and `.0011`, are positive (`REPORT.md:134-135`). The packet concludes only that "These fits do not establish a positively curved GRUT cosmology" (`REPORT.md:136-137`). This entry draws no further inference.
14. **Environment notes were stale.**
   - HEAD is now `1489773`.
   - `WORKING_THEORY.md` is tracked, not an untracked draft. Its s09 row is at `research/reality_framework/grut_closed_whole_theory_20261019/WORKING_THEORY.md:201`, and it mentions the B control at `:284` and `:333`.
   - `CAMPAIGN_LEDGER.md` in the same packet holds s09 at `:22` and `:1133`, with "Verification: adversarial re-read PENDING (reader entry only)" at `:1246`.
   - An earlier draft of this file put an s09 row at `CAMPAIGN_LEDGER.md:201`. That is wrong: `CAMPAIGN_LEDGER.md:201` is an unrelated FS-1 line.
   - All of these files are downstream of this ledger and are deliberately not used as sources.
15. **Minor citation fixes.**
    - `REPORT.md:177-180` should be `:178-180`, because line 177 is blank.
    - "protocol.json:2" in the question means the fork packet's `research/reality_framework/grut_closure_fork_20261019/protocol.json:2`, not the control's.
16. **Integrity verified (verifier, read-only).** All 14 `MANIFEST.json` file hashes match the current files. The three protocol hashes recorded in `results.json` (`protocol_sha256`, `refinement_protocol_sha256`, `curvature_transfer_protocol_sha256`) also match. I computed these with `python3 -I` and `hashlib`, without running `finalize.py`.

---

## 1. Question asked in this sector

The owner's review of Candidate 1 commit `8fdd674` opened the sector by proposing an A/B/C fork: (A) closure as a completed reading, (B) research on closure, (C) preserve the map (`research/reality_framework/grut_closure_fork_20261019/REPORT.md:1-16`, `:145-159`; `AGENT_COORDINATION.md:5796`). The packet's registered purpose reads: "Resolve the user's proposed A/B/C fork after Candidate 1 commit 8fdd674, without assuming that closure fixes a state, curvature magnitude or thermality" (`research/reality_framework/grut_closure_fork_20261019/protocol.json:2`).

Restated as the parent question: does the closure premise force any of the following beyond adopted GR/QFT and a named state prescription? The premise is "one closed whole", "nothing outside", observer bath = commutant, read either as C1 ("No boundary data of any kind") or as C2 (GRUT III P1 observer restriction).
- (i) A selection rule for the cosmological state or the lapse integration cycle.
- (ii) A new thermality result for a comoving observer in slow roll. This is the proposed D2 claim built on `T_eff=H(1-epsilon)/(2 pi)`.
- (iii) A sign or magnitude of present spatial curvature usable as a falsifier.

(`REPORT.md:153-166`; `NEXT_QUESTION.md:24-31`.)

Follow-up questions (commit `f28b49d`): does the C2 reading supply the missing state-selection condition, and does `0b6df6f` contain a stronger hidden closure condition? Has the B three-prescription control been executed? (`CLOSURE_SELECTION_RESPONSE.md:3-12`; `AGENT_COORDINATION.md:5802`.)

The labels, as the sources use them:
- **C1** = "No boundary data of any kind" (`CLOSURE_SELECTION_RESPONSE.md:157-158`).
- **C2** = observer restriction per P1, with two sub-readings: restriction as the object of prediction, or a ban on any global state (`CLOSURE_SELECTION_RESPONSE.md:7-8`, `:49-55`).
- **D2** = the proposed claim that the slow-roll coefficient is a new GRUT result (`CLOSURE_SELECTION_RESPONSE.md:90`, `:110-111`).

The owner's messages defining C1, C2 and D1–D4 are not in the checkout, and neither is the "D1–D4 closure literature brief" (`SOURCES.md:50-54`).

## 2. Assumptions / model class

- **Thermality check.** A free, massless, conformally coupled scalar in the conformal vacuum of **spatially flat** FLRW: comoving worldline, noncoincident real proper times, fixed mean time, `H(0)=a(0)=1`, and `epsilon=-Hdot/H^2` "with the remaining slow-roll hierarchy appropriately small" (`research/reality_framework/grut_closure_fork_20261019/DERIVATION.md:50-62`; `protocol.json:5`). The exact test uses the constant-epsilon power law `a(t)=(1+epsilon t)^(1/epsilon)`, `eta(t)=-(1+epsilon t)^((epsilon-1)/epsilon)/(1-epsilon)` (`DERIVATION.md:90-92`; `probe.py:55-64`). Archived v11 requires both eps and `eta = eps_dot/(H eps)` to be first-order small (`v11_quasi_ds_kms_defect.py:18`).
- **Curvature checks.** k=+1 FLRW with `Omega_K=-1/(a^2 H^2)` and c=1. Closed global de Sitter `a=cosh(H_Lambda t)/H_Lambda`, continued from a regular Euclidean hemisphere `a_E=sin(H_Lambda tau)/H_Lambda` (`DERIVATION.md:136-152`). A flat three-torus for the topology point (`DERIVATION.md:167-171`).
- **Temperature-transfer check.** It **imports** the Gibbons–Hawking geodesic temperature `H_Lambda/(2 pi)`: "This imported theorem is not rederived by the code" (`curvature_transfer_protocol.json:4`; `DERIVATION.md:177-180`; `probe.py:117`).
- **Sequestering remark.** It assumes the **original global** equation `MP^2 Gmunu=Tmunu-(1/4)gmunu <T>`, homogeneous FLRW, finite integrals and a vanishing endpoint term. Stiff singular endpoints are explicitly left for separate treatment (`DERIVATION.md:199-211`).
- **C2/P1 analysis.** It uses archived GRUT III P1 as written: an IDENTIFICATION with credit 0; the triple `(A(O), ω|_{A(O)}, σ^ω_t)`; bath = commutant `A(O)′` (= `A(O′)` where Haag duality holds); type III₁ algebras. P1 has no tensor factorisation, no reduced density matrix, no influence functional and no dynamics (`GRUT_III/01_PRINCIPLES.md:23-43`). The accounting row reads "Identifications (no content by themselves) | P1" (`:81`).
- **Arithmetic.** Python `3.12.14`, NumPy `2.3.5`, SciPy `1.17.0`, "Python Decimal, 60 digits"; no package installed (`README.md:20-23`; `results.json` key `versions`). `probe.py` was not rerun for this entry, and `results.json` is taken as recorded.
- **Registered B control model class:**
  - Metric `ds^2=-N^2 dt^2/q + q dOmega_3^2; t in [0,1]`, units `8 pi G = 1; H = 1; Lambda = 3; hbar = 1`.
  - Background action `6 pi^2 integral[-qdot^2/(4N)+N(1-H^2 q)]dt`; mode action `one half integral[q^2 zdot^2/N - N ell(ell+2) z^2]dt`.
  - The mode is a "real, massless, minimally coupled spectator with integral Y^2 dOmega_3 = 1; not an inflaton curvature perturbation".
  - Approximation: quadratic mode action on background saddles, leading semiclassical lapse saddle, and "no mode backreaction or fluctuation prefactor correction to the width".
  - Source: `research/reality_framework/grut_closure_control_20261019/protocol.json:6-14`.

## 3. Evidence (with exact numbers)

### 3a. C2 / P1: no state-selection condition (an argument, not a computation)

- **Restriction.** The map `ω_O(A)=ω(A), A∈A(O)` "exists for every state on the encompassing algebra. Consequently the instruction to restrict a state does not choose that state." The vacuum and smooth coherent excitations share the restriction and the commutant assignment but have different one-point functions: "P1 alone cannot fix even the local one-point function. This is an elementary limitation of the definition, not a new quantum-cosmology calculation" (`CLOSURE_SELECTION_RESPONSE.md:28-40`).
- **Modular KMS.** For "a faithful normal state ... [it] is automatic. It does not select the state or set a temperature in physical time." Equating the modular flow with a geometric observer evolution "would add content; P1 has no such general condition for an inflating subregion" (`CLOSURE_SELECTION_RESPONSE.md:42-47`). Archived P1 agrees: "the medium is thermal" is "a tautology", and the geometric cases are IMPORTED (`GRUT_III/01_PRINCIPLES.md:31-39`).
- **Ivo–Li–Maldacena** (JHEP 02 (2025) 124, eqs. (2)–(4)). The trace-saddle field and momentum matching conditions "follow from tracing and stationarity; they are already part of the published prescription. They are not an extra GRUT condition obtained by renaming the complement a bath." The paper leaves normalization and selection open (`CLOSURE_SELECTION_RESPONSE.md:64-71`). The packet read it but did not reproduce it: "No reproduction of the subregion solution here" (`SOURCES.md:31`).
- **Continuum caveat.** There is no intrinsic partial trace and no canonical inside/outside factorisation. A literal subregion density matrix needs a regulator or split construction (Witten 1803.04993 §5.2, §6.4–6.5, §7) (`CLOSURE_SELECTION_RESPONSE.md:73-82`).
- **Recorded verdict.** "C2 identifies the appropriate observer-state object, but P1 does not currently add a state-selection condition. There is no stronger, hidden closure condition in the `0b6df6f` derivations. That commit formulates the selection problem; it does not solve it" (`CLOSURE_SELECTION_RESPONSE.md:7-12`). Any extension "must therefore name a condition on the local state or gluing prescription that is **not** already restriction, trace matching, positivity, normalization, or modular KMS. No such extension has been derived here" (`:84-88`).
- **Second reading of C2** (no global state at all). This "would be a new construction choice: P1 itself uses a global state". It would still need a rule assigning local states and compatibility rules for overlapping algebras, and it "does not by itself choose one of the proposed lapse cycles" (`CLOSURE_SELECTION_RESPONSE.md:49-55`).
- **C1.** "'No boundary data of any kind' does not determine a lapse integration cycle; Dirichlet zero size is itself boundary data. If C1 is **defined** to be (a), its inverse-Gaussian failure is a conditional termination of that definition, not of every meaning of closure. Similarly, normalizability is an existing admissibility condition, not a new GRUT selection principle" (`CLOSURE_SELECTION_RESPONSE.md:157-162`). In the FLT literature the "inverse-Gaussian result belongs to that specified prescription" (`:120`). In this repository it appears only as the unexecuted expectation "B<0" (`grut_closure_control_20261019/protocol.json:22`).

### 3b. Slow-roll "thermality" coefficient (executed recheck, 60-digit Decimal)

- **First-order statement.** The packet rechecks `W(0,s)=W_dS(s;H_eff=1-epsilon)+O(epsilon^2)`, with `H_eff=H+Hdot/H=H(1-epsilon)=a_ddot/a_dot` (`DERIVATION.md:73-88`).
  - Scope: "a local-in-mean-time finite-separation statement"; "An integral over arbitrarily long detector times is not justified by a finite-separation expansion alone" (`DERIVATION.md:101-105`).
  - Archived v11 already states the shift as a THEOREM for "any a(t) with eps AND eta = eps_dot/(H eps) first-order small", with scope "spatially flat FRW" and "fixed-mean-time convention; omega >> eps H" (`GRUT_III/calc/v11_quasi_ds_kms_defect.py:10`, `:18`, `:22-23`, `:38`).
- **Second-order residual coefficient** `s^2/4 - s^3 coth(s/2)/12` (`DERIVATION.md:96-99`; `refinement_protocol.json:5`).
  - Already in archived v11 as "eps^2 [x^2/4 - (x^3/12) coth(x/2)]" (`v11_quasi_ds_kms_defect.py:15-16`).
  - Matched within the `0.02` relative acceptance (`refinement_protocol.json:6`; `probe.py:106`), using the epsilon `0.000625` row only (`probe.py:104`). Analytic vs numerical:

    | s | analytic | numerical |
    |---|---|---|
    | 0.5 | `0.01996887328048337` | `0.01996994944090311` |
    | 1 | `0.06967054885511226` | `0.06968735143690326` |
    | 2 | `0.1246431430004458` | `0.12488851433610203` |
    | 3 | `-0.2357806342106518` | `-0.23469266633507632` |

    Source: `results.json:202-203`, `:260-261`, `:318-319`, `:376-377`.
- **Midpoint-residual halving ratios.** The epsilon grid is `0.02, 0.01, 0.005, 0.0025` plus the refinement values `0.00125, 0.000625` (`protocol.json:6`; `refinement_protocol.json:3`).
  - s=0.5: `4.003454247312896, 4.001725772459261, 4.00086255484067, 4.000431195374955, 4.0002155772764025`.
  - s=3: `3.6864287316097104, 3.8479396999925513, 3.9250792908337377, 3.962808578995216, 3.9814705022438335`.
  - Source: `results.json:193-199`, `:367-373`.
- **Endpoint-anchored kernel `W(t,t-s)`.** Its first-order coefficient `(s/2)[2 - s coth(s/2)]` is "nonzero and depends on separation" (`DERIVATION.md:109-118`). Analytic vs numerical:

  | s | analytic | numerical |
  |---|---|---|
  | 0.5 | `-0.01037352063419957` | `-0.01035051497103085` |
  | 1 | `-0.08197670686932643` | `-0.08186823260824225` |
  | 2 | `-0.6260705709986626` | `-0.6256977223396754` |
  | 3 | `-1.9715612684213035` | `-1.970875055506366` |

  Source: `results.json:200-201`, `:258-259`, `:316-317`, `:374-375`; computed at `probe.py:97`, `:101-102`; 2% acceptance at the smallest epsilon only (`probe.py:107`).
- **Translation of both proper times.** The check `proper_time_translation_is_not_stationary_at_first_order: true` (`MANIFEST.json:104`) tests only `|change| > 0.01*epsilon` at epsilon `0.000625`, for a one-unit mean-time shift (`probe.py:98`, `:108`). The linear scaling shows in the recorded rows: `0.0007739172028123652` at `0.02` and `2.5877658050247224e-05` at `0.000625` for s=0.5 (`results.json:155`, `:190`).
- **The packet's conclusion:** "the full state is nonstationary in comoving proper time; a thermal form at a fixed midpoint does not establish a global KMS state or an exactly thermal, finite-duration detector experiment" (`DERIVATION.md:119-122`). Archived v14 had already said that the one-sided detector rate "is not a pure temperature shift in either sector" at O(eps) (`v14_massless_sector_first_order.py:10-12`). "This calculation does not recompute an entire detector response or v14's minimally coupled IR sector" (`DERIVATION.md:122-124`).
- **v14's minimally coupled sector.** "Its full mode calculation was read, not rerun, in this packet" (`REPORT.md:66-69`). Its stored claim "explicitly excludes its zero-frequency polynomial contributions and specifies the Wigner convention" (`DERIVATION.md:124-125`).
- **Prior art: Obadia,** Phys. Rev. D 78, 083532 (2008), eqs. (20)–(22). The coefficient is `T_coefficient = sqrt(H^2+2 Hdot)/(2 pi) = (H/2 pi)[1-epsilon+O(epsilon^2)]`, and "The coefficient alone is not a thermal-state theorem: the paper separately requires the stationary thermal Wightman structure" (`CLOSURE_SELECTION_RESPONSE.md:92-111`; `DERIVATION.md:127-132`; `SOURCES.md:9`).
  - "This does not establish priority for every finite-separation or minimally coupled statement. It does rule out presenting the coefficient alone as a new GRUT prediction. A publishable novelty claim remains unestablished" (`REPORT.md:73-76`).
  - "D2 is withdrawn as a claimed new result and retained only as an appropriately scoped consistency check. No novelty certificate is issued" (`CLOSURE_SELECTION_RESPONSE.md:110-111`).
  - "No literature search proves that a finite-separation slow-roll lemma has no previous publication" (`SOURCES.md:22-25`).
- **Not applicable at the present epoch.** "The present universe is NOT slow-roll: eps_0 = (3/2) Omega_m = 0.4725 (Omega_m = 0.315), eta_0 = -3 Omega_L = -2.06"; "the present-day comoving vacuum is not described by the first-order thermal form, and the adiabatic expansion is not usable at eps_0 ~ 0.5, eta_0 ~ -2" (`v11_quasi_ds_kms_defect.py:27`, `:39-40`).
- **Closure does not force thermality.** "The known no-boundary-to-Euclidean-vacuum connection is a conditional matter state construction, not a theorem that every closed universe is thermal. ... Closure alone fixes neither the observer's restriction nor its thermal spectrum. The required action, state prescription and observation must still be supplied" (`REPORT.md:85-91`).

### 3c. Curvature sign and magnitude (executed exact identity)

- **Closed de Sitter:** `Omega_K=-1/sinh^2(H_Lambda t)=-1/(e^{2N}-1)` (`DERIVATION.md:146-148`).

  | N | Omega_K | identity residual |
  |---|---|---|
  | 1 | `-0.15651764274966565` | `0.0` |
  | 5 | `-4.5401991009687765e-05` | `9e-60` |
  | 10 | `-2.061153626686912e-09` | `0.0` |
  | 60 | `-7.667648073722e-53` | `1.6e-59` |

  Source: `results.json:383`, `:394`, `:405`, `:416`; computed at `probe.py:112-115`. The prose reads: "At N5 it is already `-4.54e-5`; at N10 it is `-2.06e-9`. The spatial sections remain spheres" (`REPORT.md:107-108`); and "At N5 it is `-4.5402e-5`, and at N10 it is `-2.06115e-9`" (`DERIVATION.md:150`). The check `positive_curvature_can_be_below_1e_minus3: true` is at `results.json:444` (criterion at N=5 only, `probe.py:122`).
- **The packet's reading of this example.** "This exact background is not a fit to our hot-big-bang history; it disproves the proposed implication that a finite near-flatness bound alone excludes all closed no-boundary histories" (`DERIVATION.md:152-155`).
  - Dilution identity: `|Omega_K,f/Omega_K,i|=exp(-2N)(H_i/H_f)^2`.
  - "Without a restriction or probability distribution for N, closure supplies no minimum present curvature magnitude. A measurement consistent with flatness at `1e-3` therefore cannot exclude the entire negative-Omega_K sign class" (`REPORT.md:108-112`; `DERIVATION.md:155-166`).
  - The packet also says: "A model can be statistically disfavored by its curvature distribution even when arbitrarily small negative values remain in its support" (`DERIVATION.md:163-166`).
- **Compactness does not imply positive curvature.** The flat three-torus has `spatial_Ricci_scalar: 0` and volume `a(t)^3 Lx Ly Lz`, and is "not a solution of every no-boundary or sequestering boundary-value problem" (`results.json:428-430`; `DERIVATION.md:167-171`; `REPORT.md:95-97`).
- **Observational context (not computed in the packet).**
  - DESI DR2 table 5 gives `Omega_K=.0023 +/-.0011` (DESI+CMB, curved LambdaCDM) and `.0011 +/-.0013` (DESI+CMB+Pantheon+, evolving dark energy), both 68% intervals.
  - "These fits do not establish a positively curved GRUT cosmology."
  - The CMB-S4 cancellation is recorded as 2025.
  - "No new empirical likelihood or Euclid forecast is calculated in this packet" (`REPORT.md:129-143`).
  - Maldacena 2024: the semiclassical no-boundary weight "heavily favors too little inflation and excessive curvature. That is a problem for a specific probability prescription, not a universal lower bound from compactness" (`REPORT.md:114-118`).
  - Verifier note: by `REPORT.md:100`, k=+1 means `Omega_K<0`, and both DESI central values are positive (`REPORT.md:134-135`). The packet draws no directional conclusion, and none is drawn here.

### 3d. Flat-to-closed temperature transfer (executed; imports Gibbons–Hawking)

- **The comparison.** On the closed slicing `a_ddot/a_dot=H_Lambda coth(H_Lambda t)` and `epsilon_H=-Hdot/H^2=Omega_K`. The ratio to the geodesic temperature is `sqrt(1+|Omega_K|)`, a discrepancy of `|Omega_K|/2+O(Omega_K^2)` (`DERIVATION.md:183-191`).

  | N | relative_temperature_discrepancy | temperature_discrepancy_div_abs_Omega_K |
  |---|---|---|
  | 1 | `0.07541510253002569` | `0.4818313207709409` |
  | 5 | `2.2700737843094575e-05` | `0.49999432487995393` |
  | 10 | `1.0305768128124118e-09` | `0.4999999997423558` |
  | 60 | `3.833824e-53` | `0.49999999519265886` |

  Source: `results.json:389-390`, `:400-401`, `:411-412`, `:422-423`. The check `closed_geometry_changes_flat_temperature_formula_at_curvature_order: true` is at `results.json:445`; its criterion is `|ratio-.5|<1e-4` at N=5 only (`probe.py:123`).
- **The packet's conclusion:** "extending the flat formula to a closed saddle without accounting for curvature already fails in this exactly soluble limit ... This does not contradict v11's explicitly spatially flat result" (`DERIVATION.md:191-195`).

### 3e. Sequestering (analytic, conditional)

- **Original global sequestering.** It implies `<R>=0`. In FLRW, `∫dt a^3 R=6[a^3H]_i^f+6∫dt(ka-a^3H^2)`. "With finite integrals and a vanishing endpoint term, this would require `k integral a dt=integral a^3 H^2 dt>0` for a nonstatic history. That is a conditional dynamical result, with explicit boundary assumptions. It does not follow from the phrase 'nothing outside'. Stiff singular endpoints and their cutoff need separate treatment" (`DERIVATION.md:199-211`).
- **Local four-form version.** It "does not require finite spacetime volume; the original global restriction cannot be assigned to all sequestering theories. Neither derivation supplies a minimum observable |Omega_K|" (`DERIVATION.md:211-215`; `REPORT.md:120-127`).

### 3f. Candidate 1 tree channels (executed here; they belong to the Candidate 1 sector)

- **Six BVP solves:** N = 39944, 39943 and 19972, each at L40/tol 1e-9 and L50/tol 1e-10. All report `success: True` (`results.json` key `classical_BVP_runs`; `probe.py:72-85`).
- **Recorded values** (`results.json:137-140`): `Delta_single_emission_over_gv: -0.19999915211519692`, `Delta_equal_split_over_gv: -3067.397336564565`, `E_N_minus_E_Nminus1_over_gv: 0.8000008478848031`.
- **Refinement changes:** `1.382431946694851e-09`, `1.3751559890806675e-09` and `1.014996087178588e-09` gv (key `classical_energy_changes_over_gv`).
- **Order:** "Tree level, same original couplings. No quantum fission or lifetime claim" (`results.json:143`; `DERIVATION.md:31-48`; `REPORT.md:35-50`).

These solves have no bearing on closure.

### 3g. Integrity

- **Final run.** 13 checks, `all_checks_pass: true`, `elapsed_seconds: 4.2040241949998745`. Scope: "Classical channel arithmetic, a conformal-field finite-separation expansion, and exact geometric counterexamples. No new quantum energy, no general no-boundary theorem, no empirical fit, no novelty certification." (`results.json:447-449`).
- **Finalizer.** `finalize.py` asserts both the final pass and the retained initial failure (`finalize.py:17-18`).
- **Manifest.** `MANIFEST.json` lists 14 files and 13 executed conditions (`MANIFEST.json:4-75`, `:94-108`). Its scope reads "Byte integrity, syntax, links and recorded numerical outcomes. No independent external validation, cosmological path integral, new empirical fit or priority certificate." (`:109`).
- **Verifier's hash check.** I recomputed the hashes read-only. All 14 manifest hashes and the 3 protocol hashes in `results.json` match the current files.

### 3h. B three-prescription control: registered, unexecuted

- **Status.** The protocol declares `registered_before_control_execution: true` with parent `f28b49d` (`research/reality_framework/grut_closure_control_20261019/protocol.json:3-4`). STATUS.md reads: "The protocol was written before any control computation. The owner then steered the active construction toward structured memory and scrambling. No contour or perturbation code was executed for this packet. ... It may be resumed without treating its source expectations as numerical results. This interruption is a change of requested research direction, not a failed three-width result. No B Step 2 or 3 calculation has been performed" (`STATUS.md:1-10`).
- **Files and provenance.** Only `STATUS.md` and `protocol.json` exist in the directory. Both first entered the tree in `1997d20`, the M1 memory commit (`git log -- research/reality_framework/grut_closure_control_20261019`; `AGENT_COORDINATION.md:5808`). Git therefore holds no separately committed earlier registration against which "preserved unchanged" (`STATUS.md:5-6`) could be checked.
- **What it would compute.**
  - **Saddle selections.** For each prescription: "Find all roots of the reduced-action derivative, derive Robin reduction by independent elimination of q0, and calculate oriented upward-flow intersections with representative prescribed contours", using phase tilts and recording "unresolved or regulator-dependent intersections rather than supplying desired saddle lists" (`protocol.json:28-29`).
  - **Widths.** The width coefficient `B = Im[q1^2 zdot(1)/(N z(1))]` in `Psi proportional exp[-(B+i C) z1^2/(2 hbar)]` ("Negative B is an inverse-Gaussian coefficient, not a variance") (`:13`). The modes come from the regular Frobenius branch `z~t^(ell/2)`, "without prescribing Im K or the final width" (`:30`).
  - **Benchmark:** `q1: 40.0, ell: 2, robin_b: 1.5` (`:15`).
  - **Prescriptions** (`:16-20`):
    - (a) "Dirichlet q0=0, positive-real lapse with convergent endpoint deformation; FLT 1703.02076,1705.00192".
    - (b) "Dirichlet q0=0, whole-real lapse contour passing below zero; 1705.05340; report the Gaussian envelope of each WKB component, not a single exact Gaussian for the saddle sum".
    - (c) "Di Tucci-Lehners Robin alpha=-6 pi^2 i, beta=-i b/(3 pi^2), b=1.5, contour emanating from Nstar=3 pi^2 beta toward positive infinity in the corresponding relative cycle".
- **Pre-registered expectations, NOT results** (`:21-26`):
  - (a) "upper-right Dirichlet saddle; B<0, late-size coefficient B -> -ell(ell+1)(ell+2)/H^2".
  - (b) "two lower Dirichlet saddles; positive B envelope, B -> +ell(ell+1)(ell+2)/H^2".
  - (c) "one lower-right original-lapse saddle at admissible b; positive B with the same late-size coefficient".
  - Plus an audit of a possible "2 pi^2 convention mismatch" in 1903.06757.
- **Grids and acceptance** (`:31-38`):
  - `q1_grid [10.0, 40.0, 160.0, 10000.0]`, `ell_grid [2, 3, 5]`, `initial_time_grid [0.00004, 0.00002, 0.00001]`, `ode_tolerances [1e-10, 2e-12]`, `robin_b_checks [1.25, 1.5, 2.0]`, `mode_acceptance_relative_error 0.0000005`.
  - Large-q check: at `q1=10000` compare with "the published asymptotic value with its explicit O(1/q1) correction; do not demand equality at finite size".
  - Contour acceptance: "If only prescribed saddle widths are obtained, label Step 1 partial".
- **Decision rules** (`:40-45`):
  - Pass: "All three prescription-specific saddle selections and canonical widths reproduced in stated leading approximation. This is established physics, not a GRUT prediction or a solution of the full contour dispute".
  - Partial: "Widths agree but a contour or normalization issue remains unresolved".
  - Failure: "A disagreement survives refinement ... do not advance the selection claim".
  - Hypothesis status: "Failure of one closure definition does not disprove every possible observer/state principle. No success criterion promises a universal theory".
- **Limits.** Existing Python/NumPy/SciPy, bounded local runs, no installs (`:46`). Single-branch policy, no Step 2 or 3 (`:5`).
- **Source corrections made before any execution** (`CLOSURE_SELECTION_RESPONSE.md:113-137`):
  - (b) is the whole real axis bypassing zero below: a single inhomogeneous scalar in 1705.05340 §III, eqs. (3.20)–(3.23).
  - 1804.01102 is nonlinear biaxial Bianchi IX, "not directly a third scalar width".
  - (c) is a specified Robin family: "Arbitrary Robin and Neumann data are not synonymous with this family".
  - Yamada 2511.09621 / PRD 113, 063547 (2026) is a distinct regulated prescription that "does not establish 'stable HH only'". Its arXiv v2 "was revised 24 August 2026. Claims about its prescription must specify the version and order of limits".
- **Conventions to fix first** (`CLOSURE_SELECTION_RESPONSE.md:139-155`): "the expected sign must not be put into the numerical boundary conditions as the result being tested", and a de Sitter spectator "should not be silently called an inflaton curvature perturbation".

## 4. Answer

1. **C2 does not supply a state-selection rule.**
   - Restriction exists for every global state. Modular KMS is automatic and sets no physical-clock temperature. The Ivo–Li–Maldacena matching is already part of the published prescription.
   - "C2 as presently stated adds an observer interpretation; it does not yet supply the missing selection condition. There is no alternate derived GRUT condition to substitute for it" (`CLOSURE_SELECTION_RESPONSE.md:169-171`).
   - The no-global-state reading is a new construction that still lacks a local-state rule (`:49-55`).
   - There is no stronger hidden condition in `0b6df6f` (`:7-12`).
   - C1 does not determine a lapse cycle. If C1 is defined as (a), the inverse-Gaussian result from the literature (not reproduced here) terminates only that definition (`:157-162`).
2. **The slow-roll coefficient is not a new thermality result.**
   - `T_eff=H(1-epsilon)/(2 pi)` is the first-order term of Obadia's 2008 coefficient and an archived v11 theorem for flat FLRW.
   - The packet **rechecks** the archived second-order coefficient (`v11:15-16`) and the O(epsilon) non-stationarity of the endpoint kernel (`v14:10-12`). It adds the explicit time-domain endpoint coefficient `(s/2)[2 - s coth(s/2)]` and the 60-digit verification.
   - All of this follows from the adopted conformal vacuum in flat FLRW, not from closure. A thermal form at a fixed midpoint is neither a global KMS state nor a detector result.
   - D2 is withdrawn, no novelty certificate is issued, and priority of the finite-separation lemma is open.
   - The first-order form does not describe the present epoch (`v11:27`, `:39-40`).
   - Sources: `REPORT.md:71-76`, `:85-91`; `DERIVATION.md:50-132`; `CLOSURE_SELECTION_RESPONSE.md:106-111`; `SOURCES.md:22-25`.
3. **The curvature sign alone is not a numerical falsifier or discriminator.**
   - Compactness does not imply `k=+1`: the flat torus has Ricci scalar `0`.
   - `k=+1` supplies no lower bound on `|Omega_K|`. In closed de Sitter it is `-4.54e-5` at N5 and `-2.06e-9` at N10 (`REPORT.md:107`), and `-7.667648073722e-53` at N=60 (`results.json:416`).
   - So a finite flatness bound cannot exclude the sign class. A model's distribution over N can still be statistically disfavoured, but that distribution, along with reheating and conditioning, is a prescription choice (`REPORT.md:93-118`; `DERIVATION.md:134-171`).
   - The flat thermal formula is off on the closed saddle by `sqrt(1+|Omega_K|)` (`DERIVATION.md:173-195`).
   - The only positive-curvature derivation is conditional on original global sequestering with finite integrals, a vanishing endpoint term and separately treated singular endpoints (`DERIVATION.md:199-215`).
4. **The B control is registered and unexecuted.**
   - The model, benchmark (`q1 40.0, ell 2, robin_b 1.5`), three published prescriptions, saddle-selection and width procedure, grids, acceptance and decision rules were fixed in advance.
   - A pass would reproduce published saddle selections and spectator-mode widths. The protocol calls that "established physics, not a GRUT prediction".
   - No contour or perturbation code was run, and Steps 2–3 were not performed. Their content is not defined in any file read (`grut_closure_control_20261019/protocol.json:1-47`; `STATUS.md:1-10`).
5. **Program decision recorded:** option C (preserve the map), with one admissible B question: "does an explicit closure condition select a cosmological quantum state or integration cycle that the supplied dynamics otherwise leave undetermined? No such selection rule has been established in the material inspected here" (`REPORT.md:153-159`).

## 5. Closure verdict

**Nothing additional was forced in the stated class.** The computations and arguments that show it:

- (i) The restriction map exists for every state on the encompassing algebra. So C2/P1 selects no state and fixes no one-point function (`CLOSURE_SELECTION_RESPONSE.md:28-40`).
- (ii) The first-order temperature shift is the Obadia 2008 coefficient, which is also the archived v11 theorem. The executed finite-separation rechecks follow from the adopted conformal vacuum in flat FLRW, not from closure:
  - the archived second-order residual `s^2/4 - s^3 coth(s/2)/12`, matched within 2% at s = 0.5, 1, 2, 3;
  - the endpoint term `(s/2)[2 - s coth(s/2)]`, matched within 2%;
  - the non-stationarity under translation.

  Closure "fixes neither the observer's restriction nor its thermal spectrum" (`DERIVATION.md:50-132`; `results.json:193-203`, `:367-377`; `REPORT.md:85-91`).
- (iii) The exact identity `Omega_K=-1/(e^{2N}-1)` (`-4.5401991009687765e-05` at N=5, `-2.061153626686912e-09` at N=10) and the flat torus (Ricci scalar `0`) show that closure and compactness force neither the curvature sign nor its magnitude (`results.json:383-430`).
- (iv) In the exactly soluble closed case, with Gibbons–Hawking imported, the flat temperature formula is off by `sqrt(1+|Omega_K|)`. The discrepancy over |Omega_K| is `0.49999432487995393` at N=5 (`results.json:389-423`).
- (v) The only positive-curvature derivation is conditional on the original global sequestering equation plus boundary assumptions. That is a supplied global modification, not closure (`DERIVATION.md:199-215`).

This is not a proof of permanent impossibility: "This does not declare the program permanently impossible. It identifies the missing mathematical choice that another construction would have to resolve ... It must actually select the state, however. A name or a compact topology does not make that selection" (`REPORT.md:161-166`; `grut_closure_control_20261019/protocol.json:44`).

## 6. Classification

**mixed**

- **Recovered known physics** (recomputed in the packet):
  - the first-order shift `H_eff=H(1-epsilon)`, which matches Obadia 2008 (`DERIVATION.md:73-88`; `CLOSURE_SELECTION_RESPONSE.md:97-105`);
  - the closed-de-Sitter `Omega_K` identity and the flat-torus geometry (`DERIVATION.md:136-171`).
- **Imported (cited, not recomputed):**
  - Gibbons–Hawking `H_Lambda/(2 pi)` (`curvature_transfer_protocol.json:4`; `DERIVATION.md:177-180`);
  - the restriction and modular-KMS facts (`CLOSURE_SELECTION_RESPONSE.md:28-47`; `01_PRINCIPLES.md:31-39`);
  - Ivo–Li–Maldacena trace matching, read but not reproduced (`CLOSURE_SELECTION_RESPONSE.md:57-71`; `SOURCES.md:31`).
- **Conditional constraint:**
  - archived v11/v14 finite-separation first- and second-order statements, for a conformal (respectively minimally coupled) scalar in flat FLRW at fixed mean time, with priority open and not valid at the present epoch (`DERIVATION.md:101-105`; `v11:15-23`, `:27-40`; `SOURCES.md:22-25`);
  - original global sequestering `<R>=0 ⇒ k∫a dt>0`, under finite integrals, a vanishing endpoint term and treated singular endpoints (`DERIVATION.md:199-211`).
- **Supplied model choice:** action and matter content; the state amplitude with its cycle, boundary conditions, topologies, measure and perturbative domain; the N-distribution, potential, reheating and conditioning (`NEXT_QUESTION.md:6-39`).
- **Rejected mechanism:**
  - the sign-only curvature falsifier, "all sequestering formulations force positive spatial curvature", and D2 novelty (`CLOSURE_SELECTION_RESPONSE.md:166-168`);
  - C1 ≡ prescription (a) as a consequence of closure (`:157-162`);
  - transferring the flat formula to a closed saddle (`DERIVATION.md:191-195`);
  - "every closed universe is thermal" (`REPORT.md:85-86`);
  - compactness ⇒ k=+1 (`REPORT.md:95-97`).
- **Open research question / unexecuted task:**
  - whether any explicit closure condition selects a state or cycle "without using the desired final spectrum or curvature to make the choice". The central result is open and no such rule has been established; the sources do not assert this as a conjecture (`NEXT_QUESTION.md:4`, `:24-31`, `:48-51`; `REPORT.md:155-159`);
  - B Step 1 and Steps 2–3 (`grut_closure_control_20261019/STATUS.md:1-10`; `CLOSURE_SELECTION_RESPONSE.md:172-174`).

## 7. Remaining freedom

- **Action and matter content.** Einstein/inflaton theory or an explicit global modification. A no-boundary wavefunction for one action cannot be combined with another action's constraints without recalculation. "Do not silently identify Candidate 1's binding scalar with an inflaton" (`NEXT_QUESTION.md:6-10`).
- **The state datum:** cycle C, boundary conditions, topologies, measure and perturbative domain. "The phrase 'one closed whole' does not by itself specify this datum" (`NEXT_QUESTION.md:12-22`).
- **Selection without peeking:** whether a closure condition selects among viable cycles or states without using the desired final spectrum or curvature (`NEXT_QUESTION.md:24-31`).
- **Inputs for any present-curvature prediction:** the N-distribution, potential, reheating/history and conditioning (`NEXT_QUESTION.md:33-39`; `DERIVATION.md:162-166`).
- **The thermal bridge:** closed-geometry corrections, what a finite-duration detector actually measures, and the priority of the finite-separation/minimally coupled lemma (`NEXT_QUESTION.md:41-46`; `SOURCES.md:22-25`).
- **Second C2 reading:** a local-state rule and overlap compatibility (`CLOSURE_SELECTION_RESPONSE.md:49-55`).
- **The B control:** the Step 1 outcome (pass, partial or failure), and so whether Steps 2–3 follow (`grut_closure_control_20261019/protocol.json:40-45`).
- **The contour dispute:** FLT vs Diaz Dorronsoro, "not resolved here" (`SOURCES.md:17`).

## 8. Failures, withdrawals, corrections and unexecuted tasks retained

- **Dependency failure.** The initial `probe.py` run failed with `ModuleNotFoundError: No module named 'mpmath'` (the protocol had specified `"mpmath_digits": 60`, `protocol.json:6`). It was amended to 60-digit Decimal and the SymPy dependency was removed, with the rule "do not count a literal polynomial identity as an independent numerical test". No package was installed (`execution_amendment.json:2-8`; `REPORT.md:170-172`).
- **Coarse-epsilon criterion failure.**
  - At separation 3 the coarsest halving ratio, `3.6864`, "failed the proposed 3.7 lower acceptance bound". `results_initial.json` records `midpoint_residual_scales_quadratically: false` and `all_checks_pass: false` (`refinement_protocol.json:2`; `results_initial.json:297-301`, `:341`, `:347`).
  - The declared refinement adds epsilon `0.00125, 0.000625` and judges "the final three halving ratios against 4 with absolute tolerance .1. This tests asymptotic order rather than loosening the initial bound" (`refinement_protocol.json:3-4`; `probe.py:105`).
  - The verifier's arithmetic on the recorded ratios shows the window leaves out `3.6864287316097104` and `3.8479396999925513` at s=3 (`results.json:368-369`) and `4.121409752180324` at s=2 (`results.json:310`). All three lie more than `.1` from 4. Under the original all-ratio bound, the six-epsilon grid would still fail at s=3. The pass is a pass of the changed criterion and is retained as such.
- **Weak checks.** The endpoint and translation checks test only the smallest-epsilon row, and the translation threshold is `0.01*epsilon` (`probe.py:107-108`). Their "linear"/"first order" names rest on the recorded rows, not on the check.
- **Corrections in the archived inputs.**
  - v11 revision 3, "CORRECTION OF REVISION 2": the earlier "c ~ 1.014 and a '1%-of-epsilon frequency dependence'" were "finite-epsilon contamination", and "Richardson extrapolation of the old numbers gives c = 1 to 3e-5" (`GRUT_III/calc/v11_quasi_ds_kms_defect.py:6-13`).
  - v14 exists because "v13's zero-mode coefficient turned out to depend on the stochastic coarse-graining scale sigma" (`v14_massless_sector_first_order.py:4-5`).
  - Both belong primarily to the GRUT III sector. They are retained here because this packet's recheck rests on them.
- **Recorded negative scope result.** At the present epoch (`eps_0 = 0.4725`, `eta_0 = -2.06`), the first-order thermal form does not describe the comoving vacuum kernel, and the adiabatic expansion is not usable (`v11_quasi_ds_kms_defect.py:27-40`).
- **Withdrawn by the owner and recorded:** D2 novelty; a sign-only curvature falsifier; "the claim that all sequestering formulations force positive spatial curvature" (`CLOSURE_SELECTION_RESPONSE.md:166-168`).
- **Corrected B-control source map:** the case (b) contour and model; 1804.01102 is Bianchi IX; Di Tucci–Lehners is a specific Robin family; Yamada is a distinct prescription with a revised arXiv v2 (`CLOSURE_SELECTION_RESPONSE.md:113-137`).
- **Corrected empirical framing:** the Planck-only curvature preference is historical; CMB-S4 was cancelled in 2025; no new likelihood or Euclid forecast was computed (`REPORT.md:129-143`).
- **Unexecuted:**
  - B Step 1 (three widths, contour reproduction), Steps 2–3, the density-matrix spectrum, the KSW bound and any new observable (`CLOSURE_SELECTION_RESPONSE.md:172-174`; `grut_closure_control_20261019/STATUS.md:4-10`);
  - "No no-boundary path integral, cosmological fit, new Q-ball loop spectrum or quantum lifetime was computed. No independent external validation or publication novelty certification has occurred" (`REPORT.md:178-180`);
  - v14's mode calculation was read, not rerun (`REPORT.md:68-69`), and no entire detector response was recomputed (`DERIVATION.md:122-124`);
  - the Ivo–Li–Maldacena subregion solution was not reproduced (`SOURCES.md:31`);
  - the stiff singular sequestering endpoints were not treated (`DERIVATION.md:210-211`).
- **Unavailable:** the Laflamme 1987 full text (`SOURCES.md:12`); the exact CP-04 implementation and the D1–D4 brief, "not found in this checkout" (`SOURCES.md:50-54`).
- **Status of the interruption:** "a change of requested research direction, not a failed three-width result" (`STATUS.md:9-10`). It is neither a pass nor a failure.

## 9. What must not be promoted

- `T_eff=H(1-epsilon)/(2 pi)` as a GRUT prediction or a thermality theorem, or as a description of the present epoch. The second-order coefficient and the endpoint term are not novelties either (`REPORT.md:74-76`; `SOURCES.md:22-25`; `v11:15-16`, `:27-40`).
- "Closed universe ⇒ positive curvature", "⇒ measurable |Omega_K|" or "⇒ thermal". The closed-de-Sitter numbers are a limiting counterexample, "not a fit to our hot-big-bang history" (`DERIVATION.md:152-155`; `REPORT.md:85-86`).
- The DESI DR2 intervals `.0023 +/-.0011` and `.0011 +/-.0013` as establishing a positively curved GRUT cosmology (`REPORT.md:131-137`). By the packet's own sign convention (`REPORT.md:100`), their positive central values do not have the closed sign.
- The B control's `source_expectations_before_execution` (`B -> -/+ ell(ell+1)(ell+2)/H^2`) as results (`STATUS.md:6-7`; `protocol.json:21-26`). Likewise, a future pass must not be read as a GRUT prediction (`protocol.json:41`), and a failure must not be read as disproof of every observer/state principle (`protocol.json:44`).
- P1/C2 as a state-selection law, or modular KMS as a physical-clock temperature (`CLOSURE_SELECTION_RESPONSE.md:7-12`, `:42-47`).
- The conditional sequestering `k>0` result as a consequence of closure, or as a present `|Omega_K|` bound (`DERIVATION.md:207-215`).
- The 13 passing checks (`all_checks_pass: true`) as physical validation: "Verification counts measure executed checks; they are not evidence for a new law" (`REPORT.md:31-33`). The refined asymptotic-order pass is not a pass of the original criterion (`refinement_protocol.json:2-4`).
- The Candidate 1 tree-channel gaps as quantum channel or lifetime results (`REPORT.md:47-50`). Candidate 1's binding scalar is not an inflaton, and a de Sitter spectator is not an inflaton curvature perturbation (`NEXT_QUESTION.md:10`; `CLOSURE_SELECTION_RESPONSE.md:141-142`).
- Option C as a theorem that no responsive-vacuum law is possible: "failure of tested mechanisms is not a theorem excluding every possible responsive-vacuum law" (`REPORT.md:151`).

## 10. Sources read (this verification)

- **Fork packet** (`research/reality_framework/grut_closure_fork_20261019/`):
  - read in full: `REPORT.md` 1–181, `DERIVATION.md` 1–215, `NEXT_QUESTION.md` 1–51, `CLOSURE_SELECTION_RESPONSE.md` 1–177, `SOURCES.md` 1–54, `README.md` 1–35;
  - JSON protocols: `protocol.json` 1–11, `refinement_protocol.json` 1–8, `curvature_transfer_protocol.json` 1–7, `execution_amendment.json` 1–8;
  - scripts: `probe.py` 1–137, `finalize.py` 1–34;
  - `MANIFEST.json` 1–12, 70–110;
  - `results.json` (450 lines): summarized by key with `python3 -I -c json`, and the cited lines located with `grep -n` and printed with `sed`;
  - `results_initial.json` (350 lines): summarized by key, lines 290–350 printed.
- **Control packet:** `research/reality_framework/grut_closure_control_20261019/STATUS.md` 1–10 and `protocol.json` 1–47.
- **Handoff:** `research/reality_framework/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md` lines 24–34, 92–124 and 180–205. The row is at 122. Checked, not used as a primary source.
- **GRUT III:** `GRUT_III/01_PRINCIPLES.md` 15–50 and 76–88. The docstrings only of `GRUT_III/calc/v11_quasi_ds_kms_defect.py` 1–45 and `GRUT_III/calc/v14_massless_sector_first_order.py` 1–45.
- **Coordination log:** `AGENT_COORDINATION.md:5794-5810`.
- **Git:** `git show --stat` for `0b6df6f`, `f28b49d`, `1997d20`; `git log` on both sector directories; `git diff --stat 96a7e7e HEAD` (empty over the three cited directories); `git diff --stat 0b6df6f HEAD` (empty over the P1/v11/v14 files); ancestry checks for `96a7e7e` and `da8f925`.
- **Downstream files, line location only** (not used as sources): `research/reality_framework/grut_closed_whole_theory_20261019/WORKING_THEORY.md:201`, `:284`, `:333` and `CAMPAIGN_LEDGER.md:22`, `:1133`, `:1246`.

## 11. Handoff-row check (`NEXT_CHAT_HANDOFF.md:122`)

The row reads: "C2 as stated adds observer interpretation, not the missing state-selection condition. The slow-roll coefficient is not a new general thermality result; curvature sign without a magnitude is not a discriminator. The registered B three-prescription control remains unexecuted."

- **No contradiction.** All four claims are supported (`CLOSURE_SELECTION_RESPONSE.md:169-171`, `:106-111`; `REPORT.md:71-76`, `:93-112`; `grut_closure_control_20261019/STATUS.md:1-10`).
- **Sourcing discrepancy.** The row links only the C2 response and the B STATUS. Its thermality and curvature claims rest on executed computations in `REPORT.md`, `DERIVATION.md` and `results.json`, which the handoff links only under the Candidate 1 row (`NEXT_CHAT_HANDOFF.md:121-122`).
- **Omissions.** These are not errors, since a one-line row cannot carry everything, but the row leaves out:
  - the three owner withdrawals (`CLOSURE_SELECTION_RESPONSE.md:166-168`);
  - the coarse-epsilon failure and the change of criterion (`refinement_protocol.json:2-4`; `results_initial.json:341`, `:347`);
  - the mpmath failure (`execution_amendment.json:4`);
  - the closed-saddle temperature-transfer counterexample (`DERIVATION.md:173-195`);
  - the conditional sequestering `k>0` result (`DERIVATION.md:199-211`);
  - the C1-vs-(a) distinction (`CLOSURE_SELECTION_RESPONSE.md:157-162`);
  - the rejection of "every closed universe is thermal" (`REPORT.md:85-86`);
  - the archived v11 statement that the first-order form fails at the present epoch (`v11:27-40`).
- **Precision.** "Not a discriminator" means the sign gives no magnitude bound. The source heading is "Curvature sign is not a numerical falsifier" (`REPORT.md:93`), and "A model can be statistically disfavored by its curvature distribution even when arbitrarily small negative values remain in its support" (`DERIVATION.md:163-166`).
- **Scope note.** The handoff's broader list, "Not forced by closure alone: ... compact spatial topology, finite four-volume, positive curvature, a preferred initial state, a thermal bath ..." (`NEXT_CHAT_HANDOFF.md:98-102`), is consistent with this sector. But this sector's computations support only the positive-curvature (given compactness) and thermal-spectrum items (`REPORT.md:85-97`). The handoff itself calls those paragraphs "a construction brief, not a newly executed general theorem" (`NEXT_CHAT_HANDOFF.md:107-108`).
- **Provenance.** The row and `STATUS.md` call the control "registered". Its files entered the tree in `1997d20`, the M1 commit (the handoff itself names `1997d20` as the "original M1 packet", `:32`), not in a dedicated commit (`AGENT_COORDINATION.md:5808`). "Registered before execution" is self-declared and cannot be checked in git, though nothing contradicts it.
- **Environment.** The task text gives the checkout as `96a7e7e`, but HEAD is `1489773`. The cited directories are identical between the two.

## 12. Unread or missing

- The owner's review of `8fdd674` and the messages defining C1, C2 and D1–D4 are not in the repository. They are summarized only at `REPORT.md:1-16`, `CLOSURE_SELECTION_RESPONSE.md:3-12`, `AGENT_COORDINATION.md:5796` and `:5802`.
- The content of B Steps 2–3 is not defined in any file read. They are mentioned only at `CLOSURE_SELECTION_RESPONSE.md:5`, `:174`, `grut_closure_control_20261019/protocol.json:5` and `STATUS.md:10`.
- The "D1–D4 closure literature brief" and the exact CP-04 implementation were "not found in this checkout" (`SOURCES.md:50-54`).
- The Laflamme 1987 full text was inaccessible to the packet authors (`SOURCES.md:12`).
- External literature was not re-read for this entry: Obadia, Ivo–Li–Maldacena, FLT, Diaz Dorronsoro, Di Tucci–Lehners, Yamada, Maldacena 2024, Kaloper–Padilla, DESI DR2, Chen–Lin–Yeom, Witten and Graham. Read scopes are as the packet records them (`SOURCES.md:6-36`).
- `GRUT_III/09_TWO_TEMPERATURES.md` and `13_SYNTHESIS.md` (`SOURCES.md:44`) were not read for this entry. The v11/v14 bodies beyond their docstrings were not read either.
- `research/reality_framework/grut_closed_whole_theory_20261019/WORKING_THEORY.md` and `CAMPAIGN_LEDGER.md` are downstream of this ledger and were deliberately not used.
- `probe.py` was not rerun, and `results.json` is taken as recorded. Only the hashes were recomputed, read-only.
