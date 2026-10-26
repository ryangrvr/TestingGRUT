# Handoff: kinetic case-L threshold at the heavy-mediator spot point (state on 26 October 2026)

Brief: [BRIEF_KINETIC.md](BRIEF_KINETIC.md). Owner: Claude (issue #13 comment 6413510190). Usage limits stopped the run
after the build stage. **The three build tracks are complete; the solver, verifier and critic did not run.** There is
no f_sym and no kinetic threshold yet. Everything below is same-session AI work on one supplied model; nothing is
observational evidence, a validation or an exclusion.

## Complete

| Track | Deliverable | Main results |
|---|---|---|
| `tracks/kin_source` | `results/source_interface.json` | Codex's DECAY_SOURCE.json reproduced to ≤ 5e-9 (e = 1e-7: pairs/net 7.406624e-4, R = 40.56262, f_rest 0.0015025, ⟨p⟩/m_X 77200.9 at the onset coordinate). Broad source R = 48.50; Codex's band is 0.82 of it. Mediator debris as timed late sources (R = 1.144). All sources finish 220–7000× before onset (a_ref = 7.05e12, T_vis ≈ 44 GeV). **The release's own symmetric X/X̄ population is ~3500× the vector daughters in number**: the vectors supply energy, not pairs. Release X̄ population bracketed by five labelled alternatives (0.13–2.6 X̄ per net charge). |
| `tracks/kin_rates` | `results/rates_interface.json` | One-loop T_c = 18.1 m_X: the broken phase holds through onset (vacuum masses within 10% up to 8–11 m_X; masses(T) tabulated). Energy-dependent kernels for 2→2, s ↔ XX̄ (Bose factors separate), 14 n→2 channels by Berends–Giele recursion, with a resonant-chain (RIS) overlap subtraction. **The old threshold constants contain no double count, but the old s route overestimates ss → XX̄ by 1.6× at 0.3 m_X up to 170× at 12 m_X**, and the old 4→2 constants are valid only for T ≲ 0.03 m_X; above that the rate is the resonant chain through real s and X, far larger. s self-reactions are negligible. Codex's λ₃ and vacuum width confirmed. |
| `tracks/kin_thermo` | `results/thermo_interface.json` | Bose state functions with the explicit scalar reservoir; Codex's control confirmed (scalar fraction 0.332366 at 12 m; T/m 8.1177 → 7.3439 at R = 40.56). Versus our old proxies: same energy and charge give 1.6× the condensed charge and 1.54× the pair inventory; the MB susceptibility overestimates by 18–19% for T ≥ 3 m. Object branch Q_min = 5.2e4; E − Qω = 1148.1 GeV at Q = 2e5; the linear comparator errs by −14% to +12%. Capture kernel with the hard outgoing X (251.6 MeV) and its reverse. Formation rules (a) gas-only, (b) instant condensation, (c) fragmentation: α_max/H = 2.2e13 at onset, but the fastest-growing lump (Q ≈ 2.5e4) is below Q_min, so fragments must merge or dissolve. The two old proxies tabulated at the heavy point for R = 1–90. 23/23 checks. |

## Not done

- `tracks/kin_solve`: partial code only (`ks_common.py`, `thermo_state.py`, `stage2_moments.py`), no results. The
  relaxation stage, moments stage, the e × source × formation runs, the R scan and the controls remain.
- `tracks/kin_verify` and `tracks/kin_critic`: not started.

## What the build results already imply (argued, not computed)

- The old s route made erasure too easy at high T (up to 170×), while the reservoir raises the pair inventory (1.54×).
  Both push toward **less** erasure than the old proxies assumed, i.e. toward the energy-conserving side.
- The formation question is sharper than either proxy: fragments grow 10¹³ times faster than expansion but start
  below the smallest stationary object, so the population depends on merger kinetics, which no rule here derives.
- Since the vectors add energy and not number, the verdict will turn on how the injected heat changes the formation
  of the release's own antiparticles. Expect dependence on the formation rule and on the initial-spectrum bracket.

## To resume

Re-run the workflow script `kinetic-threshold` from run `wf_f8566571-4fd` (completed agents replay from cache); the
solver reads the three interface files above. Then write REVIEW_KINETIC.md and post on issue #13. Codex's 2cbc7f1
(linear modes) awaits review afterwards; its bookkeeping correction (a detectable Higgs portal needs 99.97–99.99%
erasure against the 1% allowance) is accepted.
