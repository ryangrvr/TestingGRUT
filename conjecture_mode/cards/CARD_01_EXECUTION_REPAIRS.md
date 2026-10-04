# CARD #1 v1 — EXECUTION REPAIR LOG

**Scope.** These are execution and infrastructure repairs only, authorized by the owner convergence ruling (`CARD_01_OWNER_CONVERGENCE_RULING.md` §1). None of them changes:
- the likelihoods, the model, the thresholds or the result logic;
- any **minimizer setting** (BOBYQA, `rhoend` default 0.05, 3 seeded starts, 0.2 spread flag — all unchanged).

The owner threshold ruling and the run configuration are not amended. Card #1 v1 remains under its frozen run configuration.

| ID | Event | Repair |
|---|---|---|
| **ER-A** | **OOM.** The first driver launch used 4 parallel single-thread workers. Each fit needs about 3.7 GB RSS (ACT DR6 lensing data plus CamSpec), and 4 concurrent fits exceeded the container's memory cgroup. The kernel log shows `Memory cgroup out of memory: Killed process …`. **58 jobs ended with rc = −9** before finishing. Only jobs started after other workers had died completed (ε = −1 × 3 starts). | Driver changed to `NW = 2` workers with `OMP_NUM_THREADS = cpu_count // NW` (= 2 threads each). This is a thread count, not a minimizer setting. The original log is kept as `d3_driver_attempt1_oom.log` (scratch). |
| **ER-B** | **Phase ran too early.** Because the killed P1 jobs returned, the driver went on to P2. It started the free-ε fits seeded from the best **available** grid point (ε = 0), not from the complete grid as the spec requires. Two of the three finished. | The two fits are **archived and excluded** in `d3_out_archive_misseeded_free/` (scratch), not deleted. A **driver guard** now prevents P2 until all 20 × 3 primary grid fits exist. The free-ε fits will be rerun under the frozen protocol with the correct best-grid seed. |
| **ER-C** | **Container restart.** The container was recycled while the session was idle (uptime reset to 0 min at 14:44 UTC), which killed the relaunched driver and its two in-flight fits (CPL s0, s1; no output written). Completed result files, the environment and the downloaded products all survived. | The driver is relaunched and resumable (completed fits are skipped). It is kept under harness-tracked supervision so the session stays active. |
| **ER-F** | **Second container restart** (≈ 23:12 UTC). Sixty-five of 66 primary fits were complete. The in-flight fit, free-ε start 2, was killed without writing output. | The driver was relaunched (resumable). The fit was rerun under the unchanged frozen protocol from the same best-grid seed (ε = −0.080) as starts 0 and 1. The driver logs are kept: `d3_driver_attempt3_restart2.log` and `d3_driver.log` (scratch). |
| **ER-D** | **Owner ruling §6.** | A **driver guard** makes the driver **stop after the PRIMARY (P1 + P2)** for owner review. The extensions (P3) and the control (P4) run only if `CARD01_RUN_EXTENSIONS=1` is set, which requires a new owner instruction. |
| **ER-E** | **Owner ruling §3–4.** | The analyzer now reports all three start χ² values, the best value, the spread and the `> 0.2` flag for every point. It applies the frozen `decide()` mechanically and labels the result **NUMERICALLY UNRESOLVED / NO ACCEPTED SCIENTIFIC VERDICT** if a load-bearing minimization is flagged. |

## Raw values inspected before this log (disclosure)

| Fit | χ²_eff | Notes |
|---|---|---|
| PRIMARY ε = 0 (ΛCDM null), start 0 | 10977.473 | Single-start smoke test, 4 threads, 473 s. H0 = 68.03, ω_c = 0.1179, BAO χ² 13.75 |
| PRIMARY ε = −1, starts 0/1/2 | 11116.361 / 11114.438 / 11119.446 | **Spread 5.0** (flagged). H0 ≈ 82, BAO χ² ≈ 135–140 |
| ARCHIVED, EXCLUDED: free-ε, seeded at 0, start 0 | 10974.070 | ε = −0.0510; used for **no statistic** |
| ARCHIVED, EXCLUDED: free-ε, seeded at 0, start 2 | 10976.550 | ε = −0.0225; used for **no statistic** |

The owner has ruled that these preliminary numbers carry no scientific interpretation.
