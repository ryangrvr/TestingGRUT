# N22-01 same-rule detector refinement: result

Protocol frozen before production at `0f2d2e4ef53f4bcc2613f5cb1373603c4551d9b4` on 9 November 2026. This follows the [research note](https://github.com/ryangrvr/TestingGRUT/pull/14#issuecomment-6733231617).

**Verdict: finite_massless_pass_with_numerical_guard.** Execution completed: True. Audit passed: True. Massive transfer: withheld under every N22-01 verdict. B2's N20 combined measured budget remains a pass (0.000887384268 < 0.001), its separate 0.000740 lattice allowance remains a miss, and its N20 massive transfer remains withheld.

The geometry coefficients were recomputed at N22 with the identical B2 reconstruction rule; the N20 coefficient array was not reused. B2 model.py and point.run were imported without edits. Length, Gaussian, coupling, mass convention, preparation rule, Hamiltonian, boundaries and all seven timings were retained. These are inherited timing checks, not fresh blind timings for N22.

N22 meets the registered finite massless numerical requirements, including the extra component margin guard. This is resolution-specific progress; it does not change N20 or remove the operator-dictionary and interacting-parent obligations.

## Frozen decisions

All nominal gates passed: **True**. Component plus full numerical reserve passed: **True**. Failed nominal gates: none.

The worst matched TV across both solves is 0.000729406759042; the nominal limit is 0.000740. Its nominal headroom is 1.0593240958e-05. Adding the frozen 0.000010 numerical guard leaves 5.93240957993e-07 headroom. The guard is an operational numerical precaution, not a proved bound on all numerical or physical errors. A positive margin is not a rigorous continuum certificate.

The largest direct target TV is 0.000730900161437. The largest matched+geometry+reserve total is 0.000740900161437, against 0.001. Maximum geometry TV is 6.08109054548e-06, against 0.00025.

| Times | Standard matched TV | Refined matched TV | Worst geometry TV | Worst direct TV | Worst additive total | Margin after reserve |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.0002170393747 | 0.0002170393747 | 4.865349938e-07 | 0.0002175259097 | 0.0002275259097 | 0.0005129606253 |
| 0,1 | 0.000265635662 | 0.000265635662 | 3.126581462e-06 | 0.0002636212652 | 0.0002787622435 | 0.000464364338 |
| 0,2 | 0.000729406759 | 0.000729406759 | 1.493402395e-06 | 0.0007309001614 | 0.0007409001614 | 5.93240958e-07 |
| 0,1,2 | 0.0007080509354 | 0.0007080509354 | 4.864652149e-06 | 0.0007058884767 | 0.0007229155876 | 2.194906456e-05 |
| 0,0.6 | 0.0002184676369 | 0.0002184676369 | 4.38965718e-06 | 0.0002175259098 | 0.0002328572941 | 0.0005115323631 |
| 0,1.7 | 0.0006507092847 | 0.0006507092847 | 1.901207893e-06 | 0.0006524370247 | 0.0006626104926 | 7.929071528e-05 |
| 0,0.6,1.7 | 0.0005826566253 | 0.0005826566253 | 6.081090545e-06 | 0.00058006057 | 0.0005987377158 | 0.0001473433747 |

## Same-size refinement and resources

Both points use the 705432-dimensional neutral sector. The standard solve uses tolerance 1e-11. The refined solve uses 1e-13, a different starting-vector seed, and split propagation intervals. Maximum standard-to-refined history TV: 1.28340914285e-13. Maximum N20-to-N22 TV across both solves: 0.000145990263031.

| Solve | Seconds | Peak job-member private memory (GiB) | Eigensolver residual | Computed low-eigenvalue gap |
|---|---:|---:|---:|---:|
| N22 | 401.453 | 1.4412 | 2.25298e-10 | 0.842006737401 |
| N22_refined | 449.438 | 1.4840 | 2.70836e-13 | 0.842006737401 |

Total registered execution: 850.969 seconds. The limits were 900 seconds per point, 2000 seconds total and 4 GiB committed memory per child process. Job peak memory includes descendant processes. The periodic heartbeat sampled the Windows virtual-environment launcher, so its small working-set values should not be mistaken for solver memory; use the job peak above. Additional worker observations are retained separately. No automatic extension or retry was performed.

Maximum reference mode-doubling TV: 8.24277868183e-11. Maximum analytic history-tail bound evaluated in floating point: 5.52895267289e-07.

The computed eigenvalue gap is reported as a diagnostic. Neither it nor the residual provides a certified spectral-separation argument or a rigorous ground-state error bound.

## Outcome-resolved scaling forecast

All signed-vector forecast comparisons accepted: **True**. The criterion, fixed before N22, was at most 10% relative half-L1 change between delta22/a22^2 and delta20/a20^2, separately for every history and solve. This is a diagnostic tolerance, not a derived law. The forecast was never subtracted from measured probabilities.

| Times | N20-based predicted TV at N22 | Standard relative signed-vector change | Refined relative signed-vector change |
|---|---:|---:|---:|
| 0 | 0.0002159327423 | 0.0051248942 | 0.0051248942 |
| 0,1 | 0.0002603849738 | 0.020165097 | 0.020165097 |
| 0,2 | 0.0007295070061 | 0.0016543768 | 0.0016543769 |
| 0,1,2 | 0.0007067576134 | 0.0087909384 | 0.0087909385 |
| 0,0.6 | 0.0002159327423 | 0.02429205 | 0.02429205 |
| 0,1.7 | 0.0006504255634 | 0.0017013975 | 0.0017013975 |
| 0,0.6,1.7 | 0.0005757885547 | 0.02026545 | 0.02026545 |

The full signed outcome vectors and vectors divided by a^2 are preserved for N10,14,18,20 and both N22 solves in assessment.json. Agreement of this one refinement with the forecast does not prove a continuum rate; histories have different finite-spacing trends.

## Checks, interpretation and next work

Preflight algebra/oscillator/N6 checks passed. Normal child execution, memory-limit rejection and timeout termination were tested before freeze. Synthetic tests distinguish nominal failure, nominal pass with unresolved margin, and a guarded pass; a sign-reversed forecast with unchanged TV is rejected.

The post-run audit reports: True. It verifies source/input hashes, distribution structure and normalization, independently recomputes TVs, signed vectors and forecast ratios, checks the saved ground states and first-event formula, and replays the frozen decision. It is an arithmetic/integrity audit using the inherited Hamiltonian and weights, not an independent full-size physical implementation.

An N22 numerical pass supports the registered finite massless experiment under the adopted lattice-to-cell-average dictionary. The physical operator identification remains an assumption requiring separate derivation or evidence. A finite miss would also need to be distinguished from a combined-total violation or a falsification of GRUT.

N22-01 contains no massive execution path. It cannot retroactively pass N20 or clear B2's withheld N20 transfer. Any later massive comparison needs a separately frozen protocol specifying its own lattice resolution and dictionary; P1's interacting-parent cutoff and overflow obligations remain open. No merge, release or Zenodo deposit was performed.

## Archive and reproduction

Source, protocol, raw JSON, telemetry, audit and checkpoint hashes are saved on codex/grut-n22-refinement-20261109. The delivered ZIP also contains the binary ground states and all source/input files needed for the offline audit. CHECKPOINTS.json identifies those binaries; SHA256.json covers every payload file except itself. The ZIP preserves the exact Windows source bytes whose hashes were recorded during execution.

Use Python 3.12, numpy==2.3.5 and scipy==1.17.0. To reproduce production, check out the source freeze in a clean Git worktree on Windows and follow PROTOCOL.md. To replay the archive audit, preserve the original packet and use a copy of raw without its existing audit.json, then run:

    python audit_n22.py --root source --output RAW_COPY

The audit refuses to overwrite its output. The source archive omits .git; production intentionally requires a clean committed Git checkout.
