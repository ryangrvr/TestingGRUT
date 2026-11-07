# Clean-checkout validation repair

The original repository entry point used a stale `../qed` relative path, overwrote its archived JSON and did not return a failing exit status when its numerical gate failed. The repaired wrapper resolves paths from its own file, requires a new output outside the archived packet, records reference hashes, and returns 0 for PASS and 1 for FAIL. The original `rel < 1e-5` acceptance rule and DMRG physics are unchanged.

The first clean attempt at frozen commit `2b13287` exposed a new wrapper bug: NumPy's Boolean scalar was not JSON serializable. [The failure log](validation_first_attempt.log) is preserved. Local source commit `29c3675` converts that flag to a Python bool. This was an implementation failure, not a physical discrepancy or a gate to be waived.

The corrected command was run from `/tmp` against a fresh detached checkout of local source commit `29c3675`, with NumPy 2.3.5, SciPy 1.17.0 and physics-tenpy 1.1.1. All four cases pass, with maximum relative gap difference `1.333e-13`. A separate reference copy with the first gap deliberately increased by 10% fails at relative difference `0.0909091` and returns exit status 1. The checkout remains clean. Both JSON results are in [validation_replays](validation_replays/).

These checks verify the wrapper and repeat the existing solver benchmark. They do not constitute an independent DMRG implementation, a full production-grid replay, or a frozen external test of an evaluator's ability to judge unfamiliar theories.

Use the command in the corrected QED `RESULTS.md`. `--case 0` is a diagnostic subset and is labeled as such. A supplied `--reference-dir` is recorded; the corrupted-reference diagnostic is never used as scientific input.
