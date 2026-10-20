# Assembly from prepared Gaussian charge clouds

The full coupled action is unchanged. Start with charge spread over a Gaussian profile and S=1 everywhere, then evolve without a later drive or damping. Widths 4,6,8 have the same supplied continuum charge. A rigid-S comparison removes mediator feedback explicitly. This tests conditional assembly; it does not derive the charged low-entropy preparation.

Sources and state choices are frozen before these cloud evolutions. `node run.cjs RUN_ID` replays one case. The rigid comparison is an ablation law, not a claim that holding a field fixed is free in a universe.

All 13 cases completed. Medium-grid width-4 and width-6 radius-12 cores retain 81.18% and 55.00% of total charge at tau=240, versus 1.55% for the matched width-6 rigid control. Width 8 remains marginal and has not earned a bound-core conclusion. See [the combined report](../REPORT.md), [summary](../SUMMARY.json), and [exact output archive](../raw_outputs.tar.gz).
