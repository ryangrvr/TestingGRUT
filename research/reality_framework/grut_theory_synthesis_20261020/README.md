# GRUT working theory and derived predictions

**Subsequent test:** [the reality audit](../grut_reality_audit_20261020/REPORT.md) supplies ten full-field evolutions with a separate Fourier/RK4 solver, physical unit conversions and an empirical-access assessment. [START_HERE.md](../grut_reality_audit_20261020/START_HERE.md) preserves the continuation checkpoint.

Read [THEORY.md](THEORY.md) for the complete candidate and [CLAUDE_HANDOFF.md](CLAUDE_HANDOFF.md) for the next controlled constructions. This is an explicit working EFT synthesis with a demonstrated classical sector; it is not an established theory of everything.

New deductions/checks under the existing reference scalar action:

- A sharp classical energy/charge bound and the energetic binding boundary a*b<1.
- Bulk coexistence density, frequency and equation of state.
- The full dense-phase linear spectrum and sound speed c/sqrt(3).
- A specified coherent mixed-branch instability and its predicted band/growth.
- A resolved two-field planar interface with surface cost about .680013489 and the corresponding large-Q surface correction.
- Exact continuum symmetry tangents that clarify how to test the earlier unconverged m=1 response.

All 10 author calculation checks pass. The scalar mechanism is in the known FLS/Q-ball family; deriving these formulas here does not establish historical novelty. The new bulk predictions have not yet been tested in full nonlinear wave or fragmentation evolutions. The prior assembly, persistence and charge-record simulations remain the [separate completed packet](../grut_candidate1_two_object_20261019/persistence/REPORT.md).

![Binding, interface and linear response](THEORY_FIGURE.svg)

Run `python3 derive_and_check.py` to reproduce [DERIVATION_RESULTS.json](DERIVATION_RESULTS.json) and [wall_profile.json](wall_profile.json). NumPy and SciPy are required. `python3 plot.py` creates the SVG and a local PNG with Matplotlib. The results record runtime versions and the calculation-source SHA-256. The interface calculation is checked for domain/tolerance sensitivity, not certified globally minimal. No outside data or parameter fitting is used.
