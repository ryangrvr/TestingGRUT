# Reproduction

Run from repository root at the packet commit or a descendant containing
its pinned inputs. No external review checkout is imported or executed by
this packet. The original cold-rate and vector outputs are reused as
explicit inputs; this is not a rerun of those larger calculations.

Recorded environment: Python3.12.14, NumPy2.3.5, SciPy1.17.0,
GCC13.3.0 on Linux. Matplotlib is needed for the figure. The C solver uses
DOP853 coefficient arrays generated from the installed SciPy distribution;
it does not download a solver or require SymPy. NumPy2 supplies trapezoid.

```bash
export GRUT_REHEATING_WORKDIR=/tmp/grut_reheating_preparation_20261026
python -B research/reality_framework/grut_reheating_preparation_20261026/code/reheating.py
python -B research/reality_framework/grut_reheating_preparation_20261026/code/verify.py
python -B research/reality_framework/grut_reheating_preparation_20261026/code/boundary.py
python -B research/reality_framework/grut_reheating_preparation_20261026/code/portals.py
python -B research/reality_framework/grut_reheating_preparation_20261026/code/plot.py
python -B research/reality_framework/grut_reheating_preparation_20261026/code/audit.py --seal
python -B research/reality_framework/grut_reheating_preparation_20261026/code/audit.py
```

The first command compiles the C solver and runs 11 standalone histories
using up to three ordinary CPU worker threads. The verification command
adds two C histories and one SciPy Radau comparison; it requires the raw
work directory created by the first command. These are numerical jobs, not
AI subagents. Portal sources and the finite boundary model can run without
the raw reheating directory. The figure uses the committed compact profiles.

`audit.py --seal` records the successfully reproduced output hashes. Do not
expect bitwise equality across numerical-library/compiler versions. First
compare physical values, controls and residuals to the committed output.
The manifest records this run's bytes, not a claim of universal bitwise
reproducibility. Plain `audit.py` checks the saved packet and its hashes
without repeating the expensive integrations.

Primary outputs:

- REHEATING_RESULTS.json: six production histories and five controls,
  parameter/preparation dictionaries, output summaries and residuals.
- INDEPENDENT_REHEATING.json: short Radau/C endpoints and long control.
- REHEATING_VERIFICATION.json: four geometry/integration groups.
- REHEATING_PROFILES.npz: compact early and late diagnostics, with column
  names in REHEATING_RESULTS.json; all dimensional conventions are explicit.
- BOUNDARY_RESULTS.json: two finite cutoffs, phases, constraints and four
  state/observation groups.
- PORTAL_RESULTS.json: 19 source diagnostics, refinement, analytic checks,
  widths, matching, scattering and the inherited included cold leakage.
- VERIFICATION_RESULTS.json: the 11 mathematical groups and the physical
  claims they do not establish.

The figure was rendered and visually inspected. Protocol choices were not
retuned after seeing yields or source failures. A NumPy-boolean console
serialization error in the verification script was corrected after the
scientific results had been saved; the corrected script was rerun. No
failed physical outcome was changed by that correction.
The packaging audit initially checked its own output link before creating
that output. Its write/check order was corrected; the completed audit then
passed. No scientific calculation or tolerance was changed.
