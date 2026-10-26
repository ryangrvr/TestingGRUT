# Reproduction and data conventions

Run from the repository root. Python3.12, NumPy2, SciPy and GCC are used;
Matplotlib produces the scientific figure. Parent e160356 and pinned
review f74fb83 are the source context. No reviewer solver is imported or
run by this packet. Exact source hashes are in SOURCES.json.

```bash
git fetch origin f74fb83071a952e888e507733aeacbfc1b4ad8b9
export GRUT_MODE_WORKDIR=/tmp/grut_reheating_modes_20261026
python -B research/reality_framework/grut_reheating_modes_20261026/code/modes.py
python -B research/reality_framework/grut_reheating_modes_20261026/code/verify.py
python -B research/reality_framework/grut_reheating_modes_20261026/code/refine.py
python -B research/reality_framework/grut_reheating_modes_20261026/code/phase_and_orbit.py
python -B research/reality_framework/grut_reheating_modes_20261026/code/final_checks.py
python -B research/reality_framework/grut_reheating_modes_20261026/code/neutral.py
python -B research/reality_framework/grut_reheating_modes_20261026/code/plot.py
python -B research/reality_framework/grut_reheating_modes_20261026/code/audit.py --seal
python -B research/reality_framework/grut_reheating_modes_20261026/code/audit.py
```

The scripts use up to three ordinary CPU workers, without AI subagents.
Independent numerical jobs may run in parallel after modes.py; dependencies
in the sequence above must be preserved. Unique temporary build directories
prevent concurrent compiler output collisions. DOP853 coefficient arrays
come from the installed SciPy distribution; Radau is the distinct short
verification method. No SymPy or external solver download is needed.

Raw binary transfer arrays are in the work directory and can be recreated.
They contain each sampled tau,19 background variables and one row-major
12-by12 transfer matrix per charged-system k. Neutral mode blocks are
2-by2. Codes preserve compact review arrays in the repository:

- MODE_PROFILES.npz: original times,39 ks, six backgrounds and eight
  diagnostic columns named in MODE_RESULTS.json. Selected exact transfer
  matrices use times selected from0,1,3,10,40,80 and k indices
  0,2,6,14,22,30,38 of the original grid.
- REFINED_PROFILES.npz: the65-node .8–6 band and the same eight columns.
- PHASE_PROFILES.npz: merged original/refined ks per background, signed
  amplitude and charge response per supplied phase radian. These displayed
  curves are float32; summaries and the underlying integrations are double.
- NEUTRAL_PROFILES.npz:136 ks, two conditional occupations and omega, plus
  the two conditional excitation-energy ratios. Occupation curves are
  float32, ratio arrays and summaries double.
- QUARTIC_ORBIT.json: matched orbit, radial closing angle/period and the
  auxiliary growth-band eigenvalue moduli.

Seven verification groups are saved in VERIFICATION_RESULTS.json. The
instantaneous and first-order WKB phi seed controls are conditional; their
particle projection is not a thermal or interacting abundance calculation.
The generic six-field Gaussian reference contains zero-point fluctuations;
its short-wavelength charge variance must not be interpreted as relic mass.

The first neutral invocation stopped because the input contained136 modes
and the generic reader limit was128. The limit was raised to256, without
changing the grid, equations or tolerance, and the calculation completed.
The initial GitHub metadata request used a short commit identifier on an
endpoint requiring a full one; the normal commit endpoint resolved it.
These were execution/packaging corrections, not changes to a physical
outcome. Build-directory isolation was improved without changing equations.

audit.py --seal refreshes the packet hashes after actual controls and
source/JSON/syntax/link checks pass. Plain audit.py verifies the saved
manifest. Cross-platform numerical runs need comparison by the recorded
physical controls and tolerances; bitwise identity is not assumed.
