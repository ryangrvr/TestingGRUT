# Reproduce the local controls

Requirements used: Python 3.12, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib. No dependency installation or external numerical service was used. Run at the repository root. The numerical workers are ordinary local processes, not separate AI authors.

Fetch the pinned completed Claude builds if absent:

```bash
git fetch origin e71730c24fe6a0b28ddd72cd30982c5b39c9aedf
```

The bootstrap extracts exact input files to a temporary work directory using git archive, leaving both commits untouched. Optional GRUT_KINETIC_INPUT_ROOT points directly at the pinned grut_formation_review_20261025 directory; GRUT_KINETIC_WORKDIR selects scratch caches. A cached conversion table uses only the certified cold forward data; its high-momentum hold recipe is an assumption, not a proven bound.

```bash
python -B research/reality_framework/grut_kinetic_continuation_20261026/code/relax_control.py
python -B research/reality_framework/grut_kinetic_continuation_20261026/code/production.py
python -B research/reality_framework/grut_kinetic_continuation_20261026/code/followup.py
python -B research/reality_framework/grut_kinetic_continuation_20261026/code/verify.py
python -B research/reality_framework/grut_kinetic_continuation_20261026/code/summarize.py
python -B research/reality_framework/grut_kinetic_continuation_20261026/code/plot.py
python -B research/reality_framework/grut_kinetic_continuation_20261026/code/audit.py --seal
```

Production uses Radau rtol=2e-9 and an active-set susceptibility ODE. The follow-up includes a BDF comparison and looser Radau controls. There are 13 production cases, 11 follow-ups and four conditional hard cascades. Object cap and energy differ deliberately: omega sets the supplied chemical cap, E/Q books fixed-Q object energy. Their physical availability is not derived.

Use moment_controls.py for the actual production method. stage2_moments.py and thermo_state.py are credited ports of Claude's unfinished code, with corrections recorded in REPORT.md. The legacy standalone stage2 integrator is disabled; only its rate/capture/Case helpers are used by production. Its retained methods are not separately certified.

All moment starts are forced full chemical/kinetic equilibrium. The hard contact control does not establish that start. No initial-spectrum bracketing, nonlinear formation, physical phase-transition matching or late nonthermal relic solver is hidden in the commands. The outputs deliberately leave the physical f_sym and kinetic threshold null.

Wall times and CPU call counts are diagnostic and need not reproduce byte-for-byte. audit.py checks source hashes, actual saved numerical checks, finite JSON, syntax, links and manifest. After regenerating data, --seal creates a new manifest; running audit.py without --seal verifies the supplied one.
