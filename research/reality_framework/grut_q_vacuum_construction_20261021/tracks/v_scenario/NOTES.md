# Track `v_scenario` (verifier of `../scenario/`)

The record is `VERIFY.md`:
- §1 is the verdict table: 33 claims, with the scenario's number, my number and the verdict.
- §2 gives the equations and the method of each check.
- §3 covers the audit of conventions, calibration, credit and consistency.
- §4 lists the inputs with provenance.
- §5 gives the reproduction command and the limits.

`results/verdicts.json` holds the same table in machine-readable form.

Reproduction: `python3 -I code/run_all.py`. Downloads go only to the scratchpad `qtheory/v_scenario_dl`; the run
takes about 6 minutes on 4 cores.

Corrections, in short:
1. The QV1-cal second-order term is c/(K0 m²) = 1.08e-120, not (4/3)c/(K0 m²) = 1.4e-120. The μ offset is not the
   fixed-μ algebra.
2. The 10 TeV threshold at fixed μ turns around at 4.51e-13 s and crunches at 9.01e-13 s, not at 4.0e-13 s and
   ~8e-13 s.
3. The deposit at x = 8, m/H = 100 is 1.077 ± 0.005 times the formula, not 1.028. The scenario's floor interferes at
   ±6% there.
4. Under P_q, the 10 TeV scalar gives δG/G = +6.8e-61, not −6.8e-61.
5. The toy dust fraction at N = 0 is 2.05%, not 6%.
6. Clarification: the "7.2e-121" offset is c/M⁴; δμ/μ0 = 2.15e-120.
7. Missing credit for the chameleon/f(R) trace kick: astro-ph/0408415 and arXiv:1304.0009.

No verdict of the scenario track changes.
