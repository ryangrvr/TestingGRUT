#!/usr/bin/env python3
"""Run every GRUT III verification script and write ../RESULTS.md.
Dependencies (see requirements.txt): v1-v3 need sympy + mpmath (v1: mpmath+sympy, v2/v3: sympy); v2 and v4 need numpy."""
import subprocess, sys, datetime, pathlib
here = pathlib.Path(__file__).resolve().parent
scripts = ["v1_ds_thermal_kernel.py", "v2_jacobson_eos.py", "v3_horizon_fluid.py",
           "v4_medium_absorptive_kernel.py", "v5_anomaly_selfconsistency.py", "v6_low_acceleration_scale.py",
           "v7_access_geometry_memory.py", "v8_residual_framework.py", "v9_induced_gravity_matching.py",
           "v10_graviton_sector_memory.py", "v11_quasi_ds_kms_defect.py"]
rows, logs, allok = [], [], True
for s in scripts:
    r = subprocess.run([sys.executable, "-I", str(here/s)], capture_output=True, text=True)
    status = "PASS" if r.returncode == 0 else "FAIL"
    if s.startswith("v6"): status = "INFO"
    allok &= (r.returncode == 0)
    rows.append((s, status)); logs.append((s, r.stdout + (("\nSTDERR:\n" + r.stderr) if r.stderr.strip() else "")))
out = [f"# GRUT III — verification results\n", f"Generated {datetime.datetime.now(datetime.UTC).strftime('%Y-%m-%d %H:%M UTC')} by `GRUT_III/calc/run_all.py` (Python {sys.version.split()[0]}; numpy, sympy, mpmath).\n",
       "| script | status |", "|---|---|"] + [f"| `{s}` | **{st}** |" for s, st in rows] + ["", "---", ""]
for s, log in logs:
    out += [f"## `{s}`", "", "```text", log.rstrip(), "```", ""]
(here.parent/"RESULTS.md").write_text("\n".join(out))
print("\n".join(f"{s}: {st}" for s, st in rows)); print("ALL PASS" if allok else "SOME FAILED")
sys.exit(0 if allok else 1)
