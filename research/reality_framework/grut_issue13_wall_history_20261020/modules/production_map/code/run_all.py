"""production_map / run_all: run a1, b1, a2 in order, then write results/headline.json with assertions and hashes.

  cd .../modules/production_map/code
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 run_all.py
"""
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
env = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2")
times = {}
if "--no-run" in sys.argv:
    try:
        times = json.load(open(os.path.join(HERE, "..", "results", "headline.json")))["runtimes_s"]
        times["note"] = "runtimes of the last full run; this headline was rebuilt with --no-run plus the fast scripts"
    except Exception:
        times = {}
else:
    for script in ("a1_restoration.py", "b1_spectator.py", "b2_spectator_bath.py", "b3_isocurvature.py",
                   "a2_producer.py"):
        t = time.time()
        subprocess.run([sys.executable, os.path.join(HERE, script)], check=True, env=env, cwd=HERE,
                       stdout=subprocess.DEVNULL)
        times[script] = time.time() - t

sys.path.insert(0, HERE)
from common import load, dump, TC_RESUMMED, Y_MAX_STRICT  # noqa: E402

a1 = load("a1_restoration.json")
b1 = load("b1_spectator.json")
a2 = load("a2_producer.json")
b2 = load("b2_spectator_bath.json")
b3 = load("b3_isocurvature.json")
A = []


def check(name, ok, value):
    A.append(dict(name=name, ok=bool(ok), value=value))


ch = a1["checks"]
check("Hartree Tc(z=1) inside resummed range", TC_RESUMMED[0] <= ch["Tc(z=1)/v Hartree [C]"] <= TC_RESUMMED[1],
      ch["Tc(z=1)/v Hartree [C]"])
check("xi_kin(phi) at z=1 within 5% of 6.6e-5 (thermal_potential)",
      abs(ch["xi_kin(phi) at z=1 [C] vs 6.6e-5 [Q thermal_potential]"] / 6.6e-5 - 1) < 0.05,
      ch["xi_kin(phi) at z=1 [C] vs 6.6e-5 [Q thermal_potential]"])
tmin = a1["minimum T_SM of a relic-compatible restored phase (z = 1) [C]"]["strict 4.373e-12 [Q v_relic]"]
check("T_SM,min (strict Y) in 2.0-2.2e6 GeV", 2.0e6 < tmin["T_SM_min_Hartree_GeV"] < 2.2e6,
      tmin["T_SM_min_Hartree_GeV"])
mb = [r for r in a1["exit grid T_SM,end(Y_N, eps) [C]"] if r["z"] < 1e-6]
check("MB closed form = Hartree/Bose exit for z < 1e-6 (to 1e-3)",
      all(abs(r["T_SM_end_MB_closed_form"] / r["T_SM_end_GeV"] - 1) < 1e-3 for r in mb), len(mb))
check("Tc(z) -> leading sqrt(zeta2/Li2) as z -> 0 (to 1e-6)",
      abs(ch["Tc(z)->Tc_leading*sqrt(zeta2/Li2(z)) as z->0 (thermal masses vanish) [C]"] - 1) < 1e-6,
      ch["Tc(z)->Tc_leading*sqrt(zeta2/Li2(z)) as z->0 (thermal masses vanish) [C]"])
cal = b1["calibration single quartic field [C]"]
check("single quartic field: J (zero-crossing) = exact adiabatic invariant to 1%",
      all(abs(x - 1) < 0.01 for x in cal["J_zc_over_exact"]), cal["J_zc_over_exact"])
dyn = b1["conformal dynamics, instant reheating, minimal couplings (units H_e) [C]"]["with slow-roll velocity"]
check("conformal energy conserved (Verlet) to 1e-4", dyn["energy_conservation"] < 1e-4, dyn["energy_conservation"])
inst = b1["instant reheating, minimal couplings [C]"]
e_an = inst["eps = coef (H/M_P)^2: analytic sum<chi^2>/6 + <U4>/3 [D]"]
e_sm = inst["eps coef from sampled initial data incl. slow roll [C]"]
check("eps coefficient: sampled vs analytic within 3%", abs(e_sm / e_an - 1) < 0.03, [e_sm, e_an])
check("overpopulation J/n_eq > 1", inst["overpopulation J/n_eq(rho) [C]"] > 1, inst["overpopulation J/n_eq(rho) [C]"])
hh = a2["check vs hidden_heating c5 row (H_e=1e9, m_I=3e9, Gamma_I = H(1e8 GeV)) [C]"]
check("producer B_min (Hubble-damped) reproduces hidden_heating c5 7.22e-14 within 10%",
      abs(hh["B_min_eff_free"] / 7.219962118404513e-14 - 1) < 0.10, hh["B_min_eff_free"])
check("producer C in [0.9, 1.0]", 0.9 <= a2["headline"]["C (Y_N = 1.5 C B T_RH/m_I)"] <= 1.0,
      a2["headline"]["C (Y_N = 1.5 C B T_RH/m_I)"])

check("quartic FP Lambda1/sqrt(lambda) = 0.08892 (b3)",
      abs(b3["check: Lambda1/sqrt(lambda) for quartic = 0.08892 [R]"] / 0.08892 - 1) < 1e-3,
      b3["check: Lambda1/sqrt(lambda) for quartic = 0.08892 [R]"])
hashes = {}
for f in sorted(os.listdir(HERE)):
    if f.endswith(".py"):
        hashes[f] = hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()

head = {
    "module": "production_map",
    "Y_N_max_strict [Q v_relic]": Y_MAX_STRICT,
    "A1": a1["headline"],
    "A2": a2["headline"],
    "B": b1["headline"],
    "B2 spectator bath: H_inf for kinetic exit / zero-mode response": [
        b2["H_inf above which a mu=0-thermalized spectator bath is kinetically thermal at its exit (xi_E >= 6.6e-5)"],
        b2["H_inf above which the gas-bound zero mode responds (A > 200 GeV)"]],
    "B3 isocurvature [E, R]": {k: dict(P_S=r["P_S_total"], Y_N_max_iso=r["Y_N_max_iso"], H_thr_floor=r["H_thr_floor"],
                                      H_thr_J=r["H_thr_J"]) for k, r in b3["rows"].items()},
    "assertions": A,
    "all_assertions_pass": all(x["ok"] for x in A),
    "runtimes_s": times,
    "sha256": hashes,
}
dump("headline.json", head)
print(json.dumps({"all_assertions_pass": head["all_assertions_pass"],
                  "failed": [x for x in A if not x["ok"]], "runtimes_s": times}, indent=1, default=float))
