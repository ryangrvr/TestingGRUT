"""Frozen deterministic family list, original root weights, complete attempt log."""
from fractions import Fraction
from pathlib import Path
import json
import time

import active as A
import controls as C

BASE = Path(__file__).resolve().parent


def families():
    out = []
    for n in [4, 5, 6]:
        for bar in [0.0, 0.5]:
            values = [bar-.25, bar+.25]
            out.append(dict(name=f"balanced_N{n}_bar{bar}", role="PRIMARY", N=n,
                            values=values, global_counts=[(n+1)//2, n//2], g=.9))
    out.append(dict(name="outlier_N5", role="SECONDARY", N=5,
                    values=[-.25, .25], global_counts=[4, 1], g=.9))
    out.append(dict(name="ramp_N4", role="SECONDARY", N=4,
                    values=[float(Fraction(-1,4)), float(Fraction(-1,12)),
                            float(Fraction(1,12)), float(Fraction(1,4))],
                    global_counts=[1]*4, g=.9))
    return out


def components(family):
    out = []
    for root, count in enumerate(family["global_counts"]):
        remaining = family["global_counts"].copy()
        remaining[root] -= 1
        values, counts = [], []
        for v, c in zip(family["values"], remaining):
            if c:
                values.append(v)
                counts.append(c)
        out.append(dict(root_type=root, weight=count/family["N"],
                        counts=counts, values=values, g=family["g"]))
    return out


def key(counts, values, g):
    return json.dumps([list(counts), list(values), g], separators=(",", ":"))


def mirror_key(counts, values, g):
    pairs = sorted((-v, c) for v, c in zip(values, counts))
    return key([c for _, c in pairs], [v for v, _ in pairs], g)


def run(use_pilot=None):
    plan = dict(base_commit="c6ac0c68480d4164ad620bbf2a4b3d8c65ac85b1",
                frozen_charter_commit="64378573e8251224dd563928667e961a0effb3da",
                observation="4*N^2", strang_budget=.04, projection_budget=.15,
                rank_ceiling=24, coefficient_ceiling=2_000_000,
                per_instance_wall_seconds=120, degree=8,
                root_weights="original global class count / N", families=families())
    (BASE/"RUN_PLAN.json").write_text(json.dumps(plan, indent=2)+"\n")
    output = dict(plan=plan, checks={}, controls=[], attempts=[], families=[],
                  all_size_formation_proved=False, external_validation=False,
                  roundoff_certified=False)
    def checkpoint():
        temporary = BASE/"RESULTS.json.partial"
        temporary.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
        temporary.replace(BASE/"RESULTS.json")
    for name in ["algebra_checks", "commutator_checks", "complete_wavefunction_checks"]:
        output["checks"][name] = getattr(C, name)()
        print("CHECK", name, "PASS", flush=True)
    # Actual full-physical-receiver uniform positive controls.
    for n, delta, g in [(4, 0.0, .9), (5, .5, 1.1)]:
        evolution = A.Evolution([n-1], [delta], 64*n*n, g)
        r = evolution.run(4*n*n)
        r["quadratic_reference"] = C.quadratic_reference(evolution, 4*n*n, delta)
        assert r["status"] == "COMPLETE"
        assert r["quadratic_reference"]["state_error"] <= r["analytic_state_error_account"]+1e-8
        output["controls"].append(r)
        print("UNIFORM", n, r["analytic_probability_lower"], flush=True)
    seen = {}
    for family in plan["families"]:
        row = dict(**family, components=[])
        for component in components(family):
            counts, values, g = component["counts"], component["values"], component["g"]
            direct, mirror = key(counts, values, g), mirror_key(counts, values, g)
            relation = "direct"
            if direct in seen:
                case = seen[direct]
                relation = "same root component reused"
            elif mirror in seen:
                case = seen[mirror]
                relation = "proved global-sign/type-relabeling identity"
            else:
                index = len(output["attempts"])
                n = sum(counts)+1
                print("START", index, family["name"], counts, values, flush=True)
                # An earlier bounded pilot used this identical numerical
                # algorithm; metadata-only fields are completed below.
                if use_pilot and counts == [2,2] and values == [-.25,.25] and g == .9:
                    case = json.loads(Path(use_pilot).read_text())
                    case["completed_physical_steps_time"] = case["completed_steps"]*case["dt"]
                    case["exponential_degree"] = 8
                    case["analytic_full_assembly_lower"] = None
                    case["incomplete_population_notice"] = "Populations are intermediate only, not observations at the frozen target time"
                    case["execution_note"] = "Reused earlier bounded S5 pilot; numerical routine unchanged, metadata added"
                else:
                    evolution = A.Evolution(counts, values, 64*n*n, g)
                    case = evolution.run(4*n*n, max_seconds=plan["per_instance_wall_seconds"],
                        progress=lambda x: print("PROGRESS", index, json.dumps(x), flush=True))
                case["case_index"] = index
                seen[direct] = case
                output["attempts"].append(case)
                print("FINISH", index, case["status"], case["analytic_probability_lower"],
                      case["stop_reason"], flush=True)
                checkpoint()
            row["components"].append(dict(**component, case_index=case["case_index"], relation=relation))
        weights = sum(c["weight"] for c in row["components"])
        assert abs(weights-1) < 1e-12
        lower = full_lower = unresolved = 0.0
        for c in row["components"]:
            case = output["attempts"][c["case_index"]]
            if case["status"] == "COMPLETE":
                lower += c["weight"]*case["analytic_probability_lower"]
                full_lower += c["weight"]*case["analytic_full_assembly_lower"]
            else:
                unresolved += c["weight"]
        row.update(root_averaged_formation_lower=lower,
                   root_averaged_full_assembly_lower=full_lower,
                   unresolved_root_weight=unresolved,
                   complete_root_ensemble=bool(unresolved == 0),
                   primary_gate_pass=bool(unresolved == 0 and lower >= .5),
                   status="COMPLETE" if unresolved == 0 else "UNRESOLVED",
                   unknown_root_lower_bound_used=0.0)
        nu, eps, t = 2.1, .25, .005
        rmin = max(1, (family["N"]+1)//2-2)
        # S4 exact formula, evaluated at the frozen near-response time.
        from math import exp, factorial
        tau = sum((nu*t)**k/factorial(k) for k in range(rmin, 50))
        d = eps/nu*tau
        kappa = .25*t-(exp(2*nu*t)-1-2*nu*t)
        row["joined_protocol"] = dict(probe_time=t, source_variance=.25,
              all_macro_records_accepted=True, actual_receiver_retained=True,
              r_min=rmin, S4_leakage_amplitude_upper=d,
              joint_accept_and_consistent_probability_lower=lower*(1-min(1,d)**2),
              near_response_per_original_trial_lower=lower*kappa,
              scope="finite family only; no all-size independent-disorder inference")
        output["families"].append(row)
        checkpoint()
        print("FAMILY", family["name"], row["status"], lower, "unresolved weight", unresolved, flush=True)
    checkpoint()
    return output


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse-pilot")
    arguments = parser.parse_args()
    run(arguments.reuse_pilot)
