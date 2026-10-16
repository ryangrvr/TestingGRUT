"""S5-P1: same physical formation test, complementary-event certificate."""
from math import exp, factorial, sqrt
from pathlib import Path
import hashlib
import json

import active as A
import campaign as P

BASE = Path(__file__).resolve().parent


def lower_bound(p, error):
    p = min(1.0, max(0.0, p))
    reverse = max(0.0, sqrt(max(0.0, p))-error)**2
    complementary = max(0.0, 1-(sqrt(max(0.0, 1-p))+error)**2)
    return max(reverse, complementary)


def run():
    original_bytes = (BASE/"RESULTS.json").read_bytes()
    original = json.loads(original_bytes)
    assert len(original["families"]) == 8
    assert original["extra_checks"]["physical_receiver_reference"]["complete_wavefunction_error"] <= original["extra_checks"]["physical_receiver_reference"]["analytic_state_error_account"]
    output = dict(charter_commit="364781b841121e3169041b121e661048fd793524",
                  state_error_envelope=.60, strang_budget=.20, projection_budget=.35,
                  probability_gate=.5, target_law_changed=False,
                  seed_changed=False, observation_time_changed=False,
                  disorder_test_list_changed=False, original_S5_results_preserved=True,
                  attempts=[], families=[], all_size_formation_proved=False,
                  roundoff_certified=False, external_validation=False)
    output["original_S5_results_sha256"] = hashlib.sha256(original_bytes).hexdigest()
    def checkpoint():
        temporary = BASE/"PROBABILITY_RESULTS.json.partial"
        temporary.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
        temporary.replace(BASE/"PROBABILITY_RESULTS.json")
    seen = {}
    for family in original["families"]:
        if family["role"] != "PRIMARY":
            continue
        row = {k:v for k,v in family.items() if k in ("name","role","N","values","global_counts","g")}
        row["components"] = []
        for component in family["components"]:
            counts, values, g = component["counts"], component["values"], component["g"]
            direct, mirror = P.key(counts, values, g), P.mirror_key(counts, values, g)
            if direct in seen or mirror in seen:
                result = seen.get(direct, seen.get(mirror))
                relation = "same/mirror component reused"
            else:
                prior = original["attempts"][component["case_index"]]
                if prior["status"] == "COMPLETE":
                    result = dict(prior)
                    relation = "completed original S5 state reused; no new trajectory"
                    result["source_S5_case"] = prior["case_index"]
                else:
                    n = family["N"]
                    index = len(output["attempts"])
                    print("P1_START", index, family["name"], counts, values, flush=True)
                    evolution = A.Evolution(counts, values, 64*n*n, g)
                    result = evolution.run(4*n*n, strang_budget=.20, compression_budget=.35,
                         max_seconds=120, progress=lambda x: print("P1_PROGRESS", index, json.dumps(x), flush=True))
                    relation = "P1 trajectory; original S5 stop retained"
                    result["source_S5_case"] = prior["case_index"]
                complete = result["status"] == "COMPLETE"
                error = result["analytic_state_error_account"]
                result["P1_probability_lower"] = lower_bound(result["macro_probability_estimate"], error) if complete else None
                result["P1_full_assembly_lower"] = lower_bound(result["full_assembly_probability_estimate"], error) if complete else None
                result["P1_gate_pass"] = bool(complete and error <= .60 and result["P1_probability_lower"] >= .5)
                result["P1_index"] = len(output["attempts"])
                seen[direct] = result
                output["attempts"].append(result)
                print("P1_FINISH", result["P1_index"], result["status"], result["P1_probability_lower"], flush=True)
                checkpoint()
            row["components"].append(dict(**{k:v for k,v in component.items() if k not in ("case_index","relation")},
                                           P1_index=result["P1_index"], relation=relation))
        low = full = missing = 0.0
        margins = []
        for component in row["components"]:
            result = output["attempts"][component["P1_index"]]
            if result["status"] != "COMPLETE":
                missing += component["weight"]
                continue
            low += component["weight"]*result["P1_probability_lower"]
            full += component["weight"]*result["P1_full_assembly_lower"]
            p, error = result["macro_probability_estimate"], result["analytic_state_error_account"]
            admissible_error = max(sqrt(p)-sqrt(.5), sqrt(.5)-sqrt(max(0.0,1-p)))
            margins.append(admissible_error-error)
        row.update(status="COMPLETE" if missing == 0 else "UNRESOLVED",
                   unresolved_root_weight=missing, root_averaged_formation_lower=low,
                   root_averaged_full_assembly_lower=full,
                   P1_gate_pass=bool(missing == 0 and low >= .5 and all(
                       output["attempts"][c["P1_index"]]["analytic_state_error_account"] <= .60
                       for c in row["components"])))
        if missing == 0 and margins and min(margins) > 0:
            n, t = family["N"], 4*family["N"]**2
            row["finite_continuous_comparison_box"] = dict(
                per_site_detuning_radius=min(margins)/(4*t*(n-1)),
                common_g_radius=min(margins)/(4*t),
                intersect_with_frozen_class=True,
                supports_distinct_detunings=True,
                all_size_nonshrinking_claim=False)
        t = .005
        rmin = max(1, (family["N"]+1)//2-2)
        tau = sum((2.1*t)**k/factorial(k) for k in range(rmin,50))
        d = .25/2.1*tau
        kappa = .25*t-(exp(4.2*t)-1-4.2*t)
        row["joined_protocol"] = dict(probe_time=t, actual_receiver_retained=True,
            all_macro_records_accepted=True, source_variance=.25, r_min=rmin,
            joint_accept_and_consistent_probability_lower=low*(1-min(1,d)**2),
            near_response_per_original_trial_lower=low*kappa,
            scope="finite instance only")
        output["families"].append(row)
        checkpoint()
        print("P1_FAMILY", family["name"], row["status"], low, "unknown root weight", missing, flush=True)
    assert hashlib.sha256((BASE/"RESULTS.json").read_bytes()).hexdigest() == output["original_S5_results_sha256"]
    return output


if __name__ == "__main__":
    run()
