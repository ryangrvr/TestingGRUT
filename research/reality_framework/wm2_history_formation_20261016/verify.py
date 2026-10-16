"""Verify frozen protocols and compute the complete, threshold-paid ledger."""
from math import exp, factorial, sqrt
from pathlib import Path
import hashlib
import json

import bounds as B
import campaign as C

BASE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_case(r, strang, projection):
    assert r["N"] == sum(r["nonroot_type_counts"])+1
    assert r["receiver_modes"] == 64*r["N"]**2
    assert r["observation_time"] == 4*r["N"]**2
    assert r["exact_history_retained"] and not r["extra_photon_cutoff"]
    assert not r["roundoff_certified"]
    assert r["maximum_active_rank"] <= 24
    assert r["maximum_coefficient_dimension"] <= 2_000_000
    assert 0 < r["completed_steps"] <= r["steps"]
    assert abs(r["dt"]*r["steps"]-r["observation_time"]) < 1e-10
    assert (r["status"] == "COMPLETE") == (r["completed_steps"] == r["steps"])
    assert abs(sum(r["populations"])-1) < 2e-9
    assert len(r["populations"]) == r["N"]
    assert r["error_components"]["strang"] <= strang+1e-10
    assert r["error_components"]["orbital_projection"] <= projection+1e-8
    assert abs(sum(r["error_components"].values())-r["analytic_state_error_account"]) < 1e-10
    if r["status"] != "COMPLETE":
        assert r["analytic_probability_lower"] is None
        assert r["incomplete_population_notice"]


def families(raw, case_audits, supplementary):
    output = []
    for family in raw["families"]:
        expected = C.components(family)
        assert len(expected) == len(family["components"])
        lower = full = missing = 0.0
        gates = []
        margins = []
        for component, target in zip(family["components"], expected):
            for name in ["root_type", "weight", "counts", "values", "g"]:
                assert component[name] == target[name]
            index = component["P1_index"] if supplementary else component["case_index"]
            result, audited = raw["attempts"][index], case_audits[index]
            actual = C.key(result["nonroot_type_counts"], result["detuning_values"], result["g"])
            direct = C.key(component["counts"], component["values"], component["g"])
            mirror = C.mirror_key(component["counts"], component["values"], component["g"])
            assert actual in (direct, mirror)
            assert result["N"] == family["N"]
            if audited["status"] != "COMPLETE":
                missing += component["weight"]
            else:
                lower += component["weight"]*audited["formation_lower"]
                full += component["weight"]*audited["full_assembly_lower"]
                p = min(1.0, max(0.0, result["macro_probability_estimate"]))
                allowance = sqrt(p)-sqrt(.5)
                if supplementary:
                    allowance = max(allowance, sqrt(.5)-sqrt(1-p))
                margins.append(allowance-audited["corrected_state_error_account"])
            gates.append(audited["corrected_state_error_account"] <= (.60 if supplementary else .20))
        assert abs(sum(c["weight"] for c in family["components"])-1) < 1e-12
        assert abs(missing-family["unresolved_root_weight"]) < 1e-12
        assert (missing == 0) == (family["status"] == "COMPLETE")
        n, time = family["N"], 4*family["N"]**2
        row = dict(name=family["name"], role=family["role"], N=n,
                   bar_delta=(sum(family["values"])/len(family["values"]) if len(family["values"]) == 2 else None),
                   status=family["status"], original_root_weights_retained=True,
                   unresolved_root_weight=missing, formation_lower=lower,
                   full_assembly_lower=full,
                   gate_pass=bool(missing == 0 and all(gates) and lower >= .5))
        if missing == 0 and margins and min(margins) > 0:
            row["finite_comparison_box"] = dict(per_site_storage_radius=min(margins)/(4*time*(n-1)),
                                               common_g_radius=min(margins)/(4*time),
                                               intersect_with_frozen_class=True,
                                               all_size_nonshrinking_claim=False)
        probe_time = .005
        r = max(1, (n+1)//2-2)
        tau = sum((2.1*probe_time)**k/factorial(k) for k in range(r, 50))
        d = .25/2.1*tau
        kappa = .25*probe_time-(exp(4.2*probe_time)-1-4.2*probe_time)
        row["joined_protocol"] = dict(source_variance=.25, probe_time=probe_time,
            no_receiver_reset=True, all_macro_records_accepted=True, r_min=r,
            accept_and_consistent_lower=lower*(1-min(1,d)**2),
            near_response_per_original_trial_lower=lower*kappa)
        output.append(row)
    return output


def audit():
    s5 = json.loads((BASE/"RESULTS.json").read_text())
    p1 = json.loads((BASE/"PROBABILITY_RESULTS.json").read_text())
    assert len(s5["families"]) == 8 and len(s5["attempts"]) == 13
    assert len(p1["families"]) == 6 and len(p1["attempts"]) == 10
    assert p1["original_S5_results_sha256"] == sha(BASE/"RESULTS.json")
    assert not p1["target_law_changed"] and not p1["seed_changed"]
    assert not p1["observation_time_changed"] and not p1["disorder_test_list_changed"]
    assert s5["plan"]["families"] == C.families()
    for r in s5["attempts"]:
        check_case(r, .04, .15)
    for r in p1["attempts"]:
        check_case(r, .20, .35)
    for r in s5["controls"]:
        check_case(r, .04, .15)
        assert r["quadratic_reference"]["state_error"] <= B.state_error(r)+1e-8
    wave = s5["checks"]["complete_wavefunction_checks"]
    assert len(wave) == 6
    assert all(x["complete_wavefunction_error"] <= x["analytic_state_error_account"]+2e-10 for x in wave)
    extra = s5["extra_checks"]
    assert len(extra["sign_symmetry"]) == 6
    assert max(x["complete_basis_population_difference"] for x in extra["sign_symmetry"]) < 1e-11
    assert extra["probability_bound"]["analytic_two_dimensional_controls"] == 10201
    assert extra["probability_bound"]["maximum_violation"] < 2e-12
    ref = extra["physical_receiver_reference"]
    assert ref["complete_physical_matrix_dimension"] == 333505
    assert ref["complete_wavefunction_error"] <= ref["analytic_state_error_account"]+2e-9
    original = [B.audited_case(r, False, .20) for r in s5["attempts"]]
    probability = [B.audited_case(r, True, .60) for r in p1["attempts"]]
    f5, fp = families(s5, original, False), families(p1, probability, True)
    return dict(S5_charter_commit="64378573e8251224dd563928667e961a0effb3da",
                P1_charter_commit="364781b841121e3169041b121e661048fd793524",
                original_S5_results_sha256=sha(BASE/"RESULTS.json"),
                original_results_preserved=True,
                threshold_correction="Universal g*T*1e-13 added where early raw records omitted it",
                largest_contact_threshold_upper=max(r["g"]*r["observation_time"]*1e-13 for r in s5["attempts"]),
                S5_case_audits=original, P1_case_audits=probability,
                S5_families=f5, P1_families=fp,
                all_size_independent_disorder_formation_proved=False,
                externally_validated=False, outward_rounded=False)


def summary(a):
    s5 = json.loads((BASE/"RESULTS.json").read_text())
    p1 = json.loads((BASE/"PROBABILITY_RESULTS.json").read_text())
    return dict(numerical_consistency="PASS",
        physical_status="Controlled finite formation; all-size independent-disorder gate OPEN",
        frozen_primary_families=6,
        original_S5_completed_primary=sum(x["gate_pass"] for x in a["S5_families"] if x["role"] == "PRIMARY"),
        P1_completed_primary=sum(x["gate_pass"] for x in a["P1_families"]),
        S5_unique_components=13, P1_unique_components=10,
        small_complete_wavefunctions=6, physical_complete_wavefunctions=1,
        physical_reference_dimension=333505,
        physical_reference_state_error=s5["extra_checks"]["physical_receiver_reference"]["complete_wavefunction_error"],
        physical_reference_error_upper=s5["extra_checks"]["physical_receiver_reference"]["analytic_state_error_account"],
        uniform_exact_controls=len(s5["controls"]),
        sign_identity_controls=6, probability_inequality_controls=10201,
        maximum_active_rank=max(r["maximum_active_rank"] for r in s5["attempts"]+p1["attempts"]),
        maximum_coefficient_dimension=max(r["maximum_coefficient_dimension"] for r in s5["attempts"]+p1["attempts"]),
        contact_threshold_paid=True, old_results_unchanged=True,
        external_independent_validation=False, outward_rounded_arithmetic=False,
        all_size_formation_proved=False)


if __name__ == "__main__":
    result = audit()
    print(json.dumps(summary(result), indent=2))
