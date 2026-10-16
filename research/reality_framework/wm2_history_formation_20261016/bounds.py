"""Interpret saved trajectories with the complete analytic error account."""
from math import sqrt

CONTACT_THRESHOLD = 1e-13


def state_error(result):
    # Early saved trajectories did not print this tiny conservative term.
    # Preserve them verbatim and include the term in the audited ledger.
    extra = (0.0 if "contact_threshold" in result["error_components"] else
             result["g"]*result["observation_time"]*CONTACT_THRESHOLD)
    return result["analytic_state_error_account"]+extra


def probability_lower(p, error, complementary=False):
    p = min(1.0, max(0.0, p))
    reverse = max(0.0, sqrt(p)-error)**2
    other = max(0.0, 1-(sqrt(1-p)+error)**2) if complementary else 0.0
    return max(reverse, other)


def audited_case(result, complementary=False, envelope=.20):
    complete = result["status"] == "COMPLETE"
    error = state_error(result)
    low = probability_lower(result["macro_probability_estimate"], error, complementary) if complete else None
    full = probability_lower(result["full_assembly_probability_estimate"], error, complementary) if complete else None
    return dict(status=result["status"], corrected_state_error_account=error,
                added_contact_threshold_upper=error-result["analytic_state_error_account"],
                formation_lower=low, full_assembly_lower=full,
                gate_pass=bool(complete and error <= envelope and low >= .5))
