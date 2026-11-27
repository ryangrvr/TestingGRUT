"""Replay and compare the two frozen-interface F4 implementations."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def replay(script, reference):
    raw = subprocess.check_output([sys.executable, str(HERE / script)])
    if raw != (HERE / reference).read_bytes():
        raise SystemExit(f"Byte-for-byte replay failed for {reference}")
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def leaves(value):
    if isinstance(value, dict):
        return sum(leaves(element) for element in value.values())
    if isinstance(value, list):
        return sum(leaves(element) for element in value)
    return 1


def main():
    production, prod_hash = replay("run_f4.py", "RESULTS.json")
    independent, ind_hash = replay("independent_check.py", "INDEPENDENT_RESULTS.json")
    fields = (
        "k", "weights", "charges", "lower_dimension", "dressed_dimension",
        "weinberg_allowed", "bY",
    )
    inputs = json.loads((HERE / "PARAMETERS.json").read_text())
    registered_k = inputs["fixtures_k"]
    if any(
        [row["k"] for row in result["fixtures"]] != registered_k
        for result in (production, independent)
    ):
        raise SystemExit("A replay omitted or reordered a frozen fixture")
    compared = 0
    for first, second in zip(production["fixtures"], independent["fixtures"]):
        for field in fields:
            if first[field] != second[field]:
                raise SystemExit(f"Independent disagreement at k={first['k']} for {field}")
            compared += leaves(first[field])
    if production["neutron_benchmark"]["absolute_comparison_scale"] != independent[
        "neutron_comparison_scale"
    ]:
        raise SystemExit("Independent disagreement on the empirical comparison scale")
    compared += 1
    print(f"{production['checks_passed']}/{production['checks_total']} production checks passed.")
    print(f"{compared} exact values agree between the two implementations.")
    print(f"Production replay  {prod_hash}")
    print(f"Independent replay {ind_hash}")


if __name__ == "__main__":
    main()
