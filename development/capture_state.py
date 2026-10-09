"""Capture existing pure provenance interfaces in a separate process."""
import argparse
from importlib.metadata import version
import json
from pathlib import Path
import platform
import sys

from integrity import source_identity


def capture(root):
    root = Path(root).resolve()
    sys.path.insert(0, str(root / "provenance"))
    import expected_red as expected
    import bankgate as bank
    import auditor  # Selected importability control; importing arbitrary tools can have side effects.
    import resident
    state = source_identity(root)
    state["declarations"] = {node: {"cases": entry["cases"],
        "live_cases": sorted(entry["enumerate"]())} for node, entry in sorted(expected.DECLARED.items())}
    state["open_passes"] = expected.open_passes()
    baseline, label = bank._resolve_baseline()
    if baseline is None:
        raise ValueError("Accepted bank baseline is missing")
    working, sources = bank._load_register()
    state["bank_inventory"] = {"baseline_label": label,
        "report": bank.bank_gate(baseline, working, sources), "held_flags": bank._load_held()}
    state["tool_versions"] = {"python": platform.python_version(),
        "packages": {name: version(name) for name in ("pytest", "iniconfig", "packaging", "pluggy", "Pygments")}}
    return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(capture(args.root), sort_keys=True))
