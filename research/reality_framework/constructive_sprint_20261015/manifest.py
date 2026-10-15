"""Hash the additive research packet and its unchanged source dependencies."""
from pathlib import Path
import hashlib
import json

P = Path(__file__).resolve().parent


def packet_files():
    for path in sorted(P.rglob("*")):
        if not path.is_file() or any(x in path.parts for x in ["__pycache__", "unpublished_local_backup"]):
            continue
        if str(path.relative_to(P)) not in {"MANIFEST.json", "VERIFICATION.json"}:
            yield path


def main():
    files = [{"path": str(f.relative_to(P)), "bytes": f.stat().st_size,
              "sha256": hashlib.sha256(f.read_bytes()).hexdigest()} for f in packet_files()]
    dependencies = []
    for relative in ["cp12_ground_connectivity/frozen_model.py",
                     "cp12_ground_connectivity/exact_destination.py",
                     "cp12_distance_history/evolution.py", "cp12_unmasked_formation/run.py"]:
        path = P.parent/relative
        dependencies.append(dict(path=relative, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    result = dict(baseline="a2f447dc37ea0a71f62c6e990a4630d727ebd85e",
                  parent="1801afee27096c03c2a3b64cc811ac0be0d2a4a7",
                  parent_tree="6e159afb2d1f08d0a9ffebccf28cb7db04585fdc",
                  publication_scope="Only additions under research/reality_framework/constructive_sprint_20261015; no prior or frozen file changes.",
                  exclusions="This manifest and root VERIFICATION.json avoid self-referential hashes; Python bytecode and unpublished local backup are not published.",
                  files=files, unchanged_dependencies=dependencies)
    (P/"MANIFEST.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(files=len(files), bytes=sum(x["bytes"] for x in files))))


if __name__ == "__main__":
    main()
