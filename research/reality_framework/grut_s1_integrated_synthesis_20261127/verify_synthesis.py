"""Check S1 document coverage, source pins, and exact archive-index replay."""
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = "d3dd471a0ad69a95cd6f0dc1eba23fc0695a39cd"
DELIVERABLES = {
    "GRUT_WORKING_THEORY.md", "MATHEMATICAL_FOUNDATIONS.md",
    "COMPLETE_EVIDENCE_MAP.md", "FUNDAMENTAL_ASSUMPTIONS.md",
    "DERIVED_RESULTS.md", "DEVELOPMENTAL_EFFECTIVE_THEORY.md",
    "FOUNDATIONAL_OBSTRUCTIONS.md", "CANDIDATE_1_STATUS.md",
    "NOVELTY_INVENTORY.md", "NEGATIVE_RESULTS_INDEX.md", "OPEN_PROBLEMS.md",
    "REPRODUCTION_INDEX.md", "EXTERNAL_RAI_REVIEW_PACKET.md", "FINAL_PROGRAM_VERDICT.md",
    "INTERNAL_CRITICAL_REVIEW.md",
}
LOCAL_LINK = re.compile(r"\]\(([^()\s]+)\)")
SEPARATOR = re.compile(r"\|\s*:?-{3,}:?\s*(?:\|\s*:?-{3,}:?\s*)+\|")
TABLE_DELIMITER = re.compile(r"(?<!\\)\|")


def git(*args):
    return subprocess.check_output(["git", *args], cwd=HERE, text=True).strip()


def main():
    missing = sorted(name for name in DELIVERABLES if not (HERE / name).is_file())
    if missing:
        raise SystemExit(f"Missing deliverables: {', '.join(missing)}")
    root = Path(git("rev-parse", "--show-toplevel"))
    old_paths = set(git("ls-tree", "-r", "--full-tree", "--name-only", BASE).splitlines())
    relative_links = 0
    table_count = 0
    for doc in sorted(HERE.glob("*.md")):
        rows = doc.read_text().splitlines()
        width = None
        for line_no, row in enumerate(rows, 1):
            if SEPARATOR.fullmatch(row):
                width = len(TABLE_DELIMITER.findall(row))
                if len(TABLE_DELIMITER.findall(rows[line_no - 2])) != width:
                    raise SystemExit(f"{doc.name}:{line_no} table header width differs")
                table_count += 1
            elif width and row.startswith("|"):
                if len(TABLE_DELIMITER.findall(row)) != width:
                    raise SystemExit(f"{doc.name}:{line_no} table row width differs")
            else:
                width = None
            for link in LOCAL_LINK.findall(row):
                if "://" in link or link.startswith("#"):
                    continue
                path = (HERE / link.split("#", 1)[0]).resolve()
                if not path.is_relative_to(root):
                    raise SystemExit(f"{doc.name}:{line_no} path outside repository")
                if not path.is_file() and str(path.relative_to(root)) not in old_paths:
                    raise SystemExit(f"{doc.name}:{line_no} unresolved link {link}")
                relative_links += 1

    archive = (HERE / "ARCHIVE_INDEX.json").read_bytes()
    if json.loads(archive)["source_commit"] != BASE:
        raise SystemExit("Archive index names the wrong source commit")
    replay = subprocess.check_output([sys.executable, str(HERE / "build_archive_index.py")], cwd=HERE)
    if replay != archive:
        raise SystemExit("Pinned archive index did not replay byte-for-byte")
    print(f"{len(DELIVERABLES)} manuscripts and review files present.")
    print(f"{table_count} tables checked, {relative_links} relative links resolved.")
    print("Pinned archive index replayed byte-for-byte.")


if __name__ == "__main__":
    main()
