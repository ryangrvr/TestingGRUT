"""Inventory the pinned S1 source tree and screen every Markdown document."""
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path, PurePosixPath

BASE = "d3dd471a0ad69a95cd6f0dc1eba23fc0695a39cd"
THEMES = {
    "universal_closure": r"closed.whole|universal closure|no.boundary|no exterior|system.complement",
    "quantum_kinematics": r"quantum state|unitarity|noncommutativity|Born rule|Hilbert space",
    "internal_time": r"internal (?:time|clock)|Page[–—-]Wootters|time.reparametri[sz]ation|relational time",
    "locality_and_causality": r"Lieb[–—-]Robinson|emergent locality|causal cone|relational locality|CP-12",
    "reconstruction_and_geometry": r"inverse (?:problem|spectrum)|isospectral|reconstruction|conformal (?:structure|geometry|covariance)",
    "matter_interactions": r"Standard Model|Yukawa|Candidate 1|scalar sector|matter sector",
    "gravity": r"gravitational constraint|Einstein[–—-]Hilbert|gravitational response|Ward identit|spin.2",
    "electric_charge": r"hypercharge|neutron charge|charge (?:dequantization|quantization)|gauge anomaly",
    "initial_conditions": r"initial condition|initial.state|primordial spectrum|inflationary (?:state|preparation)",
    "dark_matter": r"dark matter|Q.ball|relic abundance|sterile fermion|H55",
    "dark_energy": r"dark energy|vacuum energy|cosmological constant|phantom (?:cross|divid)",
    "thermodynamics": r"thermodynamic|entropy production|detailed balance|Landauer|finite.fuel",
    "physical_information": r"quantum information|mutual information|information capacity|state discrimination",
    "memory_and_identity": r"historical (?:record|identity|restor)|record.writ|memory.(?:enabled|dependent)|D7",
    "developmental_organization": r"developmental|self.organi[sz]|persistent structure",
    "predictive_adaptation": r"predictive (?:adaptation|record|memory|information)|adaptive development|learning",
    "physical_observers": r"internal observer|observer.algebra|measurement instrument|relational observab",
    "consciousness": r"consciousness|subjective awareness|experience of time",
}
DOSSIER = re.compile(
    r"(?:^|[_-])(?:verdict|ruling|adjudication|correction|synthesis|closeout|"
    r"negative_results|open_claims|no_go_ledger|claim_ledger|claim_boundary)(?:[_-]|\.)",
    re.I,
)


def git(*args, input_data=None):
    return subprocess.check_output(["git", *args], input=input_data)


def read_markdown(paths):
    object_ids = sorted(set(paths.values()))
    if not object_ids:
        raise ValueError("The pinned tree had no Markdown corpus")
    raw = git("cat-file", "--batch", input_data=("\n".join(object_ids) + "\n").encode())
    position = 0
    documents = {}
    for object_id in object_ids:
        newline = raw.index(b"\n", position)
        header = raw[position:newline].decode().split()
        if header[:2] != [object_id, "blob"] or len(header) != 3:
            raise ValueError(f"Missing or invalid Markdown blob {object_id}")
        last = newline + 1 + int(header[2])
        documents[object_id] = raw[newline + 1:last].decode("utf-8", errors="replace")
        position = last + 1
    if position != len(raw):
        raise ValueError("Unexpected trailing bytes in Git batch")
    return documents


def build():
    raw_tree = git("ls-tree", "-r", "-z", "--full-tree", BASE)
    paths = {}
    for entry in raw_tree.split(b"\0"):
        if entry:
            metadata, name = entry.split(b"\t", 1)
            paths[name.decode()] = metadata.split()[2].decode()
    extensions = Counter()
    roots = Counter()
    packets = Counter()
    markdown_paths = {}
    dossiers = []
    for name, oid in paths.items():
        path = PurePosixPath(name)
        extension = path.suffix.lower() or "[no suffix]"
        extensions[extension] += 1
        roots[path.parts[0] if len(path.parts) > 1 else "[root]"] += 1
        if path.parts[:2] == ("research", "reality_framework") and len(path.parts) > 3:
            packets[path.parts[2]] += 1
        if extension == ".md":
            markdown_paths[name] = oid
        if extension in (".md", ".json") and DOSSIER.search(path.name):
            dossiers.append(name)

    documents = read_markdown(markdown_paths)
    themes = {name: re.compile(pattern, re.I) for name, pattern in THEMES.items()}
    matches = {
        oid: [name for name, pattern in themes.items() if pattern.search(doc)]
        for oid, doc in documents.items()
    }
    theme_counts = Counter()
    for oid in markdown_paths.values():
        theme_counts.update(matches[oid])
    return {
        "source_commit": BASE,
        "source_tree": git("rev-parse", f"{BASE}^{{tree}}").decode().strip(),
        "tree_manifest_sha256": hashlib.sha256(raw_tree).hexdigest(),
        "tracked_paths": len(paths),
        "path_counts_by_extension": dict(extensions),
        "path_counts_by_root": dict(roots),
        "reality_framework_packet_path_counts": dict(packets),
        "markdown_paths_screened": len(markdown_paths),
        "distinct_markdown_blobs_screened": len(documents),
        "markdown_paths_matching_a_theme": sum(bool(matches[oid]) for oid in markdown_paths.values()),
        "theme_patterns_case_insensitive": THEMES,
        "theme_match_counts": dict(theme_counts),
        "decision_filename_pattern_case_insensitive": DOSSIER.pattern,
        "decision_dossier_candidates": sorted(dossiers),
    }


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, sort_keys=True))
