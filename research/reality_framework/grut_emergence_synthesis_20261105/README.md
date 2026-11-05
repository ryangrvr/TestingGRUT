# GRUT emergence synthesis

Start with [GRUT_WORKING_THEORY.md](GRUT_WORKING_THEORY.md), then use the
[map](EMERGENCE_MAP.svg) and [claim ledger](CLAIMS.md) to inspect its boundaries.
[SOURCES.md](SOURCES.md) gives pinned project evidence and external references.
[AUDIT.md](AUDIT.md) states what was actually checked and what was not.
[HANDOFF.md](HANDOFF.md) gives a bounded independent-review brief.

This is a reviewable, additive theory account. It is not an independently reviewed
new fundamental law or a cosmologically certified parameter point. The latest
strict M2-X result is zero passes among 135 preparation keys. No old register or
physical result is changed by this synthesis.

The finite demonstrations have analytic answers known in advance. Their 48 checks
exercise whole/part correlations, records, the Gibbs-reservoir entropy identity
and a supplied information-using controller. They do not test consciousness.

## Reproduce the packet

Python dependencies used: numpy 2.3.5, scipy 1.17.0, matplotlib (version recorded
by the environment; the image rendering is separate from the numerical checks).
From this directory:

    python -B code/witnesses.py
    python -B code/build_metadata.py
    python -B code/audit.py --source-root /path/to/pinned/TestingGRUT

The source root must contain the exact files recorded in SOURCE_MANIFEST.json at
8ab5829a69f51cf462c56b89252d78fa21760300. The audit checks their Git blob IDs and
SHA-256 hashes. It reaggregates the archived root verdict; it does not rerun the
scan engine, physical lattice or fermion solver.

Optional reader build, with Pandoc installed:

    python -B code/build_reader.py

The HTML reader contains the manuscript, diagram, claims, sources and audit.
MathML is used so no external math-rendering service is required.

## Scope and recovery

[SCOPE.md](SCOPE.md) was published before executing the reconstructed finite
checks. The earlier unsaved draft was unavailable after a workspace reset. This
packet is a reconstruction from pinned evidence, not a recovered copy of that
draft. No old unpublished local commit is counted as evidence.

The new checks and interpretive review are by the same author. No Claude or
Argon review is claimed. Their historical work retains the attribution in its
own source files.
