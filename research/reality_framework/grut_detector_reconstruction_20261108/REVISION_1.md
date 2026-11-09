# B2-01r1: startup correction before production

The first source freeze was commit 9f39ba5e721853de46b3f19017cb7320f1c05f79. Launching run.py exited during imports, before output directory creation, prechecks, or any production point:

    ImportError: cannot import name 'massive' from 'analyze'
    (.../grut_detector_matching_20261107/analyze.py)

Cause: model.py inserted the inherited module directory at sys.path[0], which shadowed this study's analyze.py. Correction: insert it at sys.path[1], retaining this study's script directory first. Run an import smoke check in addition to all existing prechecks, then commit this revision before production.

No B2 production probability was calculated or examined by this failed startup. Candidate reconstruction, timing list, thresholds, solver settings and resource caps are unchanged. The protocol identifier is B2-01r1 to preserve this correction explicitly. The next clean launch is the first numerical production attempt; this document retains the earlier startup failure.
