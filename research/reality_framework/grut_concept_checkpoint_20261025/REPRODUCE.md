# Reproduction and provenance

From the repository root:

```bash
python research/reality_framework/grut_concept_checkpoint_20261025/code/check_sources.py
```

Requires Python, NumPy, SciPy and SymPy. Runtime is a few seconds. It writes
CHECKS.json only in this packet. The checks derive the exact current sources,
compare them to full complex-field forces and independent angular finite
differences, compare the quadratic matching formula to a real-matrix
solve and a full nonlinear instantaneous minimum, and evaluate a conditional
leading-EFT radiation-envelope yield bound. No production release,
coupling scan or cosmological particle history is run.

CHECKS.json distinguishes the four mathematical control groups from the
physical/cosmological questions they do not answer. These are author checks;
no independent review of this packet is represented as completed.

SOURCES.json pins the internal files used from the immutable parent
b448f25821dd7a363d9611f5757fd53460e053d0 with Git blob and SHA-256 hashes and
their read scopes. The latest owner-supplied reviewer note is the prompt for
this checkpoint, not a completed critic report or a primary literature
verification. No new external literature claim is needed for the algebra.

MANIFEST.json hashes the packet's final files except the manifest itself.
