# CP-12 induced gauge dynamics

Virtual matter hopping generates an internal Wilson-loop coupling; no
magnetic term is added to the frozen Hamiltonian. The exact frozen-parameter
benchmark is the complete three-vertex graph. A static/heavy-matter calculation
also derives a general simple-loop coefficient and degree-homogeneity energy
term. These are a physical mechanism on a graph, not spacetime emergence.

- `DERIVATION.md`: equations, exact rational triangle coefficient, physical
  interpretation, unbounded-electric control, and limits.
- `run.py`: executable checks. It imports the two pre-existing, independently
  implemented physical transition/basis modules from
  `../../cp12_ground_connectivity/`.
- `RESULTS.json`: saved outputs, including all cutoff and parameter controls.
- `RUN_LOG.txt`: latest full execution output.
- `SOURCE_AUDIT.json`: exact file checksums and runtime versions.

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/constructive_sprint_20261015/induced_gauge/run.py
```

Dependencies are NumPy and SciPy; the tested versions are recorded in the audit.
The calculation takes a few seconds on the sprint host, performs exact rational
path arithmetic, independently checks an unpruned unregulated sixth Hamiltonian
power, and diagonalizes the full physical triangle basis through electric
cutoff 16 (638 states).

Exact identities and path coefficients are distinguished from float64
eigensolvers and analytically derived numerical error bounds. No interval
arithmetic certificate is claimed. The standalone hopping-family controls
vary `t` explicitly to test a perturbative coefficient; the main frozen-law
calculation keeps `m=t=1`, `epsilon_E=.25`, `kappa=.5`, `Gamma=.1`.

The positive Wilson coefficient favors a holonomy angle, which is conjugate to
electric flux. It does not select a definite persistent electric-flux network.
The prior all-time low-matter circulation obstruction is unchanged.
