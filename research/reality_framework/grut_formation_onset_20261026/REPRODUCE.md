# Reproduce the onset controls

Run at the TestingGRUT root with the existing Python/NumPy/SciPy/Matplotlib environment. No package installation was needed. The previous kinetic packet is an immutable dependency at34f23f8; its bootstrap supplies pinned e71730c inputs. Fetch e71730c and99aedde if absent.

```bash
git fetch origin e71730c24fe6a0b28ddd72cd30982c5b39c9aedf
git fetch origin 99aedde454c77f2c7dbabaab320a81e95a932fc7
python -B research/reality_framework/grut_formation_onset_20261026/code/circular.py
python -B research/reality_framework/grut_formation_onset_20261026/code/elliptic.py
python -B research/reality_framework/grut_formation_onset_20261026/code/inventory.py
python -B research/reality_framework/grut_formation_onset_20261026/code/roll.py
python -B research/reality_framework/grut_formation_onset_20261026/code/verify.py
python -B research/reality_framework/grut_formation_onset_20261026/code/plot.py
python -B research/reality_framework/grut_formation_onset_20261026/code/package.py --seal
```

The background uses DOP853 with rtol2e-12. The121-point monodromy scan has rtol3e-11, with continuous local peak refinement. Verification compares tighter DOP853 and Radau, a four-period transfer, independent radial quadratures, canonical identities, Hessian finite differences, a cubic/matrix spectrum comparison, and a linear-program energy bound. Two homogeneous seed controls stop at |phi|=.1v and retain linearity of the spatial perturbation. No nonlinear lattice or kinetic solver is called.

One exploratory endpoint evaluation produced a spurious growth of order sqrt(machine epsilon) by subtracting equal terms in D. Replacing it by the exact stationary identity D=2lambda_phi phi0² makes the endpoint exactly marginal; direct matrix eigenvalues retain only roundoff. The raw tiny charged k=0 Floquet exponent from a split unit Jordan block is kept in the data but is not interpreted as exponential physics. The phase/charge checks certify that symmetry sector instead.

ELLIPTIC_ORBIT.npz contains the cold orbit and neutral curvature; all other numerical summaries are JSON. There is no final particle/object population. package.py verifies actual input hashes, finite data, syntax, links and the manifest. Run without --seal to check the supplied manifest. No independent AI reviewer or subagents were used.
