# Execution record

24 November 2026. The D3 study starts from the D2 draft at
`675d991927edd520daa4b000680a15e8a5200fa7`. The construction and critical
review were assigned to separate AI agents within the same Work session.
The reduced solver and full tensor oracle were implemented separately.

- The analysis fixed the physical model, analytic gates, fault, controls
  and time rectangle before finite-system production. The design commit
  is `53fd7c9e13185e79f0983932f32807f543591e2a`. The two design files
  retain their original bytes. The construction proof was developed
  analytically, followed by the registered numerical checks. [FREEZE.json](FREEZE.json)
- Production evaluated 27 timings for the exact grid sector. A second
  propagator checked nine state vectors. The first replay exposed a type
  mismatch between a NumPy scalar and its JSON representation in the short
  control. Normalizing the scalar to a Python float resolved the mismatch.
  The same stored results then replayed with zero difference. [VERIFICATION.json](VERIFICATION.json)
- The critic checked six full physical active matrices across 54 timing
  pairs, including complete conditional energy laws in the physical basis.
  The short-chain pointer loss and the unflagged flip remain recorded.
  [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json)
- A copy of the packet was tested outside the checkout with the installed
  pinned dependencies. The reproduction includes the respective limits
  of same-code replay, the separate propagation algorithm and the small
  full-matrix implementation. External RAI and human review remain open.

Runtime was Python 3.12.14 with NumPy 2.3.5 and SciPy 1.17.0 on Linux.
OpenBLAS and OpenMP ran with one thread. Both frozen input hashes and
all report links were checked. The manifest covers the packet's
published files.
