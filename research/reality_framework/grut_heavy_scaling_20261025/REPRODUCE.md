# Reproduce the complementary controls

From the repository root:

```bash
python research/reality_framework/grut_heavy_scaling_20261025/code/check_scaling.py
```

Requires Python, NumPy, SciPy and Matplotlib. It writes RESULTS.json and
preparation_response.png in this packet. Four groups are checked:

1. Forty random complex field/velocity samples at four dilation factors:
   potential, full forces, current sources and number homogeneity.
2. Matching coefficients, entropy/yield and conditional lifetime/string
   scaling at five factors, with explicit fixed scales.
3. Ratio matching versus the dimensionful formula, independent quartic
   minima on either side of the positivity ceiling, and the conditional
   combined epsilon-bound formula.
4. Six driving frequencies, two initial preparations each, comparing exact
   filter, independent quadrature and direct real-oscillator integration.

The illustration uses no production source code from Claude. There is no
heavy three-field release, scan, quantum vector production or relic
calculation here. Same-author mathematical checks are not an independent
replication of the relayed spot point.

SOURCES.json pins the actual parent files used. RELAYED_POINT.json separates
approximate owner-relayed preliminary numbers from an immutable parameter
dictionary, which has not yet been supplied. MANIFEST.json hashes the final
packet excluding itself.
