# Reproduce the compact packet

Requires the parent bb211d1 tree, Python, NumPy, SciPy and Matplotlib.
No new installation was used. From the repository root:

```bash
python3 research/reality_framework/grut_conversion_streaming_20261023/code/streaming.py
python3 research/reality_framework/grut_conversion_streaming_20261023/code/first_conversion.py
```

The nine production runs use three grids/time steps/source-bin widths for
each of three imposed epochs. They take several minutes on the builder
runtime. Run the no-reaction diagnostic before analyzing:

```bash
PYTHONPATH=research/reality_framework/grut_conversion_streaming_20261023/code python3 - <<'PY'
import json
from streaming import ROOT,run
r=run('no_reactions_control',1e-6,800,250,.001,.005,reactions=False)
(ROOT/'NO_REACTION_CONTROL.json').write_text(json.dumps(r,indent=2)+'\n')
PY
python3 research/reality_framework/grut_conversion_streaming_20261023/code/analyze.py
python3 research/reality_framework/grut_conversion_streaming_20261023/code/verify.py
python3 research/reality_framework/grut_conversion_streaming_20261023/code/publish_report.py
python3 research/reality_framework/grut_conversion_streaming_20261023/code/manifest.py --check
```

The no-reaction diagnostic fails as a path estimator and is deliberately
fenced out. The released horizon calculation uses the separately tested
receiver, quadrature and light-cone prefix bound. A true independent
replication should challenge those derivations and recompute the source.

All production histories and time/momentum receiver tables are included
as compact JSON; there is no bulk raw-field archive in this packet.
Numeric source columns are named by `escape_source_columns`. Parameter,
time, rate and profile approximations are explicit in MATHEMATICS.md.
