"""Scientific explanatory figure from the computed finite-model results."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import ROOT

r = json.loads((ROOT / "WHOLE_PARTS_RESULTS.json").read_text())
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.7), layout="constrained")
ts = np.arange(4)
axes[0].plot(ts, np.ones(4), "o-", color="#235789", label="Whole: global purity")
axes[0].plot(ts, [x["X_purity"] for x in r["records"]], "s--", color="#c45a22", label="Charged part: local purity")
axes[0].plot(ts, [x["record_entropy_bits"] for x in r["records"]], "D:", color="#32865d", label="Record: entropy (bits)")
axes[0].set(xlabel="Internal logical clock reading", ylabel="Purity or entropy (bits)",
            title="One stationary whole, changing internal records", ylim=(-.08, 1.15), xticks=ts)
fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=3,
           fontsize=8.5, frameon=False)
eta = np.linspace(0, 1, 201)
axes[1].plot(eta, 1-eta+eta*eta/2, color="#235789", lw=2)
axes[1].scatter([x["eta"] for x in r["clock_readouts"]],
                [x["probe_purity"] for x in r["clock_readouts"]], color="#c45a22", zorder=3)
axes[1].set(xlabel="Supplied logical readout noise η", ylabel="Conditional probe purity",
            title="Clock readout changes local coherence", ylim=(.45, 1.05))
axes[1].text(.05, .55, "Same pure whole throughout\nIdeal logical tick: purity = 1", fontsize=9)
fig.suptitle("An exact finite quantum example with global charge neutrality", fontsize=12)
fig.savefig(ROOT / "WHOLE_PARTS.png", dpi=170)
plt.close(fig)
print("Wrote WHOLE_PARTS.png from measured finite-model results.")
