"""Optional publication plot, using only archived registered results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parent


def main():
    data = json.loads((ROOT / "RESULTS.json").read_text())
    cases = data["cases"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "GRUT-D2-01"})
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), layout="constrained")
    x = np.arange(len(cases))
    initial = [2**(-c["m"]) for c in cases]
    final = [c["time_rows"][-1]["order_probability"] for c in cases]
    echo = [c["bath_phase_echo"]["echo_final_order_probability"] for c in cases]
    for shift, values, color, label in [(-.25, initial, "#82909d", "Initially"),
                                        (0, final, "#117c83", "Autonomous, at 2τ"),
                                        (.25, echo, "#c57b36", "Bath phase kick at τ")]:
        axes[0].bar(x + shift, values, width=.24, color=color, label=label)
    axes[0].set_xticks(x, [str(c["n"]) for c in cases])
    axes[0].set_xlabel("Physical system spins n")
    axes[0].set_ylabel("Probability all system spins are aligned")
    axes[0].set_ylim(0, 1.05)
    axes[0].yaxis.set_major_formatter(PercentFormatter(1))
    axes[0].set_title("Order forms; a bath phase kick returns disorder", loc="left", fontsize=12)
    axes[0].legend(frameon=False, fontsize=9, loc="upper center",
                   bbox_to_anchor=(0.5, -0.19), ncols=3)
    for c in cases:
        rows = c["time_rows"]
        axes[1].semilogy([r["t"]/c["window_start"] for r in rows],
                         [1-r["order_probability"] for r in rows], label=f'n={c["n"]}')
    axes[1].axhline(1-cases[0]["certificate"]["order_lower_display"],
                    color="#39444e", linestyle="--", label="Uniform analytic ceiling")
    axes[1].set_xlabel("Time / τ  (τ = 8 Σ i²,  g = 1)")
    axes[1].set_ylabel("Probability of at least one system domain wall")
    axes[1].set_xlim(1, 2)
    axes[1].set_title("Finite-window bounds and exact mode sums", loc="left", fontsize=12)
    axes[1].legend(frameon=False, ncols=2, fontsize=9, loc="upper right")
    fig.suptitle("GRUT-D2 · supplied local alignment law and prepared internal reservoir", fontsize=13)
    fig.savefig(ROOT / "ordering.png", dpi=180)
    fig.savefig(ROOT / "ordering.svg", metadata={"Date": None})
    svg = ROOT / "ordering.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")


if __name__ == "__main__":
    main()
