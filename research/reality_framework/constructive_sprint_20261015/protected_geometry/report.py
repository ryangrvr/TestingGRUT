"""Render the same-state result summary and a compact scientific figure."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

P = Path(__file__).resolve().parent


def main():
    data = json.loads((P / "RESULTS.json").read_text())
    rows = [r for r in data["boxes"][-1]["results"] if r["time"] == 1.]
    lines = ["# Supplied geometry and matter response in the same evolving state", "",
             "The frozen CP-12 Hamiltonian preserves a supplied flux-carrying square "
             "while supporting nonzero, gauge-invariant matter-density response on "
             "that same evolving preparation. This is a conditional mechanism; the "
             "graph and initial flux are explicitly supplied.", "",
             "All 15 labelled four-edge graphs on four vertices are included. "
             "After unforced wait 0.5, a phase pulse exp(−iθn₀) is followed by lag "
             "0.5. Receiver 2 is at initial-host distance two. No loop filter is used.", "",
             "| Flux L | Host probability | Distance-two susceptibility | Truncation error | Response from escaped host |",
             "|---:|---:|---:|---:|---:|"]
    for r in rows:
        d = r["detector_responses"][2]
        lines.append(f'| {r["flux"]} | {r["host"]["value"]:.15g} | {d["response"]:.12g} | '
                     f'{d["error_bound"]:.3g} | {d["escaped_host_response"]:.3g} |')
    lines += ["", "The response is overwhelmingly carried by the original-host component. "
              "Joint host–density responses are saved separately, so rare changed-graph "
              "outcomes are not being mistaken for transmission through the retained host.", "",
              "Cutoffs |E|≤6 and |E|≤10 agree within the analytically bounded "
              "outgoing error and a reported floating-point tolerance. Every retained "
              "Hamiltonian column was independently compared with the unregulated "
              "integer transition implementation; every omitted transition was retained "
              "in the leakage operator. Bounds are ordinary floating-point evaluations, "
              "not formal interval-arithmetic certificates. Reported intervals add "
              "a separate numerical allowance of 5×10⁻¹²; this exceeds the observed "
              "cutoff discrepancies but is not a certified roundoff bound. The table "
              "separately lists the analytic electric-truncation estimate.", "",
              "Two exact results support the mechanism. A supplied host with "
              "initial edge fields (±L,∓L) obeys", "",
              r"\[1-p_g(T)\le Q\left[\frac{(\sqrt2|t|T)^{2L}}{(2L)!}\right]^2.\]", "",
              "Without any flux assumption, a fixed host also obeys "
              r"$\|U(T)\psi-e^{-iP_gHP_gT}\psi\|\le "
              r"[\Gamma\sqrt{Q(M-Q)}/M]T$. For sparse hosts this error is "
              "O(T/√N). These are supplied-host preservation results; the original "
              "neutral coherent input occupies many hosts and is not covered by "
              "the fixed-host root-mean-square estimate.", "",
              "Exact signed, unregulated Taylor calculations give graph-loss onset "
              "T⁸/259200 for L=1 and "
              "(216017/12906142470144000000)T²⁴ for L=2. Higher-flux onset "
              "certificates and the complete Gauss target catalog are in EXACT_ONSET.json.", "",
              "The test does not establish a macroscopic causal cone, autonomous "
              "geometric formation, preferred dimension, or gravity. The four-vertex "
              "host has diameter two. Initial electric energy and localized matter "
              "are physical preparation resources and are not matched to the original "
              "neutral ray.", "",
              "Reproduce:", "", "```bash",
              "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py",
              "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python exact_onset.py",
              "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python report.py", "```", ""]
    if (P / "NEUTRAL_RECORD.json").exists():
        n = json.loads((P / "NEUTRAL_RECORD.json").read_text())
        response = n["density_response"]
        lines += ["## Neutral-origin recorded-history extension", "",
                  "A separate executed prototype begins from the original neutral ray, "
                  "evolves unforced to time5, then applies a complete graph-and-circulation "
                  "diagnostic record before the same gauge-invariant source probe. "
                  "Charged full-circulation records have probability "
                  f'{n["event_probability_at_wait"]:.12g}. At lag0.5, their conditional '
                  f'host survival is {n["host_survival"]["conditional_probability"]:.12g}, '
                  f'circulation survival {n["circulation_survival"]["conditional_probability"]:.12g}, '
                  "and opposite-on-recorded-ring density susceptibility "
                  f'{response["conditional_on_initial_record_event"]:.12g} '
                  f'with interval {response["conditional_interval"]}.', "",
                  "The complete external record costs system energy "
                  f'+{n["complete_record_energy_cost"]:.12g} and measured-system entropy '
                  f'{n["complete_record_entropy_bits"]:.12g} bits. All outcomes are retained; '
                  "the rare charged branch is analyzed conditionally. This joins formation "
                  "and propagation in one finite recorded history while explicitly costing "
                  "the diagnostic intervention. It is not an autonomous geometric phase.", "",
                  "See NEUTRAL_RECORD.md for exact symmetry handling, normalization, "
                  "bounded energy-cost error control and numerical qualifications.", ""]
    (P / "REPORT.md").write_text("\n".join(lines))
    fig, ax = plt.subplots(figsize=(7.1, 4.3))
    x = np.arange(3)
    for off, r, color in zip([-.14,.14], rows, ["#245b91", "#cb6a32"]):
        v = [abs(r["detector_responses"][j]["response"]) for j in [1,2,3]]
        ax.bar(x + off, v, width=.28, label=f'Flux L={r["flux"]}', color=color)
    ax.set_yscale("log")
    ax.set_ylabel("Absolute density susceptibility")
    ax.set_xticks(x, ["Vertex 1\nDistance 1", "Vertex 2\nDistance 2", "Vertex 3\nDistance 1"])
    ax.set_title("Matter response on the same protected, supplied square")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=.2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(.5,.015,"N=4, Q=4; wait=0.5, lag=0.5. This is not an emergent causal-cone test.",
             ha="center",fontsize=9)
    fig.tight_layout(rect=[0,.045,1,1])
    fig.savefig(P / "same_state_response.png", dpi=180)
    fig.savefig(P / "same_state_response.svg")


if __name__ == "__main__":
    main()
