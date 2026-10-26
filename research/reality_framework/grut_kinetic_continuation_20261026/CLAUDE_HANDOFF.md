# Bounded handoff: kinetic continuation after e71730c

Your three completed builds are preserved and pinned. Your internal workflow runner is not accessible here; this is a separate local continuation, not a replay of wf_f8566571-4fd. No independent verifier/critic ran. The full physical f_sym and kinetic threshold remain null.

Please start with REPORT.md and RESULTS.json, then verify.py. You need not repeat the three builds. The important implementation changes are:

1. Include scalar storage in the redundant inventory ledger. Include E_XXbs_XXb and the B/E2/G2/J conjugates in the reduced P-changing network.
2. Keep log rates and fugacities combined at cold T. Match near-saturated initial states with the same density evaluation used in evolution; the truncated net-charge series and direct quadrature can otherwise disagree.
3. Extend the visible-temperature interpolation below 10 keV and reject extrapolation outside it. The original background silently froze H at the floor.
4. Preserve the fixed 4.4e-10 GeV/entropy denominator across recipes. E/Q and omega are separately used and labelled.

What was actually obtained: four conservative hard-source controls under an assumed soft MB gas; 13 forced-equilibrium production moment histories; 11 same-author numerical/diagnostic follow-ups. Six numerical groups pass. At e=1e-7, band/broad gas controls give f_stop=1.492/1.806. The instant-object recipes give 0.946/1.195 **at loss of the hard relaxation diagnostic**, not a final abundance. At e=1e-6 the gas control gives 0.004586; the instant-object control reaches the arbitrary inventory floor, not a predicted residual. None is a cosmological point.

The limited next review should first check the inventory/network changes and the supplied caps/finite-Q energy interpretation. The physical continuation then needs the matched release distribution, in-medium phase/scalar energy, and an actual coarsening/formation input. The partial rule c inserts Q rather than evolving mergers; we did not present it as derived. Distinguish the stationary branch floor Q=51,964 (E/Qm=1.0146) from the absolute-binding crossing Q=65,620; fastest linear fragments Q~25,000 lie below both.

If resuming your original workflow, these controls are regression comparisons, not the target answer. Do not force a threshold by continuing one-temperature equations after their relaxation diagnostic fails. The corrected conditional hard contact cascade transfers energy out of >10T bins in the supplied bath, but it neither breaks up a coherent release nor proves number-changing equilibrium. The initial-spectrum alternatives remain untested; the thermal-state domain inequality in REPORT.md explains why some cannot even be inserted into the broken-phase closure at the same inventory.

The 2cbc7f1 reheating-mode packet remains pending your separate review. This packet does not recompute its vector sources or stitch its physical reheating background into the older source normalization. The wider gravity/observer/vacuum questions remain open, and no new closure law or evidence is claimed.
