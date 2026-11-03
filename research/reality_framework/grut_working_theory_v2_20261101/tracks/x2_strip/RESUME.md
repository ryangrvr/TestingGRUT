# x2_strip: paused (owner's weekly limit)

**State.** The pre-registration was committed first (6d836ca). The fresh scans on the two computed nodes were
stopped partway: see the key counts in `results/scan_x2_hLF84_mid.json` and `results/scan_x2_hLF_mid.json`. The [A]
nodes (hLF55_mid, hLF270_mid) have not started. **No scan output has been analysed or inspected**, so the test is
still blind.

**To resume.** The driver skips keys already done:

    cd <scratch>/wt2/x2_strip      # or any directory
    python3 -B <repo>/research/reality_framework/grut_working_theory_v2_20261101/tracks/x2_strip/code/drive.py hLF84_mid  > d84.log &
    python3 -B <repo>/.../x2_strip/code/drive.py hLF_mid > d130.log
    # then hLF55_mid and hLF270_mid, for the extended reading only

The engine scratch goes to `wt2/x2_strip/sp_<node>`. The final histories are read from `wt2/lf_chain/sp0/hist`
(scratch). If scratch is gone, rebuild them with `lf_chain/code/build_hist.py`.

**Then analyse with exactly the PREREGISTRATION.md rule:**
- use `lf_chain/code/fp_lf.py`'s `nd()` and `root(lo, hi, 'free X')`, read-only, on the two computed nodes;
- a key is a strict M2-X point iff the root is self-consistent, 0 ≤ w ≤ 1, and 1.03×10⁻⁶ ≤ e_crit ≤ 4.12×10⁻⁶;
- report the strip (consecutive c0 per δ) or "no M2-X strip within the simulated range";
- report the extended reading and the M2-Q points separately.

**Also done in this session** (committed): the co-genesis test (`tracks/cg_cogenesis`, all four candidates fail) and
the emergence ledger v0.
