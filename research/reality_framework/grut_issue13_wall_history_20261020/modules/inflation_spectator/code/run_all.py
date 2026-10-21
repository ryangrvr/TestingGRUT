"""Run every calculation of module inflation_spectator in dependency order (deterministic; no random draws
except the fixed-seed coupling samples of the xi-beta identity check). About 4-5 minutes on 2 threads."""
import os
import time

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import xi_running
import checks_1d
import stochastic_scan
import fp2d_checks
import post_inflation
import final_sign
import classification
import plot_svg
from common import dump

if __name__ == "__main__":
    times = {}
    for name, mod in (("xi_running", xi_running), ("checks_1d", checks_1d), ("stochastic_scan", stochastic_scan),
                      ("fp2d_checks", fp2d_checks), ("post_inflation", post_inflation),
                      ("final_sign", final_sign), ("classification", classification), ("plot_svg", plot_svg)):
        t = time.time()
        mod.main()
        times[name] = round(time.time() - t, 1)
        print(name, times[name], "s", flush=True)
    dump("runtimes.json", times)
