# Reproducing the calculations

Run these commands from this packet's directory in a full checkout of the parent repository. Python 3.12, NumPy 2.3, SciPy 1.17 and Matplotlib 3.10 were used; exact versions and the optional CLASS source commit are in [ENVIRONMENT.json](ENVIRONMENT.json). The script imports a few files from sibling research packets, and [INPUT_PINS.json](INPUT_PINS.json) lists their exact Git blob IDs.

```sh
./reproduce.sh --verify
```

This checks the 18 inherited input blobs, scientific invariants and saved result arrays, and the artifact manifest. These checks are same-author consistency tests, not independent physics reproduction.

```sh
./reproduce.sh --all
```

The full run recalculates scale checks, both fixed-curvature branches, idealized phase modes, four 200-time reheating histories, both 800-time histories, five finite-difference histories, the free-streaming/PSD calculation, the gas rate, the quadrature audit, and the figure. It can take several minutes on a single core. It overwrites numerical output; CPU times and small floating-point differences make a byte-for-byte manifest match unlikely afterward. Run `python code/manifest.py` to record the new outputs and `python code/manifest.py --check` to verify them. The archived CLASS results remain as published unless CLASS is separately rerun.

For the optional CLASS calculations, use a separate checkout of [the pinned CLASS commit](https://github.com/lesgourg/class_public/tree/814ce59c380a5cccd0bed1747a384ff7bfc1f9cf) and build `make class -j2`. First regenerate the input PSDs using `python code/free_stream.py`. The following runs all successful archived cases before collecting their results; the actual-spectrum 120-node run and RK controls can be slow:

```sh
python code/class_power.py --run --class-binary /absolute/path/to/class \
  --cases cold,coldtight,coldrk,coldrk20,warm96fast20,warm96q60,warm96rk,warm96rkq60,warm96rkq120,warm96maxwell,warm96maxwellrk
```

CLASS writes temporary parameter files, logs, and spectra under this packet's ignored `.work/class_runs/`; the collector writes `class_results.json` and `class_spectra.npz`. A partial list of cases collects the runs present in that working directory, so it will replace the archived summary with that partial set; use a separate checkout or preserve the archives first. The failed or intentionally interrupted exploratory runs are described in [CLASS_ATTEMPTS.md](CLASS_ATTEMPTS.md), not reproduced by the success list. Positive small differences and inconsistent resolution controls must not be interpreted as a measured signal.

The test suite and reproducibility scripts cannot turn the chosen spectator, one-point phase-space, source, or collision assumptions into derived physical facts.
