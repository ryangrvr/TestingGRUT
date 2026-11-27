"""Matched continuum low waves, radial lattice high waves and Born graphs."""
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simpson

from branch import ARCHIVE

HERE = Path(__file__).resolve().parent
ENERGY_ARCHIVE = ARCHIVE / "grut_candidate_1_energy_20261019"
sys.path.insert(0, str(ENERGY_ARCHIVE))
import graphs
import spectral as archived_jost

PROFILES = {}


def jost_worker(job):
    name, ell, zeta, radius, tol, origin = job
    profile = PROFILES[name]
    archived_jost.SPL, archived_jost.OMEGA = profile.spline, profile.omega
    return archived_jost.run_point(ell, zeta, R=radius, tol=tol, r0=origin, label=name)


def quadrature(edges, n):
    x, w = leggauss(n)
    zeta = [(b+a)/2+(b-a)*x/2 for a, b in zip(edges, edges[1:])]
    weights = [(b-a)*w/2 for a, b in zip(edges, edges[1:])]
    return np.concatenate(zeta), np.concatenate(weights)


def low_weights(w):
    return np.array([1, 3])[:, None]*w/(2*np.pi)


def high_weights(w, cutoff):
    return (2*np.arange(2, cutoff+1)+1)[:, None]*w/(2*np.pi)


def low_energy(arr, z, w, end):
    zero = np.log(z*z/(z*z+1))
    restored = end*np.log(end*end/(end*end+1))-2*np.arctan(end)
    return float(np.sum((arr-zero)*low_weights(w)) + 4*restored/(2*np.pi))


def difference(a, b, weights):
    return math.fsum(((a-b)*weights).ravel())


def absolute_variation(a, b, weights):
    return math.fsum((abs(a-b)*weights).ravel())


class Spectrum:
    def __init__(self, protocol, output, started):
        self.settings = protocol["determinant"]
        self.limits = protocol["stopping"]
        self.deadline = started + self.limits["wall_seconds"]
        self.output = output
        self.fd_batches = self.jost_points = 0
        self.tmp = tempfile.TemporaryDirectory(prefix="grut-c1-")
        self.bin = Path(self.tmp.name) / "radial"
        self.box = self.settings["fd_box"]
        self.support = self.settings["fd_support"]
        self.cutoff = self.settings["ell_max"]
        self.low_ells = self.settings["jost_ell_primary"]
        self.env = {**os.environ, "OMP_NUM_THREADS": str(self.limits["worker_limit"])}
        subprocess.run(["g++", "-O3", "-std=c++17", "-fopenmp", str(HERE / "radial.cc"), "-o", str(self.bin)],
                       check=True, timeout=self.remaining(), env=self.env)

    def remaining(self):
        left = self.deadline-time.monotonic()
        if left <= 0:
            raise TimeoutError("Frozen production wall budget exhausted")
        return left

    def note(self, **data):
        with (self.output / "execution.jsonl").open("a") as stream:
            stream.write(json.dumps(data, default=float) + "\n")

    def jost(self, name, probes, configuration, label):
        jobs = [(name, int(ell), float(zeta), *configuration) for ell, zeta in probes]
        self.jost_points += len(jobs)
        if self.jost_points > self.limits["maximum_jost_points"]:
            raise RuntimeError("Frozen Jost-point budget exhausted")
        self.note(event="jost_start", profile=name, label=label, count=len(jobs))
        answer = {}
        with mp.get_context("fork").Pool(self.limits["worker_limit"]) as pool:
            results = pool.imap_unordered(jost_worker, jobs, chunksize=1)
            for _ in jobs:
                rec = results.next(timeout=self.remaining())
                point = rec["ell"], rec["zeta"]
                ok = rec["success"] and abs(rec["logJost_phase_mod_2pi"]) < self.settings["max_abs_jost_phase"]
                self.note(event="jost_point", stage=label, valid=bool(ok), **rec)
                if not ok:
                    raise ValueError(f"Jost failure at {name}, {label}, {point}")
                answer[point] = rec
        return [answer[int(ell), float(zeta)] for ell, zeta in probes]

    def low(self, name, zeta, configuration, label):
        probes = [(ell, z) for ell in self.low_ells for z in zeta]
        rec = self.jost(name, probes, configuration, label)
        return np.array([row["total_remainder"] for row in rec]).reshape(len(self.low_ells), len(zeta))

    def fd(self, profile, h, probes, label, *, support=None, box=None):
        self.fd_batches += 1
        if self.fd_batches > self.limits["maximum_fd_batches"]:
            raise RuntimeError("Frozen radial FD batch budget exhausted")
        box = self.box if box is None else box
        support = self.support if support is None else support
        n = round(box/h) - 1
        r = np.arange(1, n+1)*h
        fs = profile.spline(r)[[0, 2]].T
        fs[r >= support] = [0., 1.]
        request = Path(self.tmp.name) / "lattice.in"
        with request.open("w") as stream:
            stream.write(f"{n} {len(probes)} {h:.17g} {profile.omega:.17g}\n")
            np.savetxt(stream, probes, fmt=["%d", "%.17g"])
            np.savetxt(stream, fs, fmt="%.17g")
        self.note(event="fd_start", profile=profile.label, h=h, label=label, support=support, box=box, count=len(probes))
        run = subprocess.run([str(self.bin), str(request)], capture_output=True, text=True,
                             timeout=self.remaining(), env=self.env)
        if run.returncode:
            self.note(event="fd_failed", profile=profile.label, label=label, returncode=run.returncode,
                      stderr=run.stderr, last_rows=run.stdout.splitlines()[-10:])
            raise RuntimeError(f"Lattice failure: {profile.label}, {label}")
        rows = np.fromstring(run.stdout, sep=" ").reshape(-1, 3)
        if rows.shape != (len(probes), 3) or not np.all(np.isfinite(rows)) or np.min(rows[:, 2]) <= 0:
            raise RuntimeError("Invalid determinant ratios on radial grid")
        self.note(event="fd_complete", profile=profile.label, h=h, label=label, minimum_pivot_ratio=np.min(rows[:, 2]))
        return np.sum(rows[:, :2], axis=1)

    def high(self, profile, h, zeta, label, *, cutoff=None, extra=(), **kwargs):
        cutoff = self.cutoff if cutoff is None else cutoff
        probes = [(ell, z) for ell in range(2, cutoff+1) for z in zeta]
        mains = len(probes)
        out = self.fd(profile, h, probes + list(extra), label, **kwargs)
        return out[:mains].reshape(cutoff-1, len(zeta)), out[mains:]


def born_graph(profile, configuration):
    rmax, n_r, kmax, n_k, n_x = configuration
    r, k = np.linspace(0, rmax, n_r), np.linspace(0, kmax, n_k)
    f, s = profile.spline(r)[[0, 2]]
    V = graphs.potentials(f, s)
    kernel, cuts = graphs.kernels(k, profile.omega, n_x)
    vhat = graphs.transform(r, V, k)
    one = graphs.first_graph(4*np.pi*simpson(r[:, None]**2 * V, x=r, axis=0))
    two = graphs.second_graph(k, vhat, kernel)
    fourier = simpson(k[:, None]**2 * vhat**2/(2*np.pi**2), x=k, axis=0)
    spatial = 4*np.pi*simpson(r[:, None]**2 * V**2, x=r, axis=0)
    total = math.fsum(np.r_[one, two])
    return {"total": total, "tadpoles": one, "bubbles": two,
            "max_parseval_error": max(abs(fourier/spatial-1)), "min_bubble_denominator": min(cuts),
            "k": k, "vhat": vhat, "kernel": kernel}


def graph_difference(first, second):
    k = first["k"][:, None]
    if not np.array_equal(first["k"], second["k"]):
        raise ValueError("Correlated subtraction requires identical momentum nodes")
    a, b = first["vhat"], second["vhat"]
    f, g = first["kernel"], second["kernel"]
    paired = (a-b)*(a+b)*f+b*b*(f-g)
    kernel_weight = np.array([1, 2, 2, 2, 1])*k*k/(2*np.pi**2)
    deltas = -simpson(kernel_weight*paired, x=k[:, 0], axis=0)/(64*np.pi**2)
    return math.fsum(np.r_[first["tadpoles"]-second["tadpoles"], deltas])
