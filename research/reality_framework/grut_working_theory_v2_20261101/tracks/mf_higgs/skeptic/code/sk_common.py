#!/usr/bin/env python3
"""mf_higgs SKEPTIC (key mf_higgs_skeptic; Claude, multi-field pass, BRIEF_MULTIFIELD.md; GitHub issue #13 context).

Common helpers.  Every choice here is a model choice (benchmark, xi_H, criterion, plasma model); none follows from
closure.  Independent of the track's code (tracks/mf_higgs/code is NOT imported).  Shared on purpose, read-only:
  * SM running and lambda_eff: tracks/higgs_reheat/code/smpot.py (3-loop running, flat 2-loop lambda_eff, xi running,
    one-loop de Sitter potential of 1804.02020) -- the brief says not to redo the SM.
  * Histories: tracks/mf_hybrid/results/bench_HFa.json (as exported) and the mf_chain rebuilt histories
    tracks/mf_chain/results/histories/bench_it*.json (a different code path from the track's own fluid).
Units: H_i = 1 (Hubble rate at the slice sigma = sigma_c, a = 1 there), fields canonical in H_i.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SKDIR = HERE.parent
RES = SKDIR / 'results'
TRACKS = SKDIR.parent.parent
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_higgs_skeptic')
SP.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(TRACKS / 'higgs_reheat' / 'code'))
import smpot as S   # noqa: E402  read-only SM input

MPL = 2.435e18
PDG = dict(Mt=172.57, sMt=0.29, Mh=125.20, sMh=0.11, als=0.1180, sals=0.0009)
LN_K_H0 = math.log(0.05 / 2.248e-4)   # ln(k_*/(a0 H0)), H0 = 67.4 km/s/Mpc


def jdump(obj, path):
    def conv(z):
        if isinstance(z, (np.floating, np.integer)):
            return float(z)
        if isinstance(z, np.ndarray):
            return z.tolist()
        if isinstance(z, (bool, np.bool_)):
            return bool(z)
        return str(z)
    json.dump(obj, open(path, 'w'), indent=1, default=conv)


_SM = {}


def sm(Mt, Mh=PDG['Mh'], als=PDG['als']):
    k = (round(Mt, 4), round(Mh, 3), round(als, 5))
    if k not in _SM:
        _SM[k] = S.SM(k[0], Mh=k[1], als=k[2])
    return _SM[k]


class Flat:
    """flat 2-loop lambda_eff(h) table (own cache), canonical conversion at H_i, h_max, V''."""

    def __init__(self, Mt, Mh=PDG['Mh'], als=PDG['als']):
        key = SP / 'flat' / f'flat_{Mt:.4f}_{Mh:.3f}_{als:.5f}.npz'
        key.parent.mkdir(exist_ok=True)
        if key.exists():
            z = np.load(key)
            self.lh, self.le = z['lh'], z['le']
            self.lmu, self.G, self.Ixi = z['lmu'], z['G'], z['Ixi']
        else:
            s = sm(Mt, Mh, als)
            hs = np.exp(np.linspace(math.log(1e8), math.log(1e18), 1001))
            self.lh = np.log(hs)
            self.le = s.lam_eff_flat(hs)
            mus = np.exp(np.linspace(math.log(1e2), math.log(1e18), 400))
            c = s.cpl(mus)
            self.lmu, self.G, self.Ixi = np.log(mus), np.asarray(c['Gam']), np.asarray(c['Ixi'])
            self.Ixi = self.Ixi - float(s.cpl(Mt)['Ixi'])
            np.savez(key, lh=self.lh, le=self.le, lmu=self.lmu, G=self.G, Ixi=self.Ixi)
        V = np.exp(self.le) if False else self.le * np.exp(4 * self.lh) / 4
        m = self.lh < math.log(1e15)
        k = int(np.argmax(np.where(m, V, -np.inf)))
        # refine the maximum by a parabola in ln h
        x, y = self.lh[k - 1:k + 2], V[k - 1:k + 2]
        A = np.polyfit(x, y, 2)
        self.h_max = float(math.exp(-A[1] / (2 * A[0])))

    def lam(self, h):
        return np.interp(np.log(h), self.lh, self.le)

    def Gam(self, mu):
        return float(np.interp(math.log(mu), self.lmu, self.G))

    def xi(self, mu, xi_EW):
        return 1 / 6 + (xi_EW - 1 / 6) * math.exp(float(np.interp(math.log(mu), self.lmu, self.Ixi)))

    def xi_EW(self, mu, xi_mu):
        return 1 / 6 + (xi_mu - 1 / 6) * math.exp(-float(np.interp(math.log(mu), self.lmu, self.Ixi)))

    def V2_canon(self, h0, H_GeV):
        """V''(h0) (flat, V = lam_eff h0^4/4) converted to the canonical field phi = e^{Gamma(H_i)} h0, in H_i^2"""
        d = 1e-3
        V = lambda x: self.lam(x) * x ** 4 / 4
        v2 = (V(h0 * (1 + d)) - 2 * V(h0) + V(h0 * (1 - d))) / (h0 * d) ** 2
        eG = math.exp(self.Gam(H_GeV))
        return v2 / eG ** 2 / H_GeV ** 2


def load_hist(path):
    d = json.load(open(path))
    h = {k: np.array([np.nan if x is None else x for x in v], dtype=float) for k, v in d['history'].items()}
    return d, h
