#!/usr/bin/env python3
"""st_clocks_skeptic (Claude, own code): cross-check of st_clocks' N1 diagnostic (homogeneous Y_L/Y_H counts, E_X, N_X)
using CODEX's saved homogeneous fields (quiet*.npz, a different engine from st_clocks' rpreheat.c) and an own
implementation of the free-quanta decomposition about the instantaneous forced response:
  W = e^{-i arg(Sigma)/2} Y;  forced W_f: (M^2 + 2 mu |S|) Re W_f = Re J,  (M^2 - 2 mu |S|) Im W_f = Im J,
  J = -kappa X^3 e^{-i arg(Sigma)/2};  N_{H,L} = a^3 [(d dW_{re,im})^2 + w^2 (dW_{re,im})^2]/w   (dW = W - W_f).
The lambda_Y |Y|^2 Y term and d|Sigma|/dt are neglected in the forced response (both tiny here) [A].
Output: ../results/population_skeptic.json
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'results' / 'population_skeptic.json'
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
ST = json.loads((HERE.parents[1] / 'results' / 'clocks.json').read_text())['histories']
BG = json.loads((HERE.parent / 'results' / 'background_skeptic.json').read_text())['histories']
KAP, M, MU, LX = 0.0034586847211042523, 30.0, 7.159456159513794, 1e-3
mL = math.sqrt(M * M - math.sqrt(2) * MU * 40.0)


def main():
    out = {}
    for tag in ('quiet3', 'quiet5'):
        z = np.load(CODEX / (tag + '.npz'))
        t = z['times']; bg = z['bg']; a = np.exp(bg[:, 15])
        X = bg[:, 2] + 1j * bg[:, 3]; Y = bg[:, 4] + 1j * bg[:, 5]; S = bg[:, 6] + 1j * bg[:, 7]
        Xd = bg[:, 8] + 1j * bg[:, 9]; Yd = bg[:, 10] + 1j * bg[:, 11]; Sd = bg[:, 12] + 1j * bg[:, 13]
        th = np.angle(S); thd = np.imag(np.conj(S) * Sd) / np.abs(S) ** 2
        rot = np.exp(-0.5j * th)
        wH2 = M * M + 2 * MU * np.abs(S); wL2 = M * M - 2 * MU * np.abs(S)
        J = -KAP * X ** 3 * rot
        Jd = -KAP * 3 * X ** 2 * Xd * rot + (-0.5j * thd) * J
        Wf = J.real / wH2 + 1j * J.imag / wL2
        Wfd = Jd.real / wH2 + 1j * Jd.imag / wL2
        W = Y * rot; Wd = Yd * rot - 0.5j * thd * W
        dW = W - Wf; dWd = Wd - Wfd
        NH = a ** 3 * (dWd.real ** 2 + wH2 * dW.real ** 2) / np.sqrt(wH2)
        NL = a ** 3 * (dWd.imag ** 2 + wL2 * dW.imag ** 2) / np.sqrt(wL2)
        EX = a ** 4 * (np.abs(Xd) ** 2 + LX * np.abs(X) ** 4)
        sel = t >= 640
        r = dict(NYL_late=float(np.mean(NL[sel])), NYH_late=float(np.mean(NH[sel])), EX_late=float(np.mean(EX[sel])),
                 NYL_spread_late=float(np.ptp(NL[sel]) / np.mean(NL[sel])),
                 st_clocks_NYL=ST[tag]['cascade_YH_to_YL_AL']['NYL_late_own'],
                 st_clocks_NYH=ST[tag]['cascade_YH_to_YL_AL']['NYH_late_own'],
                 st_clocks_EX=ST[tag]['cascade_YH_to_YL_AL']['EX_late_own'])
        meanA = BG[tag]['sources']['YL']['mean_a_own']
        r['R_Y_homogeneous'] = (r['NYL_late'] + r['NYH_late']) * mL * meanA / r['EX_late']
        r['st_clocks_R_Y_homogeneous'] = ST[tag]['YL_energy_on_this_history_C']['R_Y_homogeneous']
        r['dictionary_NYL'] = 0.3344649635088846
        out[tag] = r
        print(tag, json.dumps(r, indent=1))
    OUT.write_text(json.dumps(out, indent=2) + '\n')


if __name__ == '__main__':
    main()
