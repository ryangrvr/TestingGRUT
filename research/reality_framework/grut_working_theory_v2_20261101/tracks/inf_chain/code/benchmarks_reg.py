#!/usr/bin/env python3
"""inf_chain (Claude, own code): registry of the benchmarks pushed through the chain.  Every entry is a model choice.

  HT6t       inf_potential's primary benchmark (tilted p = 6 hilltop) at H_i = 1.29e11 GeV, slice and holder at S1 (one
             e-fold before eps_H = 1), Gamma_I = 0.5 H_i (the engine convention the potential track used).
  HT6t_H150  the same family re-calibrated by inf_potential at H_i = 1.50e11 (its Q-matter-calibration row).
  HT6t_Send  HT6t with the holder calibrated at eps_H = 1 (mu_eff = 10i there; g_I x 3.29), slice still S1.
  T1         inf_potential's stress test (alpha-attractor T-model, alpha = 6.1e-5) at 1.29e11, S1, Gamma_I = 0.5.
  TB1, TB2   inf_reheat's quiet T-model members alpha = 0.3 (B1) and 1 (B2), here with the S1 slice/holder convention;
             V0 from A_s at N_* (iterated in chain.py); the engine runs Gamma_I = 0 (N_X is Gamma-independent for
             Gamma_I <~ 1e-3 [Q inf_reheat]); the slow decay is added in the fluid continuation.
"""
import ch_common as ch


def get_bench(key, Nstar=None, Hi_index=None):
    if key == 'HT6t':
        return ch.bench_from_inf_potential('HT6t', 0)
    if key == 'HT6t_H150':
        b = ch.bench_from_inf_potential('HT6t', 1); b.name = key; return b
    if key == 'HT6t_Send':
        b = ch.bench_from_inf_potential('HT6t', 0, Ncal=0.0); b.name = key; return b
    if key == 'T1':
        return ch.bench_from_inf_potential('T1', 0)
    if key == 'TB1':
        return ch.tmodel_bench(0.3, Nstar=Nstar, name='TB1')
    if key == 'TB2':
        return ch.tmodel_bench(1.0, Nstar=Nstar, name='TB2')
    if key.startswith('T1_cal'):
        import json, inf_common as ic
        d = json.load(open(ch.RES / ('recal_%s.json' % key)))
        m = ic.Model('T', **d['params'])
        return ch.Bench(key, m, d['V0_MP4'], Nslice=1.0, gam=0.5, note='T1 recalibrated at H_i = %.4e' % d['H_i_GeV'])
    if key.startswith('HT6t_cal'):
        # HT6t re-calibrated by chain.py at another H_i (shape stored in results/recal_<key>.json)
        import json, inf_common as ic
        d = json.load(open(ch.RES / ('recal_%s.json' % key)))
        m = ic.Model('HT', **d['params'])
        b = ch.Bench(key, m, d['V0_MP4'], Nslice=1.0, gam=0.5, note='HT6t recalibrated at H_i = %.4e' % d['H_i_GeV'])
        return b
    raise KeyError(key)
