"""Audit (b): averaged export bookkeeping (own code).

Recomputes, from the scalar trace of the rerun (spatial means q0,q1,qphi,p0,p1,pphi, rho, nQ at every sample),
the coherent outer product <v0 v0^T> with trapezoid weights in conformal (ds) and cosmic (dtau = a ds) time,
compares with the AVERAGE_<tag>.npz arrays, reconstructs charge and energy with our own index map, and measures
what squaring the time-averaged mean would lose.
Usage: python3 -B audit_export.py <AVERAGE npz> <trace.npz> <out.json>
"""
import sys, json
import numpy as np

avg_path, trace_path, out = sys.argv[1:4]
A = np.load(avg_path); T = np.load(trace_path)
s = T['s']; vals = T['values']; names = list(T['names'])
col = {n: i for i, n in enumerate(names)}
# export ordering: v = (q0, q1, p0, p1, qphi, pphi)
v = np.stack([vals[:, col[n]] for n in ['q0', 'q1', 'p0', 'p1', 'qphi', 'pphi']], axis=1)
a_m, C = 8.869218738806, 1.9857962750634526
res = {}
for lo, hi in [(14., 24.), (24., 34.)]:
    sel = (s >= lo - 1e-10) & (s <= hi + 1e-10); ss = s[sel]; vv = v[sel]
    trap = np.ones(len(ss)); trap[0] = trap[-1] = .5
    for rule in ['conformal', 'cosmic']:
        w = trap * (1. if rule == 'conformal' else (a_m + C * ss)); w = w / w.sum()
        key = f'{lo:g}_{hi:g}_{rule}'
        outer = np.einsum('t,ti,tj->ij', w, vv, vv); meanv = w @ vv; sq = np.outer(meanv, meanv)
        Mx = A[key + '_mean_outer']; Cv = A[key + '_cov']; Ck = A[key + '_covk']; U = A[key + '_potential']
        rho = float(A[key + '_rho']); Q = float(A[key + '_Q'])
        rho_tr = float(w @ vals[sel, col['rho']]); Q_tr = float(w @ vals[sel, col['nQ']])
        # own reconstruction
        Qc = Cv[:, 0, 3].sum() - Cv[:, 1, 2].sum(); Qm = Mx[0, 3] - Mx[1, 2]
        Ekin_c = .5 * (Cv[:, 2, 2] + Cv[:, 3, 3] + Cv[:, 5, 5]).sum(); Egrad_c = .5 * (Ck[:, 0, 0] + Ck[:, 1, 1] + Ck[:, 4, 4]).sum()
        Ekin_m = .5 * (Mx[2, 2] + Mx[3, 3] + Mx[5, 5])
        E = Ekin_c + Egrad_c + Ekin_m + U.sum()
        Ekin_sq = .5 * (sq[2, 2] + sq[3, 3] + sq[5, 5]); Q_sq = sq[0, 3] - sq[1, 2]
        res[key] = dict(mean_outer_vs_trace_max_rel=float(np.max(abs(Mx - outer)) / np.max(abs(Mx))),
                        rho_vs_trace_rel=abs(rho / rho_tr - 1), Q_vs_trace_rel=abs(Q / Q_tr - 1),
                        charge_reconstruction_rel=abs((Qc + Qm) / Q - 1), energy_reconstruction_rel=abs(E / rho - 1),
                        coherent_charge_share=float(Qm / Q), connected_charge_share=float(Qc / Q),
                        coherent_kinetic_energy_share=float(Ekin_m / rho), potential_share=float(U.sum() / rho),
                        connected_kinetic_plus_gradient_share=float((Ekin_c + Egrad_c) / rho),
                        squared_mean_coherent_charge_share=float(Q_sq / Q), squared_mean_coherent_kinetic_share=float(Ekin_sq / rho),
                        squared_mean_loses_charge_fraction=float((Qm - Q_sq) / Q), squared_mean_loses_energy_fraction=float((Ekin_m - Ekin_sq) / rho),
                        zero_k_bin_cov_max=float(np.max(abs(Cv[0]))), weight_ratio_end_to_start=float((a_m + C * ss[-1]) / (a_m + C * ss[0])))
# cosmic vs conformal difference of headline numbers
for lo, hi in [(14., 24.), (24., 34.)]:
    a = res[f'{lo:g}_{hi:g}_conformal']; b = res[f'{lo:g}_{hi:g}_cosmic']
    res[f'{lo:g}_{hi:g}_cosmic_minus_conformal'] = dict(coherent_charge_share=b['coherent_charge_share'] - a['coherent_charge_share'],
                                                       coherent_kinetic_energy_share=b['coherent_kinetic_energy_share'] - a['coherent_kinetic_energy_share'])
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
