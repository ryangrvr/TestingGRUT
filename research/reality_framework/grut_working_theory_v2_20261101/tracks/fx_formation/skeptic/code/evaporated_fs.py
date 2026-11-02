"""fx_formation SKEPTIC: free streaming of the free-X relic in the evaporated (e = 1e-6) branch.  NR vapour at T_s,
v ~ 1/a afterwards; lambda_fs from skeptic_checks.lambda_fs.  [E]"""
import sys, os, math, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
out_path = sys.argv[1]
sys.argv = [sys.argv[0], out_path, os.path.join(HERE, '..', '..', 'results', 'eq_history.json')]
import skeptic_checks as sc, bessel_vapour as bv
out = {}
for lab, Yh in (('heavy_e1e6', 0.2020822595873346), ('pilot_e1e6', 0.27680265413087224)):
    Tev, Ts = bv.crossings(Yh)
    mu = bv.vapour_only_mu(Yh, Ts * 1.001); n = bv.vapour(mu, Ts * 1.001)['n']; Tv = sc.Tvis_of(n)
    out[lab] = dict(T_s=Ts, T_vis_at_Ts_GeV=Tv, v_rms_at_Ts=math.sqrt(3 * Ts), lambda_fs_kpc=sc.lambda_fs(Ts, Tv),
                    fugacity_at_Ts=math.exp((mu - 1) / Ts))
    print(lab, out[lab])
json.dump(out, open(out_path, 'w'), indent=1)
