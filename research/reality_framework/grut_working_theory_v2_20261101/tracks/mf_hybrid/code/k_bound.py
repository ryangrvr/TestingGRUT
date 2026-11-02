import math, numpy as np, bench as B, reheat as RH, mfh_common as mc
R = B.run('HFa')
rec, spec, p = R['rec'], R['spec'], R['p']
jsw = len(rec['t']) - 1
Esig, Etot = float(rec['E_sig'][jsw]), float(rec['E_tot'][jsw])
class Pert(RH.Channel):
    def rate(self, rho):
        return self.Gp
cs = [Pert(c.m, c.M, c.f, dict(A=[1, 2], eps=[1, 1])) for c in R['chans_sig']]
cc = [Pert(c.m, c.M, c.f, dict(A=[1, 2], eps=[1, 1])) for c in R['chans_chi']]
st = dict(t=float(rec['t'][jsw]), lna=float(rec['lna'][jsw]), rho_sig=Esig, rho_chi=Etot - Esig, rho_r=float(rec['rho_r'][jsw]),
          s=float(rec['s_over_sstar'][jsw]), sd=float(rec['sdot_over_sstar'][jsw]), E_N=5000.0)
fl = RH.run_fluid(st, p['mP'], cs, cc, spec['H_GeV'], lambda t, rc, Rb: 9.0, t_end=3e5)
print('K with only perturbative widths after t = 0.6:', fl['K'][-1] / R['Hi'], ' (nominal', R['K'], ')  f_chi at switch', (Etot - Esig) / (Etot + rec['rho_r'][jsw]))
