"""fx_chain_skeptic: entropy per charge of the network's forced-equilibrium start versus the relaxed pre-onset
plasma (independent Hartree solver, ../hartree_indep.py).  Uses the copied network's own Bose state functions.
Prints, for starts at k x a_on with a fraction f of V(0) dumped as energy: T/m_X, s/n_Q, and the V(0) fraction
that makes the start isentropic (s/n_Q equal to the pre-onset Hartree value).  [C]"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('GRUT_KINETIC_WORKDIR', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/fx_chain_skeptic/kwork')
import m1_kin as MK
import ks_common as K
import thermo_state as TS
from scipy.optimize import brentq

MK.configure('m1')
HO = MK.HO
row = [r for r in HO['chain'] if r['label'] == 'heavy_release:tot_e1e6_band_4.97' and r['phi_in_equilibrium']][0]
R = row['R']; aon = row['onset_untilted']['a_over_aref']
HI = json.load(open(os.path.join(HERE, '..', '..', 'results', 'hartree_indep.json')))
S_pre = [r for r in HI['rows'] if abs(r['R'] - R) < 1e-6 and r['phi']][0]['S_per_net']
cap = 1 - 1e-6


def start(k, f):
    x0 = math.log(aon * k)
    rho = (K.RHO_C * (1 + R) * math.exp(-4 * x0) + f * K.V0) / K.MX ** 4
    nq = K.N_Q0 * math.exp(-3 * x0) / K.MX ** 3
    ini = TS.initial_equilibrium(nq, rho, cap, cap)
    t = ini['t']; u = ini['up']
    xp = TS.species(t, 1., u); xm = TS.species(t, 1., -u); s = TS.species(t, K.RS, 0.)
    S = xp['s'] + xm['s'] + s['s']
    return dict(k=k, f=f, t=t, s_per_net=S / nq, cond=ini['n_cond'] / nq, Xbar=xm['n'] / nq, s_quanta=s['n'] / nq,
                V0_per_charge=K.V0 / K.MX ** 4 / nq)


out = dict(S_pre_hartree=S_pre, rows=[])
print('pre-onset Hartree s/n_Q = %.5f' % S_pre)
for k in (0.84, 1.0, 1.2, 1.5, 2.0):
    for f in (0.0, 0.42, 1.0):
        r = start(k, f); out['rows'].append(r)
        print('k=%.2f f=%.2f  T/m=%.4f s/nQ=%.5f (ratio %.4f) cond=%.3f Xbar=%.3f s=%.3f V0/nQ=%.3f' % (
            k, f, r['t'], r['s_per_net'], r['s_per_net'] / S_pre, r['cond'], r['Xbar'], r['s_quanta'], r['V0_per_charge']))
    fi = brentq(lambda f: start(k, f)["s_per_net"] - S_pre, -0.5, 3.0, xtol=1e-6)
    out.setdefault('isentropic_V0_fraction', {})[str(k)] = fi
    print('   isentropic start at k=%.2f needs V(0) fraction %.3f' % (k, fi))
json.dump(out, open(os.path.join(HERE, '..', '..', 'results', 'start_entropy.json'), 'w'), indent=1)
