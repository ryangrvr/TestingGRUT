"""M1 frozen-core, radiation-dominated expansion test.

All distributions are particle counts per reference core. The inherited
collision model and source are unmodified; expansion is a conservative
first-order upwind transport in physical energy. Finite grid diffusion is
quantified, not presumed absent. No late-core radiation is included.
"""
from pathlib import Path
import argparse
import json
import os
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
RF = ROOT.parent
PARENT = RF / 'grut_m1_kinetics_20261028'
os.environ['GRUT_M1_KINETICS_WORK'] = str(ROOT / '.work' / 'inherited')
sys.path.insert(0, str(PARENT / 'code'))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from transport import Model, diagnostics, grid, initial, deposit, daughters
from setup import EREF, M, OMEGA, QREF


def dump(name, payload):
    (ROOT / name).write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')


def redshift_matrix(e):
    """dN/d ln a. Rate chosen to reproduce dE/d ln a=-(E²-1)/E.

    Lower boundary is reflecting because first energy is >1. Count that floor
    separately. All columns conserve number; each hop loses v_i= p_i²/E_i.
    """
    v = (e * e - 1.) / e
    hop = np.zeros_like(e)
    hop[1:] = v[1:] / np.diff(e)
    return diags([hop[1:], -hop], [1, 0], shape=(len(e), len(e))).tocsr(), v


def radial_source(e, tag):
    """Exploratory compact signed radial spectra, with exclusions disclosed.

    The windows are treated as an instantaneous source after t=10000/mX.
    Negative aggregate bins cannot represent positive classical occupation;
    they and unsupported bins are excluded and separately recorded. This
    deliberate prescription is not a proof of positivity or a source match.
    """
    plus = np.zeros(len(e))
    minus = np.zeros(len(e))
    audit = {}
    early_only = tag.endswith('early')
    base_tag = tag.removesuffix('early')
    windows = ('0_1000',) if early_only else ('0_1000', '1000_10000')
    for window in windows:
        p = ROOT / 'sources' / f'ringdown_{base_tag}_spectrum_{window}.npz'
        data = np.load(p)
        for species in ('plus', 'minus', 'neutral'):
            w, en = data[species+'_N_signed'], data[species+'_EN_signed']
            neg = w < 0
            good = w > 0
            energy = np.divide(en, w, out=np.zeros_like(w), where=good)
            low = 1. if species != 'neutral' else np.sqrt(6.)
            valid = good & (energy >= low) & np.isfinite(energy)
            keep = valid & (energy >= e[0]) & (energy <= e[-1])
            item = {'signed_negative_bins_count':int(neg.sum()),
                    'discarded_negative_signed_number':float(w[neg].sum()),
                    'discarded_negative_signed_energy_mX':float(en[neg].sum()),
                    'invalid_positive_number':float(w[good & ~valid].sum()),
                    'out_of_range_positive_number':float(w[valid & ~keep].sum()),
                    'out_of_range_positive_energy_mX':float(en[valid & ~keep].sum()),
                    'internal_negative_weight':float(data[species+'_N_negative_part'].sum())}
            if species == 'neutral':
                missing_n = missing_e = 0.
                for ev, weight in zip(energy[keep], w[keep]):
                    dd, fraction = daughters(e, ev)
                    plus += weight * fraction * dd
                    minus += weight * fraction * dd
                    missing_n += weight * (1-fraction)
                    missing_e += weight * ev * (1-fraction)
                item['neutral_decay_unsupported_number'] = float(missing_n)
                item['neutral_decay_unsupported_energy_mX'] = float(missing_e)
            elif species == 'plus':
                plus += deposit(e, energy[keep], w[keep])
            else:
                minus += deposit(e, energy[keep], w[keep])
            audit[window + '_' + species] = item
    return np.r_[plus, minus], audit


def clock():
    pt = json.loads((RF / 'grut_matched_candidate_20261027' / 'CANDIDATE_M1.json').read_text())
    T, g = .02, 10.75
    Y = pt['final_population_matching_hypothesis']['YQ_target']
    entropy = 2 * np.pi**2 / 45 * g * T**3
    H = np.sqrt(np.pi**2 * g / 90) * T**2 / 2.435e18
    n0 = Y * entropy / QREF
    epsilon = H * M**2 / n0
    return {'temperature_GeV': T, 'YQ': Y, 'g_rho': g, 'g_s': g,
            'mX_GeV': M, 'ncore0_GeV3': n0, 'H0_GeV': H,
            'H0_mX2_over_ncore0': epsilon, 'tau_infinity': 1 / epsilon}


class Expanding:
    def __init__(self, e, mode, angles, epsilon):
        self.e, self.n, self.epsilon = e, len(e), epsilon
        self.model = Model(e, mode, angles)
        self.R, self.v = redshift_matrix(e)
        self.Rd = self.R.toarray()
        self.work = self.v.copy()
        self.work[0] = 0.  # tracked holding boundary, not an energy sink

    def rhs(self, x, y):
        n = self.n
        multiplier = np.exp(-x) / self.epsilon
        c, omitted = self.model.rates(y[:2*n])
        z = y[:2*n]
        red = np.r_[self.R @ z[:n], self.R @ z[n:]]
        work = self.work @ (z[:n] + z[n:])
        return np.r_[multiplier * c + red, multiplier * omitted, work]

    def jac(self, x, y):
        n = self.n
        z = np.zeros((2*n+2, 2*n+2))
        z[:2*n, :2*n] = np.exp(-x) / self.epsilon * self.model.rates(y[:2*n], True)
        z[:n,:n] += self.Rd
        z[n:2*n,n:2*n] += self.Rd
        z[-1,:n] = self.work
        z[-1,n:2*n] = self.work
        # The omitted-event integral is passive, as in the inherited solver.
        return z

    def observables(self, x, y):
        n, e = self.n, self.e
        base = diagnostics(e, np.r_[y[:2*n], y[-2]])
        a = float(np.exp(x))
        base.update({'a_over_a0': a, 'temperature_GeV_if_constant_gs': .02/a,
                     'tau': float((1 - 1/a) / self.epsilon),
                     'work_redshift_mX': float(y[-1]),
                     'total_tangent_plus_work_mX': float(base['total_tangent_energy_mX'] + y[-1]),
                     'particles_on_lowest_energy_node': float(y[0] + y[n]),
                     'floor_work_rate_mX_per_lna_if_unblocked': float(self.v[0] * (y[0] + y[n]))})
        return base


def free_stream_control(e, y0, ratios):
    """Compare exact characteristics with upwind dynamics for this actual spectrum."""
    n = len(e)
    R, v = redshift_matrix(e)
    combined = y0[:n] + y0[n:2*n]
    out = []
    for a in ratios:
        approx_plus = expm_multiply(R * np.log(a), y0[:n])
        approx_minus = expm_multiply(R * np.log(a), y0[n:2*n])
        approx = approx_plus + approx_minus
        exact_E = np.sqrt(1 + (e*e-1)/(a*a))
        exact_kin = float((exact_E - 1) @ combined)
        numeric_kin = float((e - 1) @ approx)
        threshold = 1 + 2 * OMEGA
        out.append({'a_over_a0':a,
                    'exact_kinetic_mX':exact_kin,
                    'upwind_kinetic_mX':numeric_kin,
                    'relative_kinetic_error':(numeric_kin-exact_kin)/exact_kin,
                    'exact_positive_above_reverse_threshold':float(y0[:n][exact_E > threshold].sum()),
                    'upwind_positive_above_reverse_threshold':float(approx_plus[e > threshold].sum()),
                    'floor_particle_count':float(approx[0]),
                    'number_absolute_error':float(approx.sum()-combined.sum())})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shift-bins',type=int,default=32)
    ap.add_argument('--angles',type=int,default=8)
    ap.add_argument('--emax',type=float,default=8.)
    ap.add_argument('--amax',type=float,default=10.)
    ap.add_argument('--rtol',type=float,default=2e-7)
    ap.add_argument('--mode',choices=['core','resonant','full'],default='full')
    ap.add_argument('--source',choices=['parent','h0025','h0025early','h005','h005early','h01'],default='parent')
    ap.add_argument('--tag',required=True)
    args = ap.parse_args()
    t0 = time.monotonic()
    e = grid(0, args.emax, 'translated', args.shift_bins)
    source, unsupported = initial(e) if args.source == 'parent' else radial_source(e,args.source)
    initial_y = np.r_[source, 0., 0.]
    clk = clock()
    eps = clk['H0_mX2_over_ncore0']
    x1000 = -np.log(1-1000*eps)
    xvals = np.unique(np.r_[0, -np.log(1-np.array([.1,1,10,100,1000,10000,100000])*eps),
                            np.log([1.01,1.1,2.,5.,args.amax])])
    xvals = xvals[(xvals>=0) & (xvals<=np.log(args.amax)+1e-14)]
    model = Expanding(e,args.mode,args.angles,eps)
    print('built',args.tag,len(e),'energy nodes in',round(time.monotonic()-t0,2),'s',flush=True)
    result = solve_ivp(model.rhs, (0.,np.log(args.amax)), initial_y, method='BDF',
                       jac=model.jac, rtol=args.rtol, atol=1e-9, t_eval=xvals)
    if not result.success:
        raise RuntimeError(result.message)
    rows = [model.observables(x, y) for x, y in zip(result.t,result.y.T)]
    ledgers = np.array([r['total_tangent_plus_work_mX'] for r in rows])
    # Work identity is checked against the same discrete upwind operator;
    # the independent characteristic comparison checks discretization error.
    checks = {'energy_plus_work_max_absolute_drift_mX':float(np.max(abs(ledgers-ledgers[0]))),
              'energy_plus_work_max_relative_drift':float(np.max(abs(ledgers/ledgers[0]-1))),
              'minimum_saved_population':float(result.y[:2*len(e)].min()),
              'initial_unsupported_spectrum':unsupported,
              'core_event_energy_error':model.model.core_error,
              'redshift_column_sum_max':float(np.max(abs(np.asarray(model.R.sum(axis=0))))),
              'free_streaming_start_source':free_stream_control(e,source,[1.001,1.1,2.,5.,10.]),
              'omitted_events_at_end':float(result.y[-2,-1])}
    np.savez_compressed(ROOT/(args.tag+'.npz'),E=e,x=result.t,Y=result.y)
    dump(args.tag+'.json',{'configuration':vars(args),'clock':clk,
         'solver':{'success':True,'nfev':result.nfev,'njev':result.njev,'nlu':result.nlu,
                   'seconds':time.monotonic()-t0},'checks':checks,'snapshots':rows,
         'scope':'Frozen scattering profile and frequency; no evolving core excitation, formation, source after initial pulse, or entropy changes.'})
    print('done',args.tag,'nfev',result.nfev,'nlu',result.nlu,'seconds',round(time.monotonic()-t0,2),flush=True)
    for row in rows:
        print('a, tau, Nplus, Nminus, Q, Ekin, work, floor:',
              *[round(row[k],6) for k in ['a_over_a0','tau','Nplus','Nminus','Qcore','gas_kinetic_mX','work_redshift_mX','particles_on_lowest_energy_node']], flush=True)

if __name__ == '__main__':
    main()
