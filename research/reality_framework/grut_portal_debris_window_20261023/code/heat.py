#!/usr/bin/env python3
"""Old mass point: lifetime/recoil conflict and the weaker CMB heat screen."""
import importlib.util
import json
import numpy as np
from scipy.optimize import brentq
from widths import ROOT, Widths, HBAR, visible_coefficient

spec = importlib.util.spec_from_file_location('previous_debris_transport', ROOT.parent/'grut_fragmentation_transport_20261022/code/transport.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def main():
    scan = json.loads((ROOT/'SCAN_RESULTS.json').read_text())
    ref = next(r for r in scan['selected'] if r['Q']==39944 and r['lambda_phi']==.01 and r['v_phi_GeV']==1.)
    w = Widths()
    width = w.unit(np.sqrt(.02))
    coeff = visible_coefficient(1., .01, width)
    budget = ref['scalar_release_GeV_per_entropy']
    af = old.T0/.02*(old.GS0/10.75)**(1/3)

    def evaluate(kappa, gamma=1., points=16000):
        tau = HBAR/(coeff*kappa*kappa)
        r, _ = old.kernel(budget, af, gamma, tau, points)
        return dict(kappa=kappa, gamma_initial=gamma, scalar_lifetime_seconds=tau,
                    heat_proxy=r['pre_recombination_visible_heat_proxy'],
                    post_z1000_emitted_energy_over_photon_energy=r['post_z1000_emitted_energy_over_photon_energy'])

    # Restricted to the fast-decay branch of the heat condition. Arbitrarily
    # long lifetimes can also lower the pre-recombination integral, but then
    # post-recombination constraints must be calculated, not called safe.
    kheat = float(np.exp(brentq(lambda x: np.log(evaluate(float(np.exp(x)))['heat_proxy']/6e-5),
                               np.log(1e-9), np.log(1e-6))))
    rows = [evaluate(k) for k in (1e-10, 5.007066814716605e-10,
                                  ref['kappa_LZ_shape_matched_ESTIMATE'],
                                  ref['kappa_LZ_same_mass_count_ESTIMATE'], kheat)]
    rows.append(evaluate(6.209329563372114e-10, 1000.))
    k1 = ref['kappa_for_lifetime_1s_central']
    result = dict(scope='Original scalar mass and cold all-s source only; heat-deposition proxy, not a BBN calculation or all-spectra no-go.',
                  kappa_heat_fast_branch=kheat, kappa_for_1s=k1,
                  ratio_kappa_1s_to_LZ_shape_screen=k1/ref['kappa_LZ_shape_matched_ESTIMATE'],
                  ratio_event_rate_1s_to_LZ_shape_screen=(k1/ref['kappa_LZ_shape_matched_ESTIMATE'])**2,
                  accepted_signals_per_tonne_year_at_kappa_for_1s=ref['xenonnt20ty_accepted_at_kappa_1e9']/20*(k1/1e-9)**2,
                  recoil_count_times_lifetime_at_kappa_1e9=ref['xenonnt20ty_accepted_at_kappa_1e9']*HBAR/(coeff*1e-18),
                  rows=rows)
    (ROOT/'HEAT_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
