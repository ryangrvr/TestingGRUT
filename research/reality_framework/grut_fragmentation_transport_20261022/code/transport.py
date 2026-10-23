#!/usr/bin/env python3
"""Conditional debris kernels; no formation spectrum or BBN likelihood."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT.parent/'grut_formation_wall_audit_20261022'
EMPIRICAL = ROOT.parent/'grut_portal_empirical_test_20261022'
HBAR = 6.582119569e-25  # GeV s
MPL = 2.435e18
T0 = 2.7255*8.617333262e-14  # Kelvin to GeV
GS0 = 2+(7/8)*6*(4/11)
CNU = (7/8)*(4/11)**(4/3)
H0 = 67.4*1000/3.085677581491367e22  # s^-1
RHO_GAMMA0 = np.pi**2*T0**4/15
H0_GEV = H0*HBAR
OMEGA_GAMMA = RHO_GAMMA0/(3*MPL**2*H0_GEV**2)
OMEGA_R = OMEGA_GAMMA*(1+CNU*3.046)
OMEGA_M = .315
OMEGA_L = 1-OMEGA_R-OMEGA_M
PHOTON_ENERGY_PER_ENTROPY0 = 3*T0/(2*GS0)
ENTROPY_PER_PHOTON0 = (2*np.pi**2*GS0/45)/(2*1.2020569031595943/np.pi**2)
MS = np.sqrt(.02)
MD = np.sqrt(.0101)
MX = .1


def scattering_sigma_v(p_over_ms, daughter_mass, coupling_squared, identical=True, order=64):
    """Tree COM ss->yy, invariant relative velocity 2p/E in the COM frame."""
    c, weight = np.polynomial.legendre.leggauss(order)
    p = p_over_ms*MS
    energy2 = MS**2+p**2
    shat = 4*energy2
    pf = np.sqrt(energy2-daughter_mass**2)
    t_minus_md2 = -MS**2-2*p**2+2*p*pf*c
    u_minus_md2 = -MS**2-2*p**2-2*p*pf*c
    amplitude = (2*coupling_squared+6*.01*(2*coupling_squared)/(shat-MS**2)
                 +(2*coupling_squared)**2*(1/t_minus_md2+1/u_minus_md2))
    symmetry = 2 if identical else 1
    return float(np.sqrt(1-daughter_mass**2/energy2)/(16*np.pi*shat*symmetry)
                 *np.sum(weight*amplitude**2))


def hubble(a, radiation_only=False):
    # Flat adopted background, fixed late-time radiation extrapolated at very
    # early times. Entropy change is included in a_form; e+- annihilation's H
    # correction is not needed for decays at z~1e4 and is not a BBN calculation.
    return H0*np.sqrt(OMEGA_R/a**4+(0 if radiation_only else OMEGA_M/a**3+OMEGA_L))


def kernel(budget, formation_a, gamma_i, lifetime, points=16000, radiation_only=False):
    loga = np.linspace(np.log(formation_a), 0., points)
    rec = np.log(1/1001)
    if loga[0] < rec < loga[-1]:
        loga = np.sort(np.append(loga, rec))
    a = np.exp(loga)
    h = hubble(a, radiation_only)
    p_over_m = np.sqrt(gamma_i**2-1)*formation_a/a
    gamma = np.sqrt(1+p_over_m**2)
    optical_depth = cumulative_trapezoid(1/(lifetime*gamma*h), loga, initial=0.)
    survival = np.exp(-np.minimum(optical_depth, 745))
    # A moving particle's energy is gamma*m; its lab decay probability is
    # Gamma/gamma. Product = Gamma*m, without an extra gamma in this power.
    rest_budget = budget/gamma_i
    injection = rest_budget*a/PHOTON_ENERGY_PER_ENTROPY0 * survival/(lifetime*h)
    z = 1/a-1
    visibility = np.exp(-np.minimum((np.maximum(z, 0)/1.98e6)**2.5, 745))
    before_rec = loga <= rec+1e-14
    weighted = injection*visibility
    heat = trapezoid(weighted[before_rec], loga[before_rec])
    probability_density = survival/(lifetime*gamma*h)
    integrated_probability = trapezoid(probability_density, loga)
    median_a = float(np.exp(np.interp(np.log(2), optical_depth, loga))) if optical_depth[-1] >= np.log(2) else None
    median_z = None if median_a is None else 1/median_a-1
    # Momentum-redshift energy bookkeeping, all divided by initial source B.
    decay_loss = trapezoid(survival/(lifetime*h), loga)/gamma_i
    redshift_loss = trapezoid(survival*p_over_m**2/gamma, loga)/gamma_i
    remaining = survival[-1]*gamma[-1]/gamma_i
    result = dict(formation_a=float(formation_a), formation_z=float(1/formation_a-1),
                  gamma_initial=gamma_i, lifetime_seconds=lifetime,
                  emitted_energy_GeV_per_entropy=budget, scalar_yield=rest_budget/MS,
                  median_decay_z=median_z, surviving_number_fraction_today=float(survival[-1]),
                  pre_recombination_visible_heat_proxy=float(heat),
                  maximum_source_fraction_from_historical_heat_screen=float(min(1., 6e-5/max(heat, 1e-300))),
                  post_z1000_emitted_energy_over_photon_energy=float(trapezoid(injection[loga >= rec-1e-14], loga[loga >= rec-1e-14])),
                  decay_probability_integral=float(integrated_probability),
                  probability_integral_residual=float(integrated_probability-(1-survival[-1])),
                  initial_energy_fraction_decayed=float(decay_loss),
                  initial_energy_fraction_redshifted=float(redshift_loss),
                  initial_energy_fraction_remaining=float(remaining),
                  energy_continuity_residual=float(decay_loss+redshift_loss+remaining-1))
    return result, (z, injection, weighted)


def calculate(points=16000):
    parent = json.loads((PARENT/'RESULTS.json').read_text())
    case = next(c for c in parent['cases'] if c['Q'] == 39944 and c['v_phi_GeV'] == 1)
    empirical = json.loads((EMPIRICAL/'HISTORY_RESULTS.json').read_text())
    history = next(c for c in empirical['cases'] if c['Q'] == 39944 and c['v_phi_GeV'] == 1)
    yq = case['supplied_charge_yield']
    omega_bulk = (.01*.001)**.25
    budgets = dict(reference_chosen_Q=case['binding_release_GeV_per_visible_entropy'],
                   bulk_limit_at_FIXED_charge_yield=(MX-omega_bulk)*yq)
    a20 = T0/.02*(GS0/10.75)**(1/3)
    formations = dict(T20MeV=a20, T1keV=T0/1e-6, z100000=1/100001)
    lifetimes = dict(nominal=history['cold']['lifetime_seconds'],
                     enlarged_response=history['cold_production_stress_envelope']['lifetime_seconds'])
    rows = []
    curves = {}
    for bname, budget in budgets.items():
        for fname, a_form in formations.items():
            for wname, lifetime in lifetimes.items():
                for gamma in (1., 10., 1000., 1e8):
                    result, curve = kernel(budget, a_form, gamma, lifetime, points)
                    rows.append(dict(budget=bname, formation=fname, response=wname, **result))
                    if bname == 'reference_chosen_Q' and fname == 'T20MeV' and wname == 'nominal':
                        curves[f'gamma_{gamma:g}'] = np.array(curve)
    # Existing action's ss->DD threshold. XX* cancellation is threshold-only.
    lambda_phi = .01
    h2, g2 = .01**2, .1**2
    amp_dd = 4*h2-4*h2**2/lambda_phi
    amp_xx = 4*g2-4*g2**2/lambda_phi
    sigma_dd = np.sqrt(1-MD**2/MS**2)*amp_dd**2/(64*np.pi*MS**2)
    velocity_scan = [dict(p_over_ms=p,
                         sigma_v_ss_DD_GeV_minus2=scattering_sigma_v(p, MD, h2),
                         sigma_v_ss_XXstar_GeV_minus2=scattering_sigma_v(p, MX, g2, False))
                     for p in (0., .01, .1, .3, 1.)]
    ys = budgets['reference_chosen_Q']/MS
    channel_only = []
    # Separate analytic diagnostic: constant g*=g*s in each illustrative era,
    # infinite radiation domination, cold source, no inverse or other channel.
    for temp, degrees in ((.002, 10.75), (.02, 10.75), (.2, 61.75), (35., 86.25)):
        s = 2*np.pi**2*degrees*temp**3/45
        h_gev = np.sqrt(np.pi**2*degrees/90)*temp**2/MPL
        optical = ys*s*sigma_dd/h_gev
        residual = 1/(1+optical)
        channel_only.append(dict(T_formation_GeV=temp, constant_g_star=degrees,
                                 initial_chemical_rate_over_H=float(optical),
                                 scalar_number_fraction_remaining=float(residual),
                                 final_scalar_yield=float(ys*residual),
                                 generated_stable_D_yield=float(ys*(1-residual)),
                                 stable_D_rest_energy_over_reference_DM=float(MD*ys*(1-residual)/4.4e-10)))
    prompt = []
    for z in (1e4, 1e5, 1e6):
        ratio = budgets['reference_chosen_Q']/((PHOTON_ENERGY_PER_ENTROPY0)*(1+z))
        visibility = np.exp(-(z/1.98e6)**2.5)
        prompt.append(dict(redshift=z, photon_energy_injection_ratio=float(ratio),
                           blackbody_visibility=float(visibility),
                           maximum_EM_fraction_from_heat_screen=float(min(1, 6e-5/(ratio*visibility)))))
    dr = []
    for bname, budget in budgets.items():
        for name, a_form in dict(T20MeV=a20, z100000=1/100001, z10000=1/10001).items():
            dr.append(dict(budget=bname, formation=name,
                           Delta_Neff_IF_genuinely_massless_receiver=float(budget*a_form/PHOTON_ENERGY_PER_ENTROPY0/CNU)))
    stable = []
    for gamma in (1., 10., 1000., 1e8):
        u = np.sqrt(gamma**2-1)*a20*1101
        stable.append(dict(initial_gamma=gamma, p_over_m_at_z1100=float(u),
                           w_at_z1100=float(u**2/(3*(1+u**2)))))
    f2 = case['linear_condensate']['f_squared_GeV2']
    k0 = np.sqrt(.009*f2)
    qcell = yq*(2*np.pi**2*10.75*.02**3/45)*(2*np.pi/k0)**3
    n_bulk = omega_bulk*np.sqrt(.01/.001)
    bulk_radius = (3*qcell/(4*np.pi*n_bulk))**(1/3)
    wall_width = 2/MS
    lphys = 4*np.pi/k0
    # This is a direct uniform-grid resolving requirement, not a lower bound on
    # adaptive/matched methods, and thin-wall radius is only an estimate.
    grid_axis = lphys/(wall_width/4)
    return dict(scope='Conditional NR onset, collisionless scalar kernels and channel-only threshold diagnostic; no final object or formation branching prediction.',
                adopted_background=dict(T0_GeV=T0, entropy_degrees_today=GS0, H0_per_second=H0,
                    Omega_m=OMEGA_M, Omega_r=OMEGA_R, Omega_L_calibrated=OMEGA_L,
                    entropy_per_photon_today=ENTROPY_PER_PHOTON0,
                    rho_gamma_per_entropy_today_GeV=PHOTON_ENERGY_PER_ENTROPY0,
                    kernel_high_temperature_scope='Late radiation density extrapolated; not an electron-annihilation/BBN expansion solver.'),
                frozen_charge_yield=yq, energy_budgets=budgets,
                bulk_limit_fixed_yield=dict(omega_floor_GeV=omega_bulk,
                    object_energy_GeV_per_entropy=omega_bulk*yq,
                    object_DM_fraction=omega_bulk*yq/4.4e-10,
                    free_to_object_release_GeV_per_entropy=budgets['bulk_limit_at_FIXED_charge_yield']),
                energy_per_photon_today_convention_GeV=budgets['reference_chosen_Q']*ENTROPY_PER_PHOTON0,
                maximum_single_scalar_gamma_from_one_reference_objects_binding=float(case['free_to_object_release_GeV_per_object']/MS),
                spectrum=dict(m_s_GeV=MS, m_X_GeV=MX, m_D_GeV=MD,
                    hidden_massless_free_particle_present=False,
                    gamma_needed_for_D_p_over_m_1_at_z1100=float(np.sqrt(1+(1/(1101*a20))**2)),
                    stable_D_equations_of_state=stable),
                threshold=dict(amplitude_ss_DD=float(amp_dd), amplitude_ss_XXstar=float(amp_xx),
                    sigma_v_ss_DD_GeV_minus2=float(sigma_dd), finite_COM_momentum_scan=velocity_scan,
                    channel_only=channel_only),
                prompt_EM_heat_screen=prompt, massless_receiver_counterfactual=dr,
                direct_full_field_hierarchy=dict(k0_GeV=float(k0), charge_per_fastest_wavelength_cube_DIAGNOSTIC=float(qcell),
                    estimated_bulk_core_radius_GeV_inverse=float(bulk_radius),
                    mediator_wall_width_2_over_ms_GeV_inverse=float(wall_width),
                    box_length_GeV_inverse=float(lphys),
                    uniform_grid_axis_for_4_cells_across_wall=float(grid_axis),
                    uniform_grid_cells=float(grid_axis**3),
                    six_real_fields_float64_minimum_bytes=float(6*8*grid_axis**3)),
                decay_kernels=rows), curves


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=int, default=16000)
    parser.add_argument('--output', type=Path, default=ROOT/'TRANSPORT_RESULTS.json')
    args = parser.parse_args()
    result, curves = calculate(args.points)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    with args.output.with_suffix('.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=result['decay_kernels'][0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(result['decay_kernels'])
    np.savez_compressed(args.output.with_suffix('.npz'), **curves)
    print(json.dumps({k: v for k, v in result.items() if k != 'decay_kernels'}, indent=2))
    for r in result['decay_kernels']:
        if r['budget'] == 'reference_chosen_Q' and r['formation'] == 'T20MeV':
            print(r['response'], r['gamma_initial'], 'median z', r['median_decay_z'],
                  'heat', r['pre_recombination_visible_heat_proxy'],
                  'maximum source f', r['maximum_source_fraction_from_historical_heat_screen'])
