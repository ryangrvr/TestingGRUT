"""Leading-portal elastic forward model; not an experimental likelihood."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, simpson
from scipy.interpolate import PchipInterpolator
from scipy.special import erf, spherical_jn

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT.parent / "grut_unified_working_theory_20261022"
G, LAM, KAPPA = .1, .01, 1e-10
MH, VH, MN, FN = 125., 246., .939, .308
A, MA = 131, 131 * .9315
C_KMS, C_CMS = 299792.458, 29979245800.
HBARC, GEV2_CM2 = .1973269804, .3893793721e-27
V0, VE, VESC = 220., 232., 544.
YEAR = 365.25 * 86400.
N_TARGET = 1e6 / (A * 1.66054e-24)
RHO = .3
NZ = erf(VESC/V0) - 2*(VESC/V0)*np.exp(-(VESC/V0)**2)/np.sqrt(np.pi)


def speed_pdf(v):
    if v <= 0 or v >= VESC + VE:
        return 0.
    return v/(NZ*np.sqrt(np.pi)*V0*VE) * (
        np.exp(-((v-VE)/V0)**2)
        - np.exp(-min((v+VE)/V0, VESC/V0)**2))


def eta(vmin):
    x, y, z = np.asarray(vmin)/V0, VE/V0, VESC/V0
    first = erf(x+y)-erf(x-y)-4*y*np.exp(-z*z)/np.sqrt(np.pi)
    second = erf(z)-erf(x-y)-2*(z+y-x)*np.exp(-z*z)/np.sqrt(np.pi)
    return np.where(x < z-y, first, np.where(x < z+y, second, 0.))/(2*NZ*VE)


def helm(q):
    skin = .9
    r1 = np.sqrt((1.2*A**(1/3))**2-5*skin**2)
    x = np.asarray(q)*r1/HBARC
    return np.where(x == 0, 1., 3*spherical_jn(1, x)/np.where(x == 0, 1., x)) * np.exp(-.5*(np.asarray(q)*skin/HBARC)**2)


def ff(rho, s, qhat):
    density = rho*rho*(1-s*s)
    norm = simpson(density, x=rho)
    out = [simpson(density*np.sinc(q*rho/np.pi), x=rho)/norm for q in np.atleast_1d(qhat)]
    return np.asarray(out)


def parameters(mass, scale):
    mu_n = mass*MN/(mass+MN)
    mu_a = mass*MA/(mass+MA)
    gq = -KAPPA*VH*mass/(2*LAM*scale*scale)
    gn = FN*MN/VH
    sigma_n = mu_n**2*(gq*gn/MH**2)**2/np.pi*GEV2_CM2
    sigma_a = sigma_n*A*A*(mu_a/mu_n)**2
    # Upper bound from concentrating the positive source at one point.
    born_kernel_bound = mu_a*A*abs(gq*gn)/(2*np.pi*MH)
    return mu_n, mu_a, gq, gn, sigma_n, sigma_a, born_kernel_bound


def run():
    (ROOT / "data").mkdir(exist_ok=True)
    old = json.loads((PARENT / "RESULTS.json").read_text())
    energy = np.linspace(5., 40., 1401)
    q = np.sqrt(2*MA*energy*1e-6)
    output = {"status": "Exploratory, unadopted portal proposal; ideal rates, not accepted detector counts or limits.",
              "parent": "7abab860441ec6711a4d7fc4cc01c48c86073b23",
              "inputs": {"g_X": G, "lambda_phi": LAM, "kappa": KAPPA,
                         "m_h_GeV": MH, "v_H_GeV": VH, "m_N_GeV": MN,
                         "f_N": FN, "f_N_reference_uncertainty": .018,
                         "A": A, "m_A_GeV": MA,
                         "rho_local_GeV_cm3": RHO, "v0_kms": V0,
                         "vEarth_kms": VE, "vescape_kms": VESC,
                         "efficiency": 1., "background": 0.,
                         "window_keV": [5., 40.], "halo_fraction": 1.},
              "halo": {"normalization_quad": quad(speed_pdf, 0, VESC+VE, points=[VESC-VE], epsabs=1e-12)[0]},
              "profiles": [], "examples": []}
    thresholds = [0., 50., 100., 200., 400., 600., 775.]
    output["halo"]["eta_checks"] = [
        {"vmin_kms": x, "analytic": float(eta(x)),
         "quad": quad(lambda v: speed_pdf(v)/v if v else 0., x, VESC+VE,
                      points=[VESC-VE] if x < VESC-VE else None, epsabs=1e-13)[0]}
        for x in thresholds]
    for old_row in old["profiles"]:
        charge = old_row["Q"]
        if charge not in (19972, 39944, 79988):
            continue
        source = PARENT / "data" / f"profile_Q{charge}.csv"
        raw = np.loadtxt(source, delimiter=",", skiprows=1)
        rho, s = raw[:, 0], raw[:, 3]
        d2 = simpson(rho**2*(1-s*s), x=rho)
        d4 = simpson(rho**4*(1-s*s), x=rho)
        ebar_integral = 4*np.pi*LAM/G**2*d2
        probe = ff(rho, s, [0., 1e-4, .25, .5, 1., 1.5])
        coarse = ff(rho[::2], s[::2], np.linspace(0, 2, 81))
        full = ff(rho, s, np.linspace(0, 2, 81))
        spline = PchipInterpolator(rho, 1-s*s)
        dquad = quad(lambda x: x*x*spline(x), rho[0], rho[-1],
                     points=rho, epsabs=1e-8, limit=1000)[0]
        qquad = quad(lambda x: x*x*spline(x)*np.sinc(x/np.pi), rho[0], rho[-1],
                     points=rho, epsabs=1e-8, limit=1000)[0]/dquad
        radial_extra = 4*np.pi*LAM/G**3*simpson(rho*rho*(s-1)**2, x=rho)
        radial_source = 4*np.pi/G**3*simpson(
            rho*rho*(G*G*s*raw[:, 1]**2+LAM*(3*(s-1)**2+(s-1)**3)), x=rho)
        radial_from_field = -4*np.pi*2*LAM/G**3*simpson(rho*rho*(s-1), x=rho)
        output["profiles"].append({
            "Q": charge, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "Ebar_parent": old_row["Ebar"], "Ebar_profile_integral": ebar_integral,
            "normalization_relative_residual": ebar_integral/old_row["Ebar"]-1,
            "depletion_mean_rho2": d4/d2, "F_at_qhat": probe.tolist(),
            "qhat_values": [0., 1e-4, .25, .5, 1., 1.5],
            "small_q_moment_residual": float((1-probe[1])/(1e-8*d4/d2/6)-1),
            "decimation_max_absolute_difference": float(np.max(abs(full-coarse))),
            "qhat1_Pchip_quad_difference": float(qquad-probe[4]),
            "radial_static_source_zero": float(radial_source),
            "radial_source_EOM_integral_relative": float(radial_source/radial_from_field-1),
            "radial_source_extra_over_E_over_v": float(radial_extra/(old_row["Ebar"]/G)),
            "radial_source_identity_relative": float((radial_source-old_row["Ebar"]/G-radial_extra)/(old_row["Ebar"]/G))})
        for scale in (1., 10., 1000., 2.435e13):
            mass = old_row["mass_over_gv"]*G*scale
            mu_n, mu_a, gq, gn, sigma_n, sigma_a, ib = parameters(mass, scale)
            fobj = ff(rho, s, q/(G*scale))
            vmin = C_KMS*np.sqrt(MA*energy*1e-6/(2*mu_a*mu_a))
            response = helm(q)**2 * (MH**2/(MH**2+q*q))**2
            # dR/dE: nuclei/tonne * n_DM * c * d-sigma prefactor * <1/beta>.
            point = (N_TARGET*RHO/mass*C_CMS*MA*sigma_a/(2*mu_a**2)
                     * 1e-6*C_KMS*eta(vmin)*YEAR*response)
            extended = point*fobj*fobj
            npoint, nextended = simpson(point, x=energy), simpson(extended, x=energy)
            mphi2 = 2*LAM*scale*scale
            mix_theta = KAPPA*scale*VH/(MH*MH-mphi2)
            # Positive-semidefinite quartic condition for both vacuum displacements.
            stability_ratio = KAPPA**2/(4*LAM*(MH*MH/(2*VH*VH)))
            mass_rule = -KAPPA*VH*mass/(2*LAM*scale**2)
            integral_rule = -KAPPA*VH/2*(4*np.pi/(G**3*scale))*d2
            row = {"Q": charge, "v_phi_GeV": scale, "mass_GeV": mass,
                   "m_phi_GeV": math.sqrt(mphi2),
                   "R90_fm": old_row["R90_over_inverse_gv"]/(G*scale)*HBARC,
                   "g_Qh_zero": gq, "g_Nh": gn,
                   "source_integral_vs_mass_relative": integral_rule/mass_rule-1,
                   "sigma_n_zero_cm2": sigma_n, "sigma_xe_zero_cm2": sigma_a,
                   "Born_kernel_sufficient_upper_bound": ib,
                   "small_mixing_theta": mix_theta,
                   "quartic_stability_ratio": stability_ratio,
                   "ideal_events_tonne_year_5_40keV": float(nextended),
                   "formal_kappa_for_one_ideal_event_per_tonne_year": float(KAPPA/np.sqrt(nextended)),
                   "point_object_ideal_events_tonne_year_5_40keV": float(npoint),
                   "extended_to_point_window_ratio": float(nextended/npoint),
                   "F2_at_ER_5_10_30_40keV": (ff(rho, s, np.sqrt(2*MA*np.asarray([5,10,30,40])*1e-6)/(G*scale))**2).tolist()}
            # Fixed-v scaling and an independent uniform-Higgs derivative.
            dh = .001
            vp = np.sqrt(scale**2-KAPPA/LAM*VH*dh)
            vm = np.sqrt(scale**2+KAPPA/LAM*VH*dh)
            # Stable algebra avoids subtracting nearly identical large masses.
            derivative = mass/scale*(-2*KAPPA/LAM*VH)/(vp+vm)/2
            row["uniform_Higgs_derivative_vs_source_relative"] = derivative/gq-1
            output["examples"].append(row)
            if scale == 1. and charge == 39944:
                np.savetxt(ROOT/"data"/f"recoil_Q{charge}_v{int(scale)}.csv",
                           np.column_stack([energy, q, fobj, helm(q), point, extended]),
                           delimiter=",", header="E_keV,q_GeV,F_object,F_Helm,point_events_per_tonne_year_keV,extended_events_per_tonne_year_keV", comments="")
    # Cross-check equal-velocity total contact scattering by direct recoil integration.
    mass = output["examples"][0]["mass_GeV"]
    mun, mua, _, _, sn, sa, _ = parameters(mass, 1.)
    beta = 220./C_KMS
    emax = 2*mua*mua*beta*beta/MA
    output["contact_cross_section_recoil_integral_ratio"] = quad(
        lambda e: MA*sa/(2*mua*mua*beta*beta), 0, emax, epsabs=1e-60)[0]/sa
    (ROOT/"RESULTS.json").write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"profiles": len(output["profiles"]), "examples": len(output["examples"]),
                      "Q39944_examples": [x for x in output["examples"] if x["Q"] == 39944]}, indent=2))


if __name__ == "__main__":
    run()
