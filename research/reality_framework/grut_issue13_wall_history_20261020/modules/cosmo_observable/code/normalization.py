"""Part 3: tension normalization reconciliation. SI surface density <-> GeV^3, the Lazanu-Martins-Shellard (LMS)
arXiv:1505.03673 CMB bound, the Zel'dovich-Kobzarev-Okun (ZKO) estimate as stated in Saikawa arXiv:1703.02576 eqs.
(2.21)-(2.22), and Codex's dimensional estimate G sigma / H0 ~ 1e-5. Writes results/normalization.json.

These are SUBDOMINANT-network bounds (LMS adopt a LambdaCDM background and treat walls as first-order perturbations,
arXiv:1505.03673 sec. 2). They are used here only to ask how small sigma must be for walls to be tolerable; they are
not applied to the wall-dominated benchmark history.
"""
import json
import math
import os
import core as C

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

def cube_root_MeV(sig):
    return sig ** (1.0 / 3.0) * 1e3

def main():
    out = {'module': 'cosmo_observable', 'part': '3 normalization'}
    conv = {
        'hbar_c_GeV_m': C.HBARC_GEV_M, 'GeV_J': C.GEV_J, 'c_m_s': C.C_M_S,
        '1_GeV3_in_J_per_m2': C.GEV3_J_M2,
        '1_GeV3_in_kg_per_m2': C.GEV3_KG_M2,
        '1_kg_per_m2_in_GeV3': 1.0 / C.GEV3_KG_M2,
        'formula': '1 GeV^3 = (1 GeV)/(hbar c / 1 GeV)^2 energy per area; divide by c^2 for mass per area',
    }
    out['conversion'] = conv
    sb = C.SIGMA_BENCH
    out['benchmark'] = {'sigma_GeV3': sb, 'sigma_kg_m2': sb * C.GEV3_KG_M2, 'sigma_J_m2': sb * C.GEV3_J_M2,
                        'sigma_cuberoot_GeV': sb ** (1 / 3.)}

    t0_planck = C.AGE_PLANCK_GYR * C.GYR_S / C.HBAR_GEV_S      # GeV^-1 (retrospective LCDM age, for the estimates)
    H70 = C.kmsmpc_to_gev(70.0)
    lms = {}
    for label, sig_si, mev_label in (('R+M', 4.22e-9, 0.96), ('R+M+Lambda', 3.85e-9, 0.93)):
        lit = sig_si / C.GEV3_KG_M2                      # literal reading, kg/m^2
        cm2 = sig_si * 1e4 / C.GEV3_KG_M2                # reading the printed unit as kg/cm^2 (= 1e4 kg/m^2)
        lab = (mev_label * 1e-3) ** 3                    # label read as sigma^(1/3)
        lab_vev = C.sigma_tree(mev_label * 1e-3)         # label read as a vev eta at lam = lam_phi = 0.01
        d = {
            'printed_sigma_kg_m2': sig_si, 'printed_label_MeV': mev_label,
            'literal_sigma_GeV3': lit, 'literal_cuberoot_MeV': cube_root_MeV(lit),
            'kg_per_cm2_reading_sigma_GeV3': cm2, 'kg_per_cm2_reading_cuberoot_MeV': cube_root_MeV(cm2),
            'label_as_cuberoot_sigma_GeV3': lab,
            'label_as_vev_at_lam0.01_sigma_GeV3': lab_vev,
            'ratio_label_cuberoot_to_literal_cuberoot': mev_label / cube_root_MeV(lit),
            'ratio_kgcm2_cuberoot_to_label': cube_root_MeV(cm2) / mev_label,
            # LMS eq.(28) amplitude coefficient 6 pi G t0 sigma / c (dimensionless), and the matter-era wall fraction
            # rho_w / rho_bar = (1 - lambda)(A eta/V) 6 pi G sigma t0 with (1-lambda)(A eta/V) = 0.643
            'G_sigma_t0_literal': C.G_N * lit * t0_planck,
            'G_sigma_t0_kgcm2': C.G_N * cm2 * t0_planck,
            '6piGsigma_t0_literal': 6 * math.pi * C.G_N * lit * t0_planck,
            '6piGsigma_t0_kgcm2': 6 * math.pi * C.G_N * cm2 * t0_planck,
            'rho_w_over_rho_bar_today_literal': C.A_MATTER_LMS * 6 * math.pi * C.G_N * lit * t0_planck,
            'rho_w_over_rho_bar_today_kgcm2': C.A_MATTER_LMS * 6 * math.pi * C.G_N * cm2 * t0_planck,
        }
        lms[label] = d
    out['LMS_1505.03673'] = lms
    out['LMS_reconciliation'] = (
        'The printed SI values (4.22e-9, 3.85e-9 kg/m^2) convert to sigma^(1/3) = %.1f, %.1f keV, a factor %.1f below '
        'the printed MeV labels; sigma differs by ~1e4 = (m/cm)^2. Reading the unit as kg/cm^2 gives sigma^(1/3) = '
        '%.3f, %.3f MeV, within %.1f%% of both labels, and gives G sigma t0 = %.2e, consistent with the paper\'s stated '
        'agreement with the deltaT/T <~ 1e-5 estimate, whereas the literal reading gives G sigma t0 = %.1e. We therefore '
        'adopt sigma_95 = %.3e GeV^3 (R+M) and %.3e GeV^3 (R+M+Lambda) as the reconciled bound and identify the MeV label '
        'with sigma^(1/3) in natural units (not a vev). This is our reconciliation, not an erratum by the authors; '
        'the literal reading is carried as a stricter alternative.' % (
            lms['R+M']['literal_cuberoot_MeV'] * 1e3, lms['R+M+Lambda']['literal_cuberoot_MeV'] * 1e3,
            lms['R+M']['ratio_label_cuberoot_to_literal_cuberoot'],
            lms['R+M']['kg_per_cm2_reading_cuberoot_MeV'], lms['R+M+Lambda']['kg_per_cm2_reading_cuberoot_MeV'],
            100 * max(abs(lms['R+M']['ratio_kgcm2_cuberoot_to_label'] - 1), abs(lms['R+M+Lambda']['ratio_kgcm2_cuberoot_to_label'] - 1)),
            lms['R+M']['G_sigma_t0_kgcm2'], lms['R+M']['G_sigma_t0_literal'],
            lms['R+M']['kg_per_cm2_reading_sigma_GeV3'], lms['R+M+Lambda']['kg_per_cm2_reading_sigma_GeV3']))

    # ZKO estimate in the form of Saikawa 1703.02576 eq. (2.21): delta rho / rho ~ G sigma t0 <~ 1e-5
    zko = 1e-5 / (C.G_N * t0_planck)
    out['ZKO_Saikawa_eq2.21'] = {'sigma_max_GeV3': zko, 'cuberoot_MeV': cube_root_MeV(zko),
                                 'note': 'order-of-magnitude; t0 = 13.797 Gyr; ZKO 1974 itself not read (secondary via 1703.02576 sec. 2.2)'}
    # Codex: G sigma / H0 ~ 1e-5 with H0 = 70 km/s/Mpc
    cdx = 1e-5 * H70 / C.G_N
    out['Codex_G_sigma_over_H0'] = {'sigma_GeV3': cdx, 'cuberoot_MeV': cube_root_MeV(cdx),
                                    'benchmark_over_this': sb / cdx, 'quoted_value': 2.23e-9, 'quoted_ratio': 4.24e16}
    # Saikawa's coefficient check: G sigma t0 = 1e12 (sigma/TeV^3)
    out['Saikawa_eq2.21_coefficient_check'] = {'G_t0_TeV3': C.G_N * t0_planck * 1e9, 'quoted': 1e12}

    bounds = {
        'LMS_R+M+Lambda_reconciled': lms['R+M+Lambda']['kg_per_cm2_reading_sigma_GeV3'],
        'LMS_R+M_reconciled': lms['R+M']['kg_per_cm2_reading_sigma_GeV3'],
        'LMS_R+M+Lambda_literal_SI': lms['R+M+Lambda']['literal_sigma_GeV3'],
        'LMS_R+M_literal_SI': lms['R+M']['literal_sigma_GeV3'],
        'ZKO_Saikawa_Gsigma_t0_1e-5': zko,
        'Codex_Gsigma_over_H0_1e-5': cdx,
    }
    tab = {}
    for k, s in bounds.items():
        v = C.v_from_sigma(s)
        tab[k] = {'sigma_GeV3': s, 'sigma_kg_m2': s * C.GEV3_KG_M2, 'cuberoot_MeV': cube_root_MeV(s),
                  'benchmark_over_bound': sb / s, 'v_max_GeV_at_lam0.01': v, 'v_max_MeV': v * 1e3,
                  'm_X_max_MeV_gX0.1': C.G_X * v * 1e3,
                  'r_max_from_mX_r<63.93GeV': 63.93 / (C.G_X * v), 'r_max_from_mX_r<74.64GeV_stress': 74.64 / (C.G_X * v)}
    out['bounds_table'] = tab
    out['v_max_general_formula'] = 'v_max = (3 sigma_max / (2 sqrt(2 lambda_phi)))^(1/3); at lambda_phi=0.01 v_max = 2.1972 sigma_max^(1/3)'
    out['vilenkin_gravitational_time_benchmark'] = {
        't_sigma_s': 1.0 / (2 * math.pi * C.G_N * sb) * C.HBAR_GEV_S,
        'formula': 't_sigma = 1/(2 pi G sigma) (planar-wall repulsion scale, Vilenkin 1983 / Ipser-Sikivie 1984; secondary)'}
    with open(os.path.join(RES, 'normalization.json'), 'w') as f:
        json.dump(out, f, indent=1)
    return out

if __name__ == '__main__':
    o = main()
    print(json.dumps(o, indent=1))
