"""cf_state common constants and helpers (own code, 27 Oct 2026).

SI constants from scipy.constants (CODATA 2018 as shipped with scipy 1.18).
Every quoted experimental number carries its source in QUOTED below; read levels are in NOTES.md.
"""
import json
import numpy as np
from scipy import constants as sc

hbar = sc.hbar
c = sc.c
G = sc.G
kB = sc.k
u = sc.physical_constants['atomic mass constant'][0]
g_earth = 9.80665
t_P = np.sqrt(hbar * G / c**5)          # Planck time
l_P = np.sqrt(hbar * G / c**3)
m_P = np.sqrt(hbar * c / G)
GeV = 1e9 * sc.e                         # J
m0_CSL = u                               # CSL reference mass (atomic mass unit)

# Quoted experimental/theory numbers [Q], with sources (see NOTES.md for read level).
QUOTED = {
    # collapse models
    'CSL_GRW_lambda': (1e-16, 's^-1 at r_C=1e-7 m', 'GRW 1986; Carlesso et al. arXiv:2203.04231 (Nat. Phys. 18, 243 (2022))'),
    'CSL_Adler_lambda': (4e-8, 's^-1 at r_C=1e-7 m (x10^{+-2})', 'Adler 2007 via arXiv:2203.04231'),
    'CSL_bound_LPF': (3.8e-9, 's^-1 at r_C=1e-7 m', 'LISA Pathfinder, via arXiv:2203.04231'),
    'CSL_bound_cantilever': (2.0e-10, 's^-1 at r_C=1e-7 m', 'Vinante et al. PRL 125, 100404 (2020), via arXiv:2203.04231'),
    'CSL_bound_bulkheat': (3.3e-11, 's^-1 at r_C=1e-7 m', 'Adler & Vinante PRA 97, 052119 (2018), via arXiv:2203.04231'),
    'CSL_bound_Xray_GranSasso': (5.2e-13, 's^-1 at r_C=1e-7 m (white noise)', 'Donadi et al. EPJC 81, 773 (2021), via arXiv:2203.04231'),
    'CSL_bound_Majorana_mp_nuclei': (0.494, 'lambda/r_C^2 in s^-1 m^-2 (m-p CSL, coherent nuclei, white noise)', 'Arnquist et al. (Majorana) PRL 129, 080401 (2022), arXiv:2202.01343'),
    'CSL_bound_Majorana_mp_electrons': (17.4, 'lambda/r_C^2 in s^-1 m^-2 (m-p CSL, quasi-free electrons)', 'arXiv:2202.01343'),
    'CSL_bound_atomint_2026': (1.27e-5, 's^-1 at r_C=1e-5 m', 'Zhu et al. arXiv:2605.26860 (2026)'),
    'DP_R0_bound_Majorana': (2.54e-10, 'm, 95% CL lower bound', 'arXiv:2202.01343 with erratum PRL 130, 239902 (2023)'),
    'DP_R0_bound_GranSasso': (0.54e-10, 'm', 'Donadi et al. Nat. Phys. 17, 74 (2021), arXiv:2111.13490'),
    'DP_R0_upper_effective': (1e-4, 'm, "effective collapse" upper bound (somewhat arbitrary)', 'Figurato et al. NJP 26, 113004 (2024), arXiv:2406.18494'),
    # clocks and interferometers
    'Sr_clock_freq': (429228004229873.0, 'Hz (87Sr 1S0-3P0)', 'BIPM recommended value [R]'),
    'Bothwell_coherence': (37.0, 's Ramsey contrast decay constant', 'Bothwell et al. Nature 602, 420 (2022), arXiv:2109.12238'),
    'Bothwell_extra_decoherence_time': (52.8, 's beyond 118 s natural lifetime', 'arXiv:2109.12238'),
    'Bothwell_frac_uncertainty': (7.6e-21, 'fractional frequency measurement uncertainty', 'arXiv:2109.12238'),
    'Young_tweezer_coherence': (48.0, 's single-atom coherence', 'Young et al. Nature 588, 408 (2020), arXiv:2004.06095'),
    'Fein_mass_u': (25000.0, 'u (beyond 25 kDa), up to 2000 atoms', 'Fein et al. Nat. Phys. 15, 1242 (2019) via edoc.unibas.ch/75500'),
    'Fein_lambda_dB': (53e-15, 'm', 'same'),
    'Fein_length': (2.0, 'm Talbot-Lau interferometer', 'same'),
    'Fein_visibility_fraction': (0.90, '> 90% of expected visibility', 'same'),
    'Fein_grating_period': (266e-9, 'm (LUMI grating period)', '[R] not re-read'),
    'Bose_mass': (1e-14, 'kg', 'Bose et al. PRL 119, 240401 (2017), arXiv:1707.06050'),
    'Bose_radius': (1e-6, 'm', 'same'),
    'Bose_dx': (250e-6, 'm superposition size', 'same'),
    'Bose_tau': (2.5, 's', 'same'),
}


def dump(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, complex):
            return [o.real, o.imag]
        raise TypeError(type(o))
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
