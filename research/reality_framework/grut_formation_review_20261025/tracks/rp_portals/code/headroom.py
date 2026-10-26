"""rp_portals / headroom.py (Claude, own code; thresholds and existing debris quoted from our own tracks).

Item 5 of the task: the raw-source window between detectability and hidden-energy injection at the heavy point.
For each case-L threshold variant R_L on record [Q own: hw_critic/results/thresholds.json] and each inflationary-vector
estimate R_A at e = 1e-7 (Codex band 40.5626 [Q], own modes 48.47 [Q hw_gates], GMR formula 64.0 [Q hw_critic]),
headroom = R_L - R_med - R_A (R_med = 1.1437, mediator debris [Q own]).  The portal sources must satisfy
R_thermal(eps_BX) <= headroom (momentum-resolved rows of freezein.json, R ∝ eps^2 at small eps) and
R_h(b) = R_h(5e-11) (b/5e-11)^2 <= headroom (higgs_source.json).  Raw-source statements only: no f_sym, no kinetics.
Also the rest-mass-only statements (f_rest <= 0.01) and the detector-level couplings of the earlier proposal [Q 8b0d96b].

Usage: python3 -B headroom.py <freezein.json> <higgs_source.json> <out.json>
"""
import json, math, sys, pathlib
import numpy as np

def main(fz_path, hg_path, out):
    fz = json.load(open(fz_path)); hg = json.load(open(hg_path))
    root = pathlib.Path(__file__).resolve().parents[2]
    thr = json.load(open(root/'hw_critic/results/thresholds.json'))['thresholds']
    R_med = 1.1437
    R_A = {'codex_band_40.56': 40.56262288525944, 'own_modes_48.47': 48.47, 'GMR_formula_64.0': 64.0}
    # R_thermal(eps): momentum-resolved rows, log-log interpolation (monotone up to eps ~ 1e-7)
    rows = [r for r in fz['eps_rows'] if r['eps_BX'] <= 1e-7]
    le = np.log([r['eps_BX'] for r in rows]); lR = np.log([r['momentum']['R_thermal'] for r in rows])
    lf = np.log([r['f_rest_unprocessed'] for r in rows])
    def eps_for_R(R):
        if R <= 0: return 0.0
        if R > math.exp(lR[-1]): return None
        return float(math.exp(np.interp(math.log(R), lR, le)))
    def eps_for_f(f):
        return float(math.exp(np.interp(math.log(f), lf, le)))
    Rh0 = hg['R_h_gstar106p75']; b0 = 5e-11
    table = []
    for tname, RL in thr.items():
        for aname, RA in R_A.items():
            head = RL - R_med - RA
            table.append(dict(threshold=tname, R_L=RL, vector=aname, R_A=RA, headroom=head,
                              eps_BX_max=eps_for_R(head), b_Hphi_max=(b0*math.sqrt(head/Rh0) if head > 0 else 0.0)))
    # small-eps scaling constants
    r10 = [r for r in fz['eps_rows'] if r['eps_BX'] == 1e-10][0]
    kR = r10['momentum']['R_thermal']/1e-20     # R_thermal ~ kR eps^2
    kf = r10['f_rest_unprocessed']/1e-20
    # detector-level couplings of the earlier proposal (Q = 39944, v_phi = 1 GeV object) [Q 8b0d96b REPORT]
    kappa90 = dict(XENON1T_observed=2.4973e-9, XENONnT_20ty_projection_median=5.0071e-10,
                   LZ_same_mass_estimate=7.7003e-10, LZ_shape_matched_estimate=5.9950e-10)
    # rescale to the heavy-point object (E_Q = 17064.5 GeV, v = 0.924): rate per kappa^2 ∝ E_Q/(lambda v^2)^2 at fixed halo
    # density and form factor (same lambda_phi): factor (17064.5/3639.24) (1/0.924089637)^4; [E]
    resc = (17064.537/3639.24)*(1/0.924089637)**4
    b_detect = {k: v/2 for k, v in kappa90.items()}
    b_detect_heavy = {k: v/2/math.sqrt(resc) for k, v in kappa90.items()}
    # Codex scalar reference sigma_p(0) at b = 5e-11 for the 17 TeV object [Q]
    sigma_p_b5em11 = 3.927560953586381e-47
    best = max(t['headroom'] for t in table)
    res = dict(thresholds=thr, R_med=R_med, R_A=R_A, table=table,
               scaling=dict(R_thermal_per_eps2=kR, f_rest_per_eps2=kf, R_h_per_b2=Rh0/b0**2, f_rest_h_per_b2=hg['f_rest_unprocessed_gstar106p75']/b0**2),
               rest_mass_only=dict(eps_BX_for_f_rest_0p01=eps_for_f(0.01), b_for_f_rest_0p01=hg['b_for_f_rest_equal_0p01']),
               best_case=dict(headroom=best, eps_BX_max=eps_for_R(best), b_max=b0*math.sqrt(best/Rh0)),
               detector_level=dict(kappa_90_quoted=kappa90, b_equiv_proposal_object=b_detect,
                                   rescale_factor_rate_per_kappa2_heavy_object=resc, b_equiv_heavy_object_estimate=b_detect_heavy,
                                   codex_sigma_p0_cm2_at_b5em11=sigma_p_b5em11,
                                   R_h_at_b_detect={k: Rh0*(b/b0)**2 for k, b in b_detect_heavy.items()},
                                   f_rest_at_b_detect={k: hg['f_rest_unprocessed_gstar106p75']*(b/b0)**2 for k, b in b_detect_heavy.items()}),
               # heavy-photon contact check (Codex vector_scattering, Xe131): own arithmetic
               vector_contact_check=contact_check())
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps(res, indent=1))

def contact_check():
    MX = 0.1*0.924089637375362; obj = 200000*MX*0.9233161353726057
    mA = 6e-7*40*1.3334961939833566e11
    eEM = math.sqrt(4*math.pi/137.035999084)
    mXe = 131*0.931494; mu = obj*mXe/(obj+mXe)
    G = 200000*54*1e-7*eEM/mA**2
    sigma = mu*mu/math.pi*G*G*0.3893793721e-27
    flux = 0.4/obj*8e7; atoms = 1e6/131*6.02214076e23
    return dict(sigma_cm2_eps_gamma1=sigma, codex=1.8550147325805794e-51, events_per_tonne_year=flux*atoms*sigma*365.25*86400,
                codex_events=5.046454876915115e-13)

if __name__ == '__main__':
    main(*sys.argv[1:4])
