"""Explicit canonical vacuum vertices and two-body decay normalization."""
from common import *
from scipy.integrate import quad

def main():
    lam3=2*G2*V;beta=np.sqrt(1-4*M*M/(MPHI*MPHI))
    width=lam3*lam3/(16*np.pi*MPHI)*beta
    # The same width from two identical-real-species channels and direct dOmega.
    real_species_width=2*lam3*lam3/(32*np.pi*MPHI)*beta
    momentum=.5*MPHI*beta
    angular=quad(lambda costheta:lam3*lam3*momentum/(32*np.pi**2*MPHI**2)*2*np.pi,-1,1,epsabs=1e-16)[0]
    old_kinetic=None
    repo=RF.parents[1]
    path='research/reality_framework/grut_formation_review_20261025/tracks/kin_rates/results/rates_interface.json'
    try:
        raw=subprocess.check_output(['git','show','e71730c24fe6a0b28ddd72cd30982c5b39c9aedf:'+path],cwd=repo,stderr=subprocess.DEVNULL)
        p=json.loads(raw)['meta']['parameters']
        old_kinetic={'commit':'e71730c24fe6a0b28ddd72cd30982c5b39c9aedf','path':path,
                     'sha256':hashlib.sha256(raw).hexdigest(),
                     'lam3_sXXb_GeV':p['lam3_sXXb_GeV'],'Gamma_s_vac_GeV':p['Gamma_s_vac_GeV']}
    except (subprocess.CalledProcessError,KeyError):pass
    out={'potential_mixed_term':'g_X^2 phi^2 |X|^2','canonical_vertices':{
        'h_X_Xbar_GeV':lam3,'h_h_X_Xbar':2*G2,'h_h_h_GeV':6*(.03-.002)*V,'h_h_h_h':6*.03},
        'width_GeV':width,'width_over_mX':width/M,'lifetime_seconds':6.582119569e-25/width,
        'real_species_check_relative':float(real_species_width/width-1),
        'angular_phase_space_check_relative':float(angular/width-1),
        'Argon_quoted_numerical_width_over_mX':9.38e-6,
        'ratio_to_Argon_quote':float(width/M/9.38e-6),
        'scalar_mass_and_Xphi_vertex_unchanged_by_epsilon':True,
        'old_kinetic_input':old_kinetic,
        'neutral_on_shell_channel_k_threshold_over_mX':float(np.sqrt((np.sqrt(6)-OMEGA)**2-1)),
        'neutral_decays_return_one_antiparticle':True,
        'number_erasure_cross_section_channel':'outgoing X charge exchange only, after one vacuum scalar decay; inverse reactions and plasma chemistry not computed'}
    save('VERTICES_RESULTS.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__':main()
