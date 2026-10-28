"""Conditional rate/expansion diagnostics, with every population assumption explicit."""
from common import *

def main():
    part=json.loads((ROOT/'PARTIAL_RESULTS.json').read_text())
    resonance=ROOT/'RESONANCE_CHECKS.json'
    if resonance.exists():
        rate=json.loads(resonance.read_text())['rows'][0]['controls'][-1]['resolved_average']
    else:rate=part['source_averages'][0]['sigma_v_X_covered_mXminus2']
    sv=rate/M**2
    YQ=POINT['final_population_matching_hypothesis']['YQ_target'];Mred=2.435e18
    rows=[]
    for T,g in [(1.,100.),(.02,10.75),(1e-3,10.75),(1e-6,3.36)]:
        gs=g if T>1e-5 else 3.91
        entropy=2*np.pi**2/45*gs*T**3;H=np.sqrt(np.pi**2*g/90)*T*T/Mred
        nobj=YQ*entropy/QREF
        rows.append({'visible_temperature_GeV':T,'gstar':g,'gstar_entropy':gs,
                     'assumed_eta':1.,'assumed_monocharged_objects_Q':QREF,
                     'Gamma_charge_exchange_over_H_with_frozen_source_spectrum':nobj*sv/H})
    examples=[]
    for F in (.5,1.):
        from scipy.optimize import brentq
        S=brentq(lambda s:B*(s*s-1)*s-3*C*(s-1)**2+F*F*s,0,1)
        om=np.sqrt(S*S+A*F*F);nQ=om*F*F/G2*M**3
        g=100.;T=(45/(2*np.pi**2*g)*nQ/YQ)**(1/3)
        H=np.sqrt(np.pi**2*g/90)*T*T/Mred
        examples.append({'supplied_F':F,'illustrative_gstar':g,'visible_T_from_supplied_density_and_YQ_GeV':T,
                         'H_over_mX':H/M,'scale_factor_change_over_300_mX_inverse':300*H/M})
    source=json.loads((ROOT/'SOURCE_RESULTS.json').read_text())['rows'];charge_trade=[]
    for r in source:
        total=(r['Nplus_out']+r['Nminus_out']+2*r['Nphi_out_net_vacuum_reference'])/QREF
        charge_trade.append({'pulse':r['tag'],'measured_Nplus_plus_Nminus_plus_two_Nphi_over_Q':total,
            'eta_if_all_antiparticles_convert_once_and_no_free_X_rejoins':1-total,
            'assumptions':'single object source; phi decays once, all emitted antiX -> X, no further sources, subsequent X capture omitted; not a formed population'})
    save('CONDITIONAL_CLOCKS.json',{'source_averaged_sigma_v_GeVminus2_covered':sv,'eta1_rate_controls':rows,
        'local_expansion_examples':examples,'charge_trade_controls':charge_trade,
        'actual_capture_history':None,'actual_reheating_background_used_here':False,
        'caution':'Sources, momentum redshift, evolving nobj/Q, repeated scattering and inverse reactions must be evolved. Large conditional Gamma/H is not a final erasure verdict.'})
    print(json.dumps({'sigma_v_GeVminus2_covered':sv,'rate_controls':rows,'expansion_examples':examples,'charge_trade':charge_trade},indent=2))

if __name__=='__main__':main()
