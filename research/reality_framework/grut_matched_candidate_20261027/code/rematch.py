"""Recalibrate the actual changed-profile candidate, retaining earlier controls."""
import math,json
import numpy as np
from match import ROOT,old,run,save


def main():
    prior=json.loads((ROOT/'MATCHED_RESULTS.json').read_text())
    deformation=json.loads((ROOT/'DEFORMATION_RESULTS.json').read_text())
    epq=deformation['primary']['mass_GeV']/200000
    D=prior['target_DM_per_visible_entropy_GeV'];required=D/epq
    exe=old.build();rows=[];arrays={}
    for beta,gamma in ((3.,.5),(5.,1.)):
        tag=f'b{beta:g}_g{gamma:g}'
        base=next(x for x in prior['rows'] if x['tag']==tag+'_baseline')
        scale=(required/base['absolute_yield_magnitude'])**(2/3)
        row,arr=run(exe,'M1_'+tag+'_matched',beta,gamma,scale)
        row['target_mass_ratio']=row['absolute_yield_magnitude']*epq/D
        rows.append(row)
        keep=(arr[:,0]<=30.000001)|(np.arange(len(arr))%100==0)
        arrays[tag]=arr[keep]
        print('DEFORMED_MATCH',tag,scale,row['target_mass_ratio'],flush=True)
    central=rows[0]
    for delta in (-.001,.001,-.0005,.0005):
        row,_=run(exe,'M1_phase_'+str(delta),3.,.5,central['uniform_scale'],delta=delta,samples=8000)
        rows.append(row)
    tight,_=run(exe,'M1_tighter',3.,.5,central['uniform_scale'],tol=2e-12,samples=8000)
    rows.append(tight)
    def slope(step):
        pos=next(x for x in rows if x['tag']=='M1_phase_'+str(step))['absolute_yield_magnitude']
        neg=next(x for x in rows if x['tag']=='M1_phase_'+str(-step))['absolute_yield_magnitude']
        return (math.log(pos)-math.log(neg))/(2*step)
    response=prior['phase_response'].copy();response['dlnY_dphase']=slope(.001)
    response['half_step_dlnY_dphase']=slope(.0005)
    response['late_population_response_held_fixed']=True
    response['full_isocurvature_transfer_computed']=False
    response['required_phase_sqrt_power_if_linear_transfer']=math.sqrt(response['allowed_charge_isocurvature_power'])/abs(slope(.001))
    ratio=old.INP['light_sector']['vPhi_GeV']/central['physical_H_GeV']
    response['early_deformation_cubic_coefficient_epsilon_v_over_H']=.002*ratio
    response['early_background_force_epsilon_v3_over_H3']=.002*ratio**3
    np.savez_compressed(ROOT/'BACKGROUND_M1.npz',**arrays)
    save('DEFORMED_MATCHED_RESULTS.json',{'scope':'actual changed ground profile and homogeneous early-history matching; formation not computed',
        'target_YQ_for_all_reference_objects_no_pairs':required,'target_DM_per_visible_entropy_GeV':D,
        'object_energy_per_charge_GeV':epq,'rows':rows,'columns':prior['columns'],'phase_response':response,
        'early_deformation_terms_neglected':response['early_deformation_cubic_coefficient_epsilon_v_over_H'],
        'physical_final_abundance':None})


if __name__=='__main__':main()
