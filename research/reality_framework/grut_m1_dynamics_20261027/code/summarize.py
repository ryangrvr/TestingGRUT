"""Preserve conditional outcomes; do not manufacture a cosmological abundance."""
from common import *
from scipy.optimize import linear_sum_assignment
from analysis_io import read_analysis

def match(a,b,L):
    if not a or not b:return []
    x=np.array([v['centre'] for v in a]);y=np.array([v['centre'] for v in b])
    d=(x[:,None,:]-y[None,:,:]+L/2)%L-L/2;distance=np.linalg.norm(d,axis=-1)
    ia,ib=linear_sum_assignment(distance)
    return [(int(i),int(j),float(distance[i,j])) for i,j in zip(ia,ib) if distance[i,j]<8.]

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases']
    data=read_analysis()['cases'];cases={}
    for tag,d in data.items():
        rows=d['rows'];last=rows[-1];late=[r for r in rows if r['t']>=last['t']-60]
        cases[tag]={'initial_E_over_Q':runs[tag]['parameters']['initial_E_over_Q'],
                    'final_E_over_Q':last['E']/last['Q'],'energy_drift':runs[tag]['energy_relative_drift'],
                    'charge_drift':runs[tag]['charge_relative_drift'],
                    'max_axial_UV75_fraction':max(r['axial_UV75_fraction_kinetic_gradient'] for r in rows),
                    'max_negative_phi_fraction':max(r['phi_negative_volume_fraction'] for r in rows),
                    'last_contrast':last['charge_contrast'],
                    'catalogs':[{k:c[k] for k in ('threshold','count','mask_charge_fraction')} for c in last['catalogs']],
                    'mask_fraction_is_final_eta':False}
        if 'central_spheres' in last:
            cases[tag]['central_spheres']=last['central_spheres']
            cases[tag]['late_radius12_charge_envelope']=[min(r['central_spheres'][1]['Q'] for r in late),max(r['central_spheres'][1]['Q'] for r in late)]
            cases[tag]['finite_object_energy_vs_reference']=last['E']/last['Q']/.9373999905262296-1
        else:
            tracks=[]
            for left,right in zip(late[:-1],late[1:]):
                aa=left['catalogs'][1]['objects'];bb=right['catalogs'][1]['objects']
                links=match(aa,bb,last['L'])
                tracks.append({'t_from':left['t'],'t_to':right['t'],'matched_candidates':len(links),
                               'distance_max':max([v[2] for v in links],default=0.),
                               'largest_charge_relative_changes':[float(bb[j]['Q_mask']/aa[i]['Q_mask']-1) for i,j,_ in links[:5]]})
            cases[tag]['late_candidate_links']=tracks
            cases[tag]['last_compact_candidates']=sum(o['compact'] for o in last['catalogs'][1]['objects'])
            cases[tag]['largest_mask_charges']=[o['Q_mask'] for o in last['catalogs'][1]['objects'][:8]]
    comparison={}
    for a,b in [('cold64','cold48'),('cold64','cold64_half'),('dilute64','dilute48'),('object64_perturbed','object96_perturbed'),('ellipse64','ellipse96'),('ellipse64','ellipse64_half')]:
        if a not in data or b not in data:continue
        x=data[a]['rows'][-1];y=data[b]['rows'][-1]
        comparison[a+'_vs_'+b]={'final_energy_per_charge_relative_difference':(y['E']/y['Q'])/(x['E']/x['Q'])-1,
            'charge_contrast_relative_difference':y['charge_contrast']/x['charge_contrast']-1,
            'mask4_charge_fraction_difference':y['catalogs'][1]['mask_charge_fraction']-x['catalogs'][1]['mask_charge_fraction']}
        if a==b.removesuffix('_half') and b.endswith('_half'):
            pa=runs[a]['snapshots'][-1];pb=runs[b]['snapshots'][-1]
            qa=np.asarray(snapshot(pa)[3]);qb=np.asarray(snapshot(pb)[3])
            comparison[a+'_vs_'+b]['final_field_RMS_relative_difference']=float(np.linalg.norm(qa-qb)/np.linalg.norm(qb))
    exterior=json.loads((ROOT/'EXTERIOR_RESULTS.json').read_text())['cases']
    ext={}
    for tag,c in exterior.items():
        rows=c['rows'];late=[r for r in rows if r['t']>=160 and r['cut_radius']>=16]
        ext[tag]={'late_Nminus_over_Q_envelope':[min(r['Nminus_over_original_Q'] for r in late),max(r['Nminus_over_original_Q'] for r in late)],
                  'late_interaction_fraction_max':max(r['absolute_interaction_over_free_energy'] for r in late),
                  'late_absolute_taper_energy_fraction_max':max(abs(r['taper_gradient_energy_fraction']) for r in late),
                  'scope':'classical, weak-field exterior on this periodic object disturbance only; not cosmic p'}
    result={'cases':cases,'numerical_comparisons':comparison,'exterior':ext,
            'physical_final_eta':None,'physical_final_p':None,'physical_f_sym':None,
            'actual_object_charge_distribution':None,
            'frozen_reference_mass_criterion':'.0667805*(1-eta)+2.13356*p<=.01',
            'global_charge_ledger_extension':'For positive objects plus free gas, eta=Qobj/Qnet may exceed 1 when net charge in the exterior is negative. Require Nplus/Qnet=1-eta+p>=0; the same mass identity applies with actual object E(Q)/Q. This is not a new tolerance or a computed cosmological composition.'}
    save('RESULTS.json',result)
    print(json.dumps({'histories':len(cases),'comparisons':comparison,'exterior':ext},indent=2))

if __name__=='__main__':main()
