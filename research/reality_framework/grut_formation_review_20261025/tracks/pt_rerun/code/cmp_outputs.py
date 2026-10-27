"""Compare regenerated outputs of one IR12 case with published copies (reproduction audit)."""
import json,sys,numpy as np
pub,new,tag=sys.argv[1:4]
out={}
def flat(x,p=''):
    if isinstance(x,dict):
        for k,v in x.items(): yield from flat(v,p+'/'+k)
    elif isinstance(x,list) and x and isinstance(x[0],(dict,list)):
        for i,v in enumerate(x): yield from flat(v,p+f'[{i}]')
    else: yield p,x
for win in ['14_24','24_34']:
    A=json.load(open(f'{pub}/WINDOW_{tag}_{win}.json'));B=json.load(open(f'{new}/WINDOW_{tag}_{win}.json'))
    fa=dict(flat(A));fb=dict(flat(B));worst=0;wk=None;rows={}
    for k,a in fa.items():
        b=fb.get(k)
        if a is None or isinstance(a,str):continue
        a=np.array(a,float);b=np.array(b,float)
        if a.size==0:continue
        d=float(np.max(abs(a-b)/np.maximum(abs(a),1e-12)))
        if d>worst:worst=d;wk=k
    for key in ['coherent_negative_frequency_action','wave_negative_frequency_action','total_negative_frequency_action','spectral_charge_ratio','equal_time_full_charge_ratio','omitted_k_gt8_Hartree_negative_reference','low_wave_fraction_kcenters_le1p6','fundamental_shell_wave_fraction','frequency_resolution']:
        rows[key]=dict(published=A[key],rerun=B[key],abs_diff=abs(A[key]-B[key]))
    out[win]=dict(max_relative_field_diff=worst,where=wk,headline=rows)
    sa=np.load(f'{pub}/SPECTRUM_{tag}_{win}.npz');sb=np.load(f'{new}/SPECTRUM_{tag}_{win}.npz')
    P=sa['power_by_shell'];R=sb['power_by_shell']
    out[win]['spectrum_max_abs_diff_over_max']=float(np.max(abs(P-R))/np.max(P))
a=np.load(f'{pub}/AVERAGE_{tag}.npz');b=np.load(f'{new}/AVERAGE_{tag}.npz');av={}
for f in a.files:
    x=a[f];y=b[f];av[f]=float(np.max(abs(x-y))/max(np.max(abs(x)),1e-300))
out['AVERAGE_max_rel_diff_per_array']=max(av.values());out['AVERAGE_worst_array']=max(av,key=av.get)
print(json.dumps(out,indent=1))
