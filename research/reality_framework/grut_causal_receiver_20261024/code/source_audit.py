#!/usr/bin/env python3
"""Audit emitted classical fields without filling unresolved modes with particles."""
import hashlib,importlib.util,json
from pathlib import Path
import numpy as np
from run import ROOT,gradient

OLD=ROOT.parent/'grut_soft_receiver_20261023/code/trace_audit.py'
spec=importlib.util.spec_from_file_location('prior_trace_method',OLD)
method=importlib.util.module_from_spec(spec);spec.loader.exec_module(method)


def main():
    rows=[]
    names=['P025','P025_time','P0125','P025_remote_time']
    if (ROOT/'P00625_RESULTS.json').exists():names.append('P00625')
    for name in names:
        run=json.loads((ROOT/f'{name}_RESULTS.json').read_text())
        trace=np.load(ROOT/'data'/f'{name}_trace.npz')['trace']
        assert hashlib.sha256(trace.tobytes()).hexdigest()==run['raw_numeric_sha256']['trace']
        variants=[method.analyze(trace,run['parameters'],alpha) for alpha in [.025,.05,.1]]
        variants.append(method.analyze(trace,run['parameters'],.05,6000))
        rows.append(dict(name=name,core60_late_mean_fraction=run['core60_late_mean_fraction'],
            maximum_layer_charged_mass2_shift=run['maximum_layer_charged_mass2_shift'],variants=variants))
    same=json.loads((ROOT/'P025_time_RESULTS.json').read_text())
    remote=json.loads((ROOT/'P025_remote_time_RESULTS.json').read_text())
    near_table=np.load(ROOT/'data/P025_time_final.npz')['final']
    far_table=np.load(ROOT/'data/P025_remote_time_final.npz')['final'][:len(near_table)]
    r=near_table[:,0];h=same['parameters']['h'];inside=r<300
    ne=r*(near_table[:,1]+1j*near_table[:,2]);fa=r*(far_table[:,1]+1j*far_table[:,2])
    nv=r*(near_table[:,4]+1j*near_table[:,5]);fv=r*(far_table[:,4]+1j*far_table[:,5])
    overlap=np.vdot(fa[inside],ne[inside]);phase=overlap/abs(overlap)
    raw=float(sum(abs(ne[inside]-fa[inside])**2)/sum(abs(ne[inside])**2))
    aligned=float(sum(abs(ne[inside]-phase*fa[inside])**2)/sum(abs(ne[inside])**2))
    out=dict(scope='New classical-field source audit. Conditional Fourier particle numbers do not certify unresolved soft modes or a cosmological source.',
        rows=rows,boundary_radius_control=dict(names=['P025_time','P025_remote_time'],same_h=h,same_dt=same['parameters']['dt'],same_end_time=same['parameters']['tmax'],
            core_charge_relative_shift=float(remote['core60_late_mean_fraction']/same['core60_late_mean_fraction']-1),
            raw_complex_field_L2_squared_relative_difference_inside_300=raw,
            explicitly_global_phase_aligned_L2_squared_relative_difference_inside_300=aligned,
            global_phase_alignment_radians=float(np.angle(phase))),
        caution='Vacuum reflection is controlled; nonlinear matching and finite-time/bound-core mode separation remain distinct questions. No missing energy is renormalized into particles.')
    (ROOT/'SOURCE_AUDIT_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for row in rows:
        v=row['variants'][1]
        print(row['name'],'conditional incoming/outgoing',v['conditional_incoming_over_outgoing'],
            'local-cut fraction',v['conditional_controlled_energy_fraction'],
            'conditional anti',v['conditional_all_band_antiparticle_number_fraction'],flush=True)
    print('Radius control',out['boundary_radius_control'],flush=True)


if __name__=='__main__':main()
