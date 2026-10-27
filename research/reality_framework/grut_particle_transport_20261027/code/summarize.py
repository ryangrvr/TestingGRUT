"""Combine the completed histories without promoting field actions to particles."""
from pt_common import *
from pt_production import CASES

def main():
    rows=[];traces={};maxE=maxQ=maxEX=maxQX=0.
    for case in CASES:
        tag=case['tag'];r=read(f'RUN_{tag}.json');w0=read(f'WINDOW_{tag}_14_24.json');w1=read(f'WINDOW_{tag}_24_34.json')
        maxE=max(maxE,r['energy_max_relative_drift']);maxQ=max(maxQ,r['charge_max_relative_drift'])
        maxEX=max(maxEX,max(x['energy_relative'] for x in r['export_checks']));maxQX=max(maxQX,max(x['charge_relative'] for x in r['export_checks']))
        t=np.load(WORK/tag/'trace.npz');traces[tag+'_s']=t['s'];traces[tag+'_tau']=t['tau'];traces[tag+'_values']=t['values'];traces['names']=t['names']
        row=dict(case=case,initial_q_RMS_over_mean=r['initial'].get('q_RMS_over_mean'),windows=[])
        for w in [w0,w1]:
            row['windows'].append(dict(s=w['window_s'],tau=w['window_tau'],coherent_negative_action=w['coherent_negative_frequency_action'],
                wave_negative_action=w['wave_negative_frequency_action'],total_negative_action=w['total_negative_frequency_action'],
                spectral_charge_ratio=w['spectral_charge_ratio'],reference_mass_means={m:w['equal_time'][f'ref_{m}']['hann2_mean'] for m in ['0','7','12','20','26','H']},
                omitted_Hartree_reference_above8=w['omitted_k_gt8_Hartree_negative_reference'],
                fundamental_shell_wave_fraction=w['fundamental_shell_wave_fraction'],low_wave_fraction_kcenters_le1p6=w['low_wave_fraction_kcenters_le1p6']))
        rows.append(row)
    np.savez_compressed(ROOT/'FIELD_TRACES.npz',**traces)
    mainrow=next(r for r in rows if r['case']['tag']=='base48');half=next(r for r in rows if r['case']['tag']=='base48_half')
    time_error=max(abs(a['total_negative_action']-b['total_negative_action']) for a,b in zip(mainrow['windows'],half['windows']))
    coarse=next(r for r in rows if r['case']['tag']=='IR12_s1_w1');fine=next(r for r in rows if r['case']['tag']=='IR12_fine_s1_w1')
    refinement=[dict(s=a['s'],relative_total_action_difference=abs(a['total_negative_action']/b['total_negative_action']-1),relative_wave_action_difference=abs(a['wave_negative_action']/b['wave_negative_action']-1)) for a,b in zip(coarse['windows'],fine['windows'])]
    responses=[read(p.name) for p in sorted(ROOT.glob('RESPONSE_s*.json'))]
    response_summary=dict(sampling=[dict(s0=r['s0'],seed=r['probe_seed']) for r in responses],
        maximum_canonical_sum_deviation=max(abs(x['Z_sum']-1) for r in responses for x in r['rows']),
        count_dominant_IR_anomalous_normal_ratio_range=[min(x['anomalous_over_normal_RMS'] for r in responses for x in r['rows'][:3]),max(x['anomalous_over_normal_RMS'] for r in responses for x in r['rows'][:3])],
        interpretation='Normal and conjugate responses mix strongly in the infrared; finite-probe and phase variation shown. This calls for matrix/coherent matching, not an inference that quasiparticles are impossible.')
    result=dict(rows=rows,conservation=dict(maximum_energy_relative_drift=maxE,maximum_charge_relative_drift=maxQ,
        averaged_export_energy_relative_error=maxEX,averaged_export_charge_relative_error=maxQX),
        time_refinement_maximum_total_action_absolute=time_error,IR_spatial_refinement=refinement,response=response_summary,
        negative_action_ranges_by_window=[dict(s=mainrow['windows'][j]['s'],minimum=min(r['windows'][j]['total_negative_action'] for r in rows),maximum=max(r['windows'][j]['total_negative_action'] for r in rows)) for j in range(2)],
        physical_nbar=None,physical_f_sym=None,physical_kinetic_threshold=None,full_gauged_cosmology_evolved=False,
        status='Computed regulated classical action budget and field-source export. Cosmic preparation, quantum/particle/chemical matching and later binding remain unresolved.')
    save('RESULTS.json',result)
    table=['# Classical action budget','', 'Each entry uses its own measured net current. These are Hann-windowed classical field actions, not a late physical antiparticle abundance.','',
        '| Case | Grid, box | s=14–24: coherent + waves = total | s=24–34: coherent + waves = total |','|---|---|---|---|']
    for r in rows:
        c=r['case'];v=[]
        for w in r['windows']:v.append(f"{w['coherent_negative_action']:.4f} + {w['wave_negative_action']:.4f} = **{w['total_negative_action']:.4f}**")
        table.append(f"| {c['tag']} | {c['N']} cubed, L={c['L']:g} | {v[0]} | {v[1]} |")
    (ROOT/'RESULT_TABLE.md').write_text('\n'.join(table)+'\n')
    source=dict(status='conditional_field_input',coherent_and_connected_reservoirs_separate=True,
        action_statistic='Lab-frame |frequency| times classical temporal power; no universal on-shell assignment.',
        action_ranges=result['negative_action_ranges_by_window'],reference_mass_scan='See RESULTS.json rows/windows/reference_mass_means.',
        coherent_matching_control='COHERENT_MATCHING.json; supplied Hamiltonian deformations, not the cosmic transition.',
        covariance_files=[f"AVERAGE_{r['case']['tag']}.npz" for r in rows],response_files=[p.name for p in sorted(ROOT.glob('RESPONSE_s*.npz'))],
        phase_space_gas_distribution=None,coherent_antiparticle_population=None,physical_nbar=None,physical_f_sym=None,
        abundance_denominator='Final net-charge object population energy sum n(Q)E(Q) at fixed target yield; not computed here.',
        symmetric_allowance=.01,selected_cosmic_quantum_state=False,full_model_included=False,
        next_dependency='Coupled coherent/connected quantum chemical evolution and physical mass/binding history; no arbitrary decision to include or discard coherent counter-rotation.')
    save('KINETIC_FIELD_SOURCE.json',source)
    print(json.dumps({k:result[k] for k in ['conservation','time_refinement_maximum_total_action_absolute','IR_spatial_refinement','negative_action_ranges_by_window','response']},indent=2),flush=True)

if __name__=='__main__':main()
