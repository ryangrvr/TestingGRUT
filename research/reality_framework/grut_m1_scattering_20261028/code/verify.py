"""Bounded numerical checks; a passing kernel is not a cosmological verdict."""
from common import *
from scattering import scatter
from vacuum_kernels import matrix_pair,pair_cross_sections
from numpy.polynomial.legendre import leggauss

def main():
    vertices=json.loads((ROOT/'VERTICES_RESULTS.json').read_text())
    linear=json.loads((ROOT/'SCATTERING_RESULTS.json').read_text())
    continuum=json.loads((ROOT/'CONTINUUM_CHECKS.json').read_text())
    partial=json.loads((ROOT/'PARTIAL_RESULTS.json').read_text())
    source=json.loads((ROOT/'SOURCE_RESULTS.json').read_text())
    repro=json.loads((ROOT/'SOURCE_REPRODUCTION.json').read_text())
    weights=json.loads((ROOT/'SOURCE_WEIGHT_CHECKS.json').read_text())
    packets=json.loads((ROOT/'PACKET_RESULTS.json').read_text())
    groups=[]
    def add(name,ok,**extra):groups.append({'name':name,'pass':bool(ok),**extra})
    add('vacuum_vertex_and_decay_normalization',abs(vertices['angular_phase_space_check_relative'])<1e-12 and abs(vertices['real_species_check_relative'])<1e-12
        and abs(vertices['old_kinetic_input']['Gamma_s_vac_GeV']/vertices['width_GeV']-1)<1e-12,
        width_over_mX=vertices['width_over_mX'],old_kinetic_already_agrees=True)
    cc=linear['controls'];unit=max(r['unitarity_error'] for r in linear['rows'])
    forward_recheck=max(abs(scatter(r['k'])[key]-r[key]) for r in linear['rows'] for key in ('P_X','P_antiparticle','P_phi'))
    add('coupled_flux_and_null_controls',unit<1e-9 and max(r['P_X']+r['P_phi'] for r in cc)<1e-18,
        monopole_unitarity=unit,controls=cc,forward_response_unchanged_by_inverse_generalization=forward_recheck)
    inv=json.loads((ROOT/'INVERSE_RESULTS.json').read_text())
    reciprocity=max(r['matrix_reciprocity_error'] for r in inv['full_matrix_checks'])
    invunit=max(r['matrix_unitarity_error'] for r in inv['full_matrix_checks'])
    balance=max(r['k_squared_detailed_balance_relative_error'] for r in inv['rows'])
    add('full_open_channel_unitarity_and_reverse_detailed_balance',reciprocity<1e-10 and invunit<1e-10 and balance<1e-9 and forward_recheck<1e-12,
        reciprocity_error=reciprocity,full_matrix_unitarity_error=invunit,inverse_cross_section_balance_error=balance)
    difference=max(max(abs(r['P_X_difference']),abs(r['P_phi_difference'])) for r in continuum['rows'])
    add('separate_continuum_collocation',difference<.002 and max(r['continuum']['unitarity_error'] for r in continuum['rows'])<1e-8,
        largest_absolute_probability_difference_FD025= difference,
        scope='same author, separate continuum method; FD2 truncation remains reported')
    grid=max(abs(r['h025_sigma']/r['base_sigma']-1) for r in partial['refinements'])
    tail=max(abs(r['l48_sigma']/r['base_sigma']-1) for r in partial['refinements'])
    remote=max(abs(r['far1000_sigma']/r['base_sigma']-1) for r in partial['refinements'])
    add('angular_sum_and_source_quadrature',grid<.005 and tail<1e-7 and remote<1e-4
        and max(abs(r['relative_quadrature_change']) for r in weights['rows'])<1e-4,
        spacing_change=grid,l32_to_l48_change=tail,exterior_matching_change=remote)
    resonance=json.loads((ROOT/'RESONANCE_CHECKS.json').read_text());resrows=resonance['rows']
    padding_change=max(abs(r['controls'][-1]['resolved_average']/r['controls'][-2]['resolved_average']-1) for r in resrows)
    source_change=max(abs(r['controls'][-1]['resolved_resonance_finer_h_average']/r['controls'][-1]['resolved_average']-1) for r in resrows)
    add('resolved_low_k_resonance_and_rate_padding',padding_change<.001 and source_change<.001,
        padding16_to32_change=padding_change,resonance_mesh_source_average_change=source_change,
        all_possible_resonances_certified=False)
    padding=[r['Nminus_out'] for r in source['padding_controls']]
    add('source_reproduction_and_flux_transform',max(r['archived_history_max_scaled_difference'] for r in repro['rows'])<1e-12
        and max(r['net_out_vs_recorded_probe_relative'] for r in source['rows'])<1e-4 and np.ptp(padding)/np.mean(padding)<.002,
        pulse_history_reproduction='exact at printed precision',padding_relative_range=float(np.ptp(padding)/np.mean(padding)),
        largest_source_FFT_time_identity=max(r['FFT_time_charge_identity_absolute'] for r in source['rows']))
    edrift=max(r['energy_receiver_error'] for r in packets['runs'].values())
    qdrift=max(r['charge_receiver_absolute_error_in_quanta'] for r in packets['runs'].values())
    add('nonlinear_packet_receiver',edrift<1e-4 and qdrift<1e-6,max_energy_error=edrift,max_absolute_charge_error=qdrift)
    rows=packets['rows'];response=max(max(abs(x) for x in r['absolute_response_differences']) for r in rows)
    number=max(abs(r['measured_number_flux_sum']-1) for r in rows)
    shifts=max(abs(r['late_core_charge_shift']-r['charge_shift_from_outgoing_channels'])/r['incoming_number_measured'] for r in rows if r['late_core_charge_shift'] is not None)
    free=next(r for r in rows if r['tag']=='free10')
    add('finite_packets_match_linear_channels',response<.003 and number<.001 and shifts<.0002 and free['measured_P_X_anti_phi'][0]<1e-6,
        max_absolute_probability_difference=response,max_number_flux_error=number,max_core_shift_error_per_incident=shifts,
        caution='coarse/fine packet PX differs by about4% relative; no sub-percent full packet continuum certification')
    nodes,w=leggauss(192);errors=[]
    kernel=json.loads((ROOT/'VACUUM_KERNELS.json').read_text())
    for r in kernel['rows']:
        s=r['sqrt_s_over_mX']**2
        angular=np.dot(w,np.array([matrix_pair(s,z,.002)**2 for z in nodes]))
        p=np.sqrt(s/4-1);f=np.sqrt(s/4-6)
        sigma=angular/(64*np.pi*s)*f/p
        errors.append(abs(sigma/r['sigma_XbarX_to_phiphi_mXminus2']-1))
    add('limited_vacuum_kernel_quadrature_and_balance',max(errors)<1e-10 and max(abs(r['detailed_balance_relative_error']) for r in kernel['rows'])<1e-12,
        Gauss_vs_adaptive_quadrature=max(errors),full_medium_operator_complete=False)
    cold=json.loads((ROOT/'COLD_COMPLETION_RESULTS.json').read_text())['rows']
    add('conditional_cold_charge_and_finite_mass_ledger',all(abs(r['conditional_final_core_charge']+r['free_positive_X']-QREF)<1e-9 and r['bound'] and r['max_BVP_residual']<3e-9 for r in cold),
        conditional_mass_ratios=[r['cold_mass_over_reference'] for r in cold],
        processing_assumptions_have_not_been_dynamically_established=True)
    out={'groups':groups,'passed':sum(r['pass'] for r in groups),'total':len(groups),
         'independent_author_review':False,'actual_final_abundance':None,'actual_capture_history':None,
         'complete_quantum_capture_operator':False,'fixed_M1_couplings_retuned':False}
    save('VERIFICATION.json',out);print(json.dumps({'passed':out['passed'],'total':out['total'],'abundance':None},indent=2))
    if out['passed']!=out['total']:raise RuntimeError('Verification group failed')

if __name__=='__main__':main()
