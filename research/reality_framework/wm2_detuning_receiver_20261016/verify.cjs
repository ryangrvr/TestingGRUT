"use strict";
function verifyEvidence(S,r,a){
  const assert=S.assert;
  function finite(x){if(typeof x==="number")assert(Number.isFinite(x));else if(Array.isArray(x))x.forEach(finite);else if(x&&typeof x==="object")Object.values(x).forEach(finite);}
  finite(r);assert(r.status==="PASS");assert(JSON.stringify(r.constants)===JSON.stringify(S.constants()));
  assert(r.graph.length===3&&r.reduction.length===12&&r.continuum.length===90&&r.finite_receiver.length===27);
  assert(r.Schur.length===3&&r.source.length===5&&r.propagation.length===2&&r.Q_direct.length===2);
  for(const x of r.graph){
    assert(x.full_one_probe_dimension===x.reachable_basis_states&&x.full_one_probe_dimension===x.independent_register_states);
    assert(x.root_moving_covariance_checks===x.full_one_probe_dimension);
    assert(x.PQ_norm_squared_upper_checked<=.0625+1e-14&&x.emission_append_norm_squared_upper_checked<=.25+1e-14);
    assert(x.minimum_virtual_word_length_minus_distance===2&&x.metric_matches_path_including_all_protected_geometry_moves);
  }
  for(const x of r.reduction)for(const z of x.rows){
    assert(z.full_wavefunction_error<2e-10&&Math.abs(z.full_norm-1)<2e-10);
    assert(Math.abs(z.macroscopic_probability-z.binomial_probability)<2e-10);
  }
  for(const x of r.continuum){
    assert(x.delta>=-.5&&x.delta<=.5&&x.g>=.9&&x.g<=1.1);
    for(const z of x.rows){
      assert(Math.abs(z.density_mass-1)<2e-8&&z.pump_survival<.25);
      assert(z.macroscopic_probability>=.5&&z.complete_record_entropy_bits>=0);
      assert(z.t>=4*x.N*x.N&&z.t<=8*x.N*x.N);
    }
  }
  for(const x of r.finite_receiver){
    assert(x.receiver_modes===64*x.N*x.N&&x.one_particle_dimension===x.receiver_modes+1&&x.initial_receiver==="vacuum");
    for(const z of x.rows){
      assert(z.macroscopic_probability>=.5&&z.pump_survival<.25&&z.finite_vs_continuum_amplitude_error<2e-8);
      assert(z.unitarity_error<1e-9&&Math.abs(z.complete_energy_ledger_error)<2e-7&&z.complete_geometry_record_entropy_bits>=0);
      assert(Math.abs(z.mean_geometry_record_work)<=1.1+1e-8);
    }
  }
  for(const x of r.quadrature_refinement)assert(x.amplitude_change<1e-9);
  for(const x of r.Schur)for(const z of x.rows){
    assert(z.eta>0&&z.Sigma_Frobenius_upper<=z.operator_bound+1e-10&&z.solve_residual<1e-9);
    assert(Math.abs(z.operator_bound-1/(8*z.eta)*(1+5*z.eta/21)**(-z.min_distance))<1e-14);
  }
  for(const x of r.source){
    assert(x.variance===.25&&x.field_state.includes("no reservoir reset"));
    for(const z of x.rows){
      assert(z.pulse_identity_error<1e-9&&z.bulk_prefix_failure_probability<=z.bulk_prefix_bound+1e-10);
      assert(z.mismatch_probability<=z.mismatch_Duhamel_upper+1e-10);
      if(z.near_response_lower>0)assert(z.chi_all_storage_and_receiver_modes[1]>=z.near_response_lower-1e-10);
      for(const [k,v] of Object.entries(z.receiver_response_bounds))assert(Math.abs(z.chi_all_storage_and_receiver_modes[Number(k)])<=v+1e-10);
    }
  }
  for(const x of r.propagation)for(const z of x.rows)assert(z.all_receiver_Frobenius_upper<=z.operator_Dyson_bound+1e-10);
  for(const x of r.Q_near_resonance){
    assert(x.N>=64&&x.Q_residual_norm<=x.analytic_Q_residual_upper&&x.coupled_real_pole_window_radius<2.5);
    assert(x.contact_occupation_exactly_zero&&!x.real_gap_uniform_in_N);
  }
  for(const x of r.Q_direct)assert(x.direct_vs_formula_error<2e-8&&x.direct_contact_amplitude===0);
  assert(a.explicit_closed_energy_receiver&&a.macroscopic_formation.startsWith("PASS IN CLASS"));
  assert(a.real_energy_complement_gap.startsWith("FAIL")&&a.source_response_front&&a.all_receiver_storage_geometry_states_in_bounds);
  for(const k of ["isolated_low_energy_consistency_band","irreversible_lock","selected_dimension","original_CP12_matter_universality",
    "generic_local_perturbation_robustness","full_WM2_success","all_WM2_architectures_ruled_out","frozen_CP12_modified","S2_modified",
    "S1_modified","WM1_modified","SCOUT1_reopened","empirical_prediction","earned_quotient_value","large_CP12_trajectory","external_independent_validation"])assert(a[k]===false,k);
  const flat=(xs)=>xs.flatMap(x=>x.rows);
  return {numerical_consistency:"PASS",physical_acceptance:"Restricted operational stability interval PASS; isolated real-energy consistency-band gate FAIL; full WM2 NOT ESTABLISHED",
    complete_register_configuration_states:r.graph.reduce((s,x)=>s+x.full_one_probe_dimension,0),
    reverse_primitive_checks:r.graph.reduce((s,x)=>s+x.reverse_primitive_checks,0),
    virtual_exit_distance_checks:r.graph.reduce((s,x)=>s+x.virtual_exit_distance_checks,0),
    largest_full_one_probe_matrix:Math.max(...r.graph.map(x=>x.full_one_probe_dimension)),
    full_wavefunction_reduction_checks:flat(r.reduction).length,
    largest_full_vacuum_matrix:Math.max(...r.reduction.map(x=>x.full_vacuum_configuration_dimension)),
    continuum_parameter_sets:r.continuum.length,continuum_observations:flat(r.continuum).length,
    finite_resource_scaled_parameter_sets:r.finite_receiver.length,finite_resource_scaled_observations:flat(r.finite_receiver).length,
    largest_finite_resource_scaled_one_particle_matrix:Math.max(...r.finite_receiver.map(x=>x.one_particle_dimension)),
    complex_Schur_blocks:flat(r.Schur).length,finite_phase_pulse_checks:flat(r.source).length,
    full_H_propagation_blocks:flat(r.propagation).length,Q_packet_formula_checks:r.Q_near_resonance.length,Q_packet_direct_checks:r.Q_direct.length,
    largest_direct_receiver_modes:Math.max(...r.Q_direct.map(x=>x.receiver_modes)),
    max_full_wavefunction_error:Math.max(...flat(r.reduction).map(x=>x.full_wavefunction_error)),
    max_finite_continuum_amplitude_difference:Math.max(...flat(r.finite_receiver).map(x=>x.finite_vs_continuum_amplitude_error)),
    max_quadrature_refinement_difference:Math.max(...r.quadrature_refinement.map(x=>x.amplitude_change)),
    max_complex_solve_residual:Math.max(...flat(r.Schur).map(x=>x.solve_residual)),
    max_pulse_identity_error:Math.max(...flat(r.source).map(x=>x.pulse_identity_error)),
    max_complete_energy_ledger_error:Math.max(...flat(r.finite_receiver).map(x=>Math.abs(x.complete_energy_ledger_error))),
    execution:r.execution,arithmetic:r.arithmetic,external_independent_validation:false};
}
if(typeof module!=="undefined"&&module.exports)module.exports=verifyEvidence;
