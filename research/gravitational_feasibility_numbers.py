"""Arithmetic for one apparatus assessment, not experimental data or a law test.

Run from any directory: python research/gravitational_feasibility_numbers.py
Outputs a deterministic JSON record. Source configurations are kept separate.
"""
import json
import math


def calculate():
    h = 6.62607015e-34  # exact SI, J s
    k_b = 1.380649e-23  # exact SI, J/K
    hbar = h / (2 * math.pi)
    modes = []
    for label, frequency, quality, temperature, reported_floor in [
        ("mode_3", 50.59, 3.8e6, 7.1e-3, 2.7e5),
        ("mode_4", 67.98, 5.5e6, 6.6e-3, 5.7e5),
    ]:
        omega = 2 * math.pi * frequency
        modes.append({
            "mode": label,
            "source": "https://arxiv.org/pdf/2605.28479v2",
            "source_locations": "pp. 3-5, Figs. 3-4, Eqs. (4)-(5)",
            "input_frequency_Hz": frequency,
            "input_bare_Q": quality,
            "input_feedback_mode_temperature_K": temperature,
            "input_reported_feedback_floor_phonons_approximate": reported_floor,
            "quantum_energy_temperature_scale_K": h * frequency / k_b,
            "energy_equivalent_occupation": k_b * temperature / (h * frequency),
            "bare_energy_damping_rate_s_inv": omega / quality,
            "bare_energy_ringdown_s_NOT_quantum_coherence_time": quality / omega,
            "zero_point_motion_m_envelope": [
                math.sqrt(hbar / (2 * mass * omega))
                for mass in (3.8e-7, 3.5e-7)
            ],
            "mass_envelope_status": "Text reports 0.35 mg on p.2 and 0.38 mg on p.5; not a statistical confidence interval",
            "noise_ASD_product_improvement_for_floor_1_with_same_linear_feedback_model": reported_floor,
            "noise_PSD_product_improvement_for_floor_1_with_same_model": reported_floor ** 2,
        })
    noise_asd = 0.5e-15
    signal = 30e-18
    return {
        "assessment_date_utc": "2026-10-11",
        "verdict": "NOT YET FEASIBLE",
        "scope": "Documented Leiden apparatus; gravity-mediated coherence/entanglement, not ordinary force detection",
        "constants": {"h_J_s": h, "hbar_J_s": hbar, "k_B_J_per_K": k_b},
        "formulas": {
            "energy_equivalent_occupation": "k_B*T_mode/(h*f); does not certify a thermal quantum state",
            "zero_point_motion": "sqrt(hbar/(2*m*2*pi*f))",
            "bare_energy_ringdown": "Q/(2*pi*f); feedback damping excluded",
            "linear_feedback_floor": "n_min=sqrt(S_F*S_x_det)/(2*hbar), one-sided PSDs, published model",
            "classical_force_integration": "T=(z*noise_ASD/signal)^2 for known-phase sinusoid, one-sided locally white stationary noise",
        },
        "cooling_configuration_2026": modes,
        "gravity_configuration_2024": {
            "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10889343/",
            "source_locations": "Results, Discussion, Fig.2; arXiv:2303.03545v2 includes supplements",
            "input_test_mass_kg": 0.43e-6,
            "input_resonance_Hz": 26.7,
            "input_force_noise_ASD_N_per_sqrt_Hz": noise_asd,
            "input_reported_gravitational_signal_N": signal,
            "conditional_white_noise_force_5sigma_integration_s": (5 * noise_asd / signal) ** 2,
            "integration_is_quantum_forecast": False,
            "integration_has_stationarity_or_systematic_certificate": False,
        },
        "missing_inputs": ["prepared source coherence and branch geometry", "protocol-specific coherence lifetime", "second quantum-controlled mass", "joint quantum readout with errors", "background coupling limits in the intended geometry"],
        "quantum_signature_shot_count": None,
        "reason_shot_count_undefined": "No apparatus-specific signal contrast, prepared state or joint readout distribution established",
        "independent_scientific_review": "PENDING",
        "new_GRUT_law": False,
        "card3_authorized": False,
    }


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2, allow_nan=False))
