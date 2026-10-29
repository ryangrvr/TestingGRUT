"""Create a compact numerical ledger without treating missing energy as lost."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from expansion import clock, EREF, QREF

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('tags',nargs='+');args=p.parse_args()
    f=json.loads((ROOT/'FAMILY.json').read_text())['nodes']
    mass=CubicHermiteSpline([a['Q'] for a in f],[a['energy_mX'] for a in f],[a['Omega'] for a in f])
    rows=[]
    for tag in args.tags:
        d=json.loads((ROOT/(tag+'.json')).read_text());a=d['snapshots'][0];b=d['snapshots'][-1]
        source=d['configuration']['source']
        if source=='parent':
            pulse=json.loads((ROOT.parent/'grut_m1_kinetics_20261028'/'PHYSICS_CHECKS.json').read_text())['source_energy_ledger']['pulse_initial_lattice_energy_mX']
        else:pulse=json.loads((ROOT/'sources'/f'ringdown_{source.removesuffix("early")}.json').read_text())['initial_pulse_energy_mX']
        residual=pulse-float(mass(a['Qcore']))-a['gas_energy_mX']
        cold=b['stationary_cold_excess_fraction']
        physical=pulse-b['work_redshift_mX']+b['interpolated_generator_energy_defect_mX']
        components=100*(cold+(b['gas_kinetic_mX']+residual)/EREF)
        total=100*(physical/EREF-1)
        rows.append({'tag':tag,'source':source,'a_over_a0':b['a_over_a0'],'Qcore':b['Qcore'],
                     'core_charge_percent':100*b['Qcore']/QREF,'Nplus':b['Nplus'],'Nminus':b['Nminus'],
                     'cold_excess_percent':100*cold,'kinetic_mX':b['gas_kinetic_mX'],
                     'unassigned_initial_energy_mX':residual,'unassigned_energy_percent':100*residual/EREF,
                     'redshift_work_mX':b['work_redshift_mX'],
                     'generator_defect_mX':b['interpolated_generator_energy_defect_mX'],
                     'total_excess_including_held_unassigned_percent':total,
                     'component_sum_mismatch_percentage_points':components-total,
                     'floor_particles':b['floor_particles'],'omitted_events':b['omitted_events']})
    c=clock();m=c['mX_GeV'];N=c['ncore0_GeV3']/m**3;H=c['H0_GeV'];spacing=N**(-1/3)
    s=json.loads((ROOT/'sources'/'ringdown_h0025.json').read_text())['signed_spectrum_diagnostics']
    nh=sum(s[w]['neutral']['signed_number'] for w in ('0_1000','1000_10000'))
    eh=sum(s[w]['neutral']['frequency_times_number_energy_mX'] for w in ('0_1000','1000_10000'))/nh
    width=.00018756589919939715;g=eh/np.sqrt(6)
    # Union bound on causal balls when an initially localized pulse is emitted
    # simultaneously from one core per reference volume. Neglect the small
    # initial core radius relative to 10000/mX; overlaps only reduce coverage.
    cover=N*4*np.pi/3*10000**3
    clockinfo={'ncore_over_mX3':N,'H_over_mX':H/m,'reference_separation_mX_inverse':spacing,
      'light_crossing_reference_separation_optical_time':N*spacing,'ringdown_10000_optical_time':N*10000,
      'ringdown_10000_Hubble_fraction':H/m*10000,'representative_neutral_energy_mX':eh,
      'representative_neutral_mean_decay_distance_mX_inverse':np.sqrt(g*g-1)/width,
      'representative_neutral_mean_decay_distance_over_separation':np.sqrt(g*g-1)/width/spacing,
      'synchronous_pulse_t10000_causal_volume_fraction_upper_bound':cover,
      'corresponding_minimum_average_density_enhancement_inside_covered_volume':1/cover,
      'meaning':'Mean energy inserted into a single-particle decay-length expression; actual decay distances and momenta have a distribution.'}
    out={'reference_energy_mX':EREF,'rows':rows,'clock_and_geometry':clockinfo,
         'scope':'Supplied pulse and density, quasistatic cores, isotropic homogeneous transport; source remnant held fixed as a bookkeeping component. Neither the residual nor a 1% mass gate is physically resolved.'}
    (ROOT/'SUMMARY.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
