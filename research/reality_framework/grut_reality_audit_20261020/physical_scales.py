"""Map inherited illustrative parameters to SI units; this is not calibration."""
from pathlib import Path
import hashlib
import json
import math

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'grut_candidate1_two_object_20261019/source_observables.json'


def main():
    source=json.loads(SOURCE.read_text())
    # SI defining constants and the measured 2022 CODATA G, checked at NIST.
    c=299792458.
    h=6.62607015e-34
    eV=1.602176634e-19
    G=6.67430e-11
    hbar=h/(2*math.pi)
    GeV=1e9*eV
    reduced_Planck_kg=math.sqrt(hbar*c/(8*math.pi*G))
    reduced_Planck_GeV=reduced_Planck_kg*c*c/GeV
    hbar_GeV_s=hbar/GeV
    hbarc_GeV_m=hbar*c/GeV
    g=.1

    def map_scale(v_GeV,label):
        mX=g*v_GeV
        mass_GeV=v_GeV/g*source['Ebar']
        mass_kg=mass_GeV*GeV/(c*c)
        radius=source['R90_bar']*hbarc_GeV_m/mX
        time_unit=hbar_GeV_s/mX
        return dict(label=label,v_GeV=v_GeV,v_over_reduced_Planck=v_GeV/reduced_Planck_GeV,
                    mX_GeV=mX,mphi_GeV=math.sqrt(.02)*v_GeV,mD_GeV=math.sqrt(.01+.0001)*v_GeV,
                    soliton_mass_GeV=mass_GeV,soliton_mass_kg=mass_kg,R90_m=radius,
                    dimensionless_time_unit_seconds=time_unit,tau240_seconds=240*time_unit,
                    tau160_to240_interval_seconds=80*time_unit,
                    compactness_Rs_over_R90=2*G*mass_kg/(c*c*radius),
                    geometric_identical_core_area_over_mass_cm2_per_g=4*math.pi*radius*radius/mass_kg*10,
                    geometric_warning='Geometric contact-area proxy only. This is not a calculated velocity/phase-dependent transfer cross section or an exclusion test.')

    output=dict(status='Conditional unit conversion of model inputs. No masses, sizes or retention times have been observed or fitted.',
                constants=dict(c_m_per_s=c,h_J_s=h,eV_J=eV,G_SI=G,
                               G_relative_standard_uncertainty=.00015/6.67430,
                               source='https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
                               reduced_Planck_mass_GeV=reduced_Planck_GeV),
                scale_independent_relations=dict(mphi_over_mX=math.sqrt(2),mD_over_mX=math.sqrt(1.01),
                                                soliton_mass_over_mX=source['Ebar']/(g*g),
                                                R90_times_mX_over_hbar_c=source['R90_bar'],
                                                tau240_in_reference_R90_light_crossings=240/source['R90_bar'],
                                                tau160_to240_in_reference_R90_light_crossings=80/source['R90_bar'],
                                                reference_R90_m_if_tau240_equals_one_second=c*source['R90_bar']/240,
                                                E_over_mXQ=source['Ebar_over_Qbar']),
                original_benchmark=map_scale(1e-5*reduced_Planck_GeV,'Inherited illustrative v/M_P=1e-5'),
                alternative_scale_example=map_scale(10.,'Illustration only: choose mX=1 GeV; no data select this value'),
                source_observables_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                limitations=['The measured G uncertainty enters the original physical scale; printed digits are computational output, not experimental precision.',
                             'The flat dimensionless scalar tests do not determine v or the new couplings independently.',
                             'No scalar-Higgs portal is present in the displayed leading benchmark potential. All ordinary-sector communication there is through the metric.',
                             'Allowed portals, higher operators and gravity loops in the EFT require declared matching and their own quantitative analysis.',
                             'A finite simulated duration is not a measured or predicted lifetime.'])
    (ROOT/'PHYSICAL_SCALES.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
