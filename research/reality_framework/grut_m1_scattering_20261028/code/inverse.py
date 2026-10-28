"""Direct reverse-channel calculation and reciprocity/detailed-balance checks."""
from common import *
from scattering import scatter

def main():
    rows=[];full=[]
    for k in (.4,1.,1.5,2.):
        energy=np.sqrt(1+k*k);kx=np.sqrt((energy+2*OMEGA)**2-1)
        forward=0.;reverse=0.;reverse_phi=0.;maxrec=0.
        for ell in range(33):
            f=scatter(k,ell=ell);b=scatter(k,ell=ell,incoming_index=0)
            forward+=(2*ell+1)*f['P_X'];reverse+=(2*ell+1)*b['P_antiparticle']
            reverse_phi+=(2*ell+1)*b['P_phi']
            maxrec=max(maxrec,abs(f['P_X']-b['P_antiparticle']))
        sf=np.pi/k**2*forward;sr=np.pi/kx**2*reverse
        vf=k/energy;vx=kx/(energy+2*OMEGA)
        rows.append({'reference_antiparticle_k':k,'reverse_incoming_positive_k':float(kx),
            'sigma_forward_mXminus2':float(sf),'sigma_reverse_mXminus2':float(sr),
            'sigma_positive_X_to_phi_mXminus2':float(np.pi/kx**2*reverse_phi),
            'sigma_positive_X_to_antiparticles_after_one_phi_decay_mXminus2':float(sr+np.pi/kx**2*reverse_phi),
            'partial_probability_reciprocity_error':maxrec,
            'k_squared_detailed_balance_relative_error':float(abs(k*k*sf/(kx*kx*sr)-1)),
            'per_particle_reverse_over_forward_rate_same_core_density':float(sr*vx/(sf*vf))})
        # Full open-channel matrix: retain phases and check column orthogonality.
        channels=[0,1]+([2] if energy+OMEGA>np.sqrt(6) else [])
        columns=[]
        for channel in channels:
            r=scatter(k,incoming_index=channel,return_amplitudes=True)
            amp=np.array([complex(*x) for x in r['amplitudes_real_imag']]);columns.append(amp[channels])
        matrix=np.array(columns).T
        full.append({'reference_k':k,'channels':channels,
                     'matrix_unitarity_error':float(np.max(abs(matrix.conj().T@matrix-np.eye(len(channels))))),
                     'matrix_reciprocity_error':float(np.max(abs(matrix-matrix.T)))})
    result={'scope':__doc__,'rows':rows,'full_matrix_checks':full,
        'inverse_threshold_incoming_positive_energy_over_mX':float(1+2*OMEGA),
        'inverse_threshold_incoming_positive_momentum_over_mX':float(np.sqrt((1+2*OMEGA)**2-1)),
        'warning':'fast converted X can regenerate antiX before cooling; forward Gamma/H alone is not erasure',
        'actual_forward_minus_inverse_kinetic_history':None}
    save('INVERSE_RESULTS.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
