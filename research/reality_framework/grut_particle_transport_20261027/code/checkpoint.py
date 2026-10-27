"""Prepare a self-describing owner checkpoint; saving/versioning is separate."""
from pt_common import *

def main():
    tags=['base48','IR12_fine_s1_w1','IR24_s1_w1','IR24_s2_w1'];fields={}
    metadata=dict(s=34.,tau=float(clock(34.)),field_order=['chi_X_real','chi_X_imag','chi_phi','pi_X_real','pi_X_imag','pi_phi'],
        precision='float64',couplings=dict(lambda_X=LX,g_X_squared=G2,lambda_phi=LP,nu_over_H_i=0.,gauge_e=0.),
        units=dict(H_i_GeV=float(E.M.H0),conformal_time='s=H_i eta, with s=0 at tau=20',
            q='(sqrt2 Re(aX), sqrt2 Im(aX), a phi)/H_i',p='d q/ds; physical conformal momentum equals H_i^2 p',
            coordinate='Comoving coordinates in H_i^-1; L is the periodic box size.'),
        radiation_clock=dict(tau0=20.,a0=8.869218738806,C=1.9857962750634526,tau='tau0+a0*s+C*s^2/2',a='a0+C*s'),
        parent_commit='7d9c93bc1dffea59228765eed41f65baa4a9692c',review_commit='fd71d3a416db929252561d2cf296e915dd710f37',
        stepper_sha256=sha(ROOT/'code/fieldstep.cpp'),
        scope='Regulated classical scalar-limit states, not a quantum state, particle distribution, bound population or complete cosmological history.',cases={})
    for tag in tags:
        a=np.load(WORK/tag/'state_s34.npy');assert np.all(np.isfinite(a))
        fields[tag+'_fields']=a;metadata['cases'][tag]=read(f'RUN_{tag}.json')
    fields['metadata_json']=json.dumps(metadata)
    out=REPO.parent/'GRUT_particle_transport_checkpoint_20261027.npz';np.savez_compressed(out,**fields)
    info=dict(path=str(out),bytes=out.stat().st_size,sha256=sha(out),cases=tags,s=34.,field_order=metadata['field_order'])
    (REPO.parent/'particle_checkpoint_upload.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info),flush=True)

if __name__=='__main__':main()
