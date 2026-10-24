#!/usr/bin/env python3
"""Add an exact discrete extraction interface without changing field dynamics."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def build():
    code=(ROOT/'code/radial_cfs.cpp').read_text()
    def replace(old,new):
        nonlocal code
        assert code.count(old)==1,(old,code.count(old));code=code.replace(old,new)
    replace('double er=0,qr=0;','double er=0,qr=0,ef=0,qf=0;')
    replace('int n;double h,','int probe_cut=0;\nint n;double h,')
    replace('d.er=0;d.qr=0;','d.er=0;d.qr=0;d.ef=0;d.qf=0;')
    old='\n}\nvoid stage(State&out'
    flux='''
 // Exact currents for the symmetric seven-point Hamiltonian, across a node cut.
 const std::array<double,3>coupling={87./64,-3./32,1./576};
 for(int i=std::max(1,probe_cut-2);i<=probe_cut;i++)for(int ell=1;ell<=3;ell++){
  int j=i+ell;if(j<=probe_cut||j>=n)continue;
  double coefficient=4*PI/h*coupling[ell-1];
  for(int f=0;f<3;f++)d.ef+=.5*coefficient*(u.a[f][i]*u.a[3+f][j]-u.a[3+f][i]*u.a[f][j]);
  d.qf+=coefficient*(u.a[0][i]*u.a[1][j]-u.a[1][i]*u.a[0][j]);
 }
}
void stage(State&out'''
    replace(old,flux)
    replace('u.qr+=dt/6*(k1.qr+2*k2.qr+2*k3.qr+k4.qr);',
        'u.qr+=dt/6*(k1.qr+2*k2.qr+2*k3.qr+k4.qr);\n u.ef+=dt/6*(k1.ef+2*k2.ef+2*k3.ef+k4.ef);u.qf+=dt/6*(k1.qf+2*k2.qf+2*k3.qf+k4.qf);')
    replace('double E=0,Q=0,core=0,massdelta=0,Fmax=0;','double E=0,Q=0,core=0,massdelta=0,Fmax=0,Eprobe=0,Qprobe=0;')
    replace('b.E+=weight*e;b.Q+=weight*q;',
        'b.E+=weight*e;b.Q+=weight*q;if(i<=probe_cut){b.Eprobe+=weight*e;b.Qprobe+=weight*q;}')
    replace('Q60,max_layer_mass2_shift,Fmax\\n','Q60,max_layer_mass2_shift,Fmax,E_inside_cut,Q_inside_cut,E_flux_cut,Q_flux_cut\\n')
    replace('auto b0=budget(u);','probe_cut=ip;\n auto b0=budget(u);')
    replace("<<b.massdelta<<','<<b.Fmax<<'\\n';", "<<b.massdelta<<','<<b.Fmax<<','<<b.Eprobe<<','<<b.Qprobe<<','<<u.ef<<','<<u.qf<<'\\n';")
    target=ROOT/'code/radial_flux.cpp';target.write_text(code)
    return target


if __name__=='__main__':print(build())
