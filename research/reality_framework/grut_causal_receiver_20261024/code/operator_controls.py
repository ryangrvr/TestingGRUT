#!/usr/bin/env python3
import io,json,subprocess,tempfile
from pathlib import Path
import numpy as np
from scipy.linalg import eigvals
ROOT=Path(__file__).resolve().parents[1]


def main():
    with tempfile.TemporaryDirectory(prefix='grut-cfs-operator-') as tmp:
        binary=Path(tmp)/'controls'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'code/operator_controls.cpp'),'-o',str(binary)],check=True)
        local=json.loads(subprocess.check_output([str(binary)],text=True));rows=[]
        for count,h in [(64,1.),(128,.5)]:
            for mass in [1.,np.sqrt(6.)]:
                raw=subprocess.check_output([str(binary),'matrix',str(count),str(h),str(mass)],text=True)
                matrix=np.loadtxt(io.StringIO(raw),delimiter=',');spectrum=eigvals(matrix)
                for dt in [.0625,.03125,.015625]:
                    x=dt*spectrum;amplification=1+x+x*x/2+x**3/6+x**4/24
                    rows.append(dict(nodes=count-1,h=h,mass=float(mass),dt=dt,
                        state_dimension=len(matrix),maximum_real_eigenvalue=float(max(spectrum.real)),
                        maximum_RK4_eigen_amplification=float(max(abs(amplification)))))
    out=dict(scope='Actual C++ RHS energy/charge and gradient checks; complete finite-grid linear eigenspectra at two representative resolutions. This is not a nonlinear or all-resolution stability theorem.',
        local=local,spectral_controls=rows)
    (ROOT/'OPERATOR_CONTROLS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':main()
