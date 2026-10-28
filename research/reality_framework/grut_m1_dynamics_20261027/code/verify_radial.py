"""Independent DOP853 finite-volume evolution and receiver/refinement checks."""
from radial import *
from scipy.integrate import solve_ivp

def main():
    N=12;h=.5;L=N*h;r=(np.arange(N)+.5)*h
    vol=4*np.pi/3*((r+h/2)**3-(r-h/2)**3)
    face=4*np.pi*(r+h/2)**2/h;face[-1]*=2
    matrix=np.zeros((N,N));boundary=np.zeros(N)
    for i in range(N):
        matrix[i,i]-=face[i]/vol[i]
        if i+1<N:matrix[i,i+1]+=face[i]/vol[i]
        else:boundary[i]+=face[i]/vol[i]
        if i:
            matrix[i,i]-=face[i-1]/vol[i];matrix[i,i-1]+=face[i-1]/vol[i]
    q=np.array([.4*np.cos(np.pi*r/(2*L)),.05*np.sin(np.pi*r/L),.8+.2*(r/L)**2])
    p=np.array([.03*np.sin(np.pi*r/L),-.2*q[0],np.zeros(N)])
    def rhs(t,y):
        q=y[:3*N].reshape(3,N);p=y[3*N:].reshape(3,N);x,z,s=q;nn=x*x+z*z
        vp=np.array([(s*s+.1*nn)*x,(s*s+.1*nn)*z,(3*(s*s-1)+nn)*s-.6*(s-1)**2])
        acceleration=q@matrix.T-vp;acceleration[2]+=boundary
        return np.r_[p.ravel(),acceleration.ravel()]
    initial=np.r_[q.ravel(),p.ravel()]
    exact=solve_ivp(rhs,(0,1),initial,method='DOP853',rtol=2e-13,atol=2e-13).y[:,-1]
    src=ROOT/'code/radial.cpp';key=hashlib.sha256(src.read_bytes()).hexdigest()[:16];exe=WORK/('radial_'+key)
    if not exe.exists():raise RuntimeError('Run radial.py first')
    errors=[]
    for dt in (.02,.01,.005):
        d=WORK/('verify_radial_'+str(dt));d.mkdir(exist_ok=True)
        np.savetxt(d/'initial.csv',np.c_[r,q.T,p.T],delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
        subprocess.run([str(exe),str(d/'initial.csv'),str(d/'state'),str(h),str(dt),str(L),'1','2','4','0'],check=True,stdout=subprocess.DEVNULL)
        a=np.loadtxt(d/'state_final.csv',delimiter=',',skiprows=1)
        got=np.r_[a[:,1:4].T.ravel(),a[:,4:7].T.ravel()]
        errors.append(float(np.max(abs(got-exact))))
    ratios=[errors[i]/errors[i+1] for i in range(2)]
    radial=json.loads((ROOT/'RADIAL_RESULTS.json').read_text())['cases']
    ladder=json.loads((ROOT/'EXCITATION_RESULTS.json').read_text())['cases']
    base=radial['radial_main'];fine=radial['radial_fine'];half=radial['radial_half'];boundary=radial['radial_boundary']
    out={'DOP853_errors':errors,'second_order_ratios':ratios,
         'main_vs_fine_core_charge_relative':fine['final_core20_Q']/base['final_core20_Q']-1,
         'main_vs_half_core_charge_relative':half['final_core20_Q']/base['final_core20_Q']-1,
         'main_vs_far_absorber_core_charge_relative':boundary['final_core20_Q']/base['final_core20_Q']-1,
         'receiver_energy_time_refinement_ratio':base['receiver_energy_balance_relative_error']/half['receiver_energy_balance_relative_error'],
         'max_receiver_charge_balance':max(x['receiver_charge_balance_relative_error'] for x in [*radial.values(),*ladder.values()]),
         'max_receiver_energy_balance':max(x['receiver_energy_balance_relative_error'] for x in [*radial.values(),*ladder.values()]),
         'pulse40_negative_emission_relative_grid_difference':ladder['pulse40_fine']['absorbed_Q_over_original']/ladder['pulse40']['absorbed_Q_over_original']-1}
    out['pass']=bool(min(ratios)>3.9 and max(ratios)<4.1 and errors[-1]<1e-4
                     and abs(out['main_vs_fine_core_charge_relative'])<.001
                     and abs(out['main_vs_far_absorber_core_charge_relative'])<1e-5
                     and out['max_receiver_charge_balance']<1e-10 and out['max_receiver_energy_balance']<.001)
    save('RADIAL_CHECKS.json',out);print(json.dumps(out,indent=2))
    if not out['pass']:raise RuntimeError('Radial checks failed')

if __name__=='__main__':main()
