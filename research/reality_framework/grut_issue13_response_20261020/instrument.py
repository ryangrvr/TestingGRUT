#!/usr/bin/env python3
"""Add read-only field diagnostics to the frozen solver; replay six original cases.

From repository root: python research/reality_framework/grut_issue13_response_20261020/instrument.py
Generated compiler products stay outside the tracked packet. No evolution statement
or original output schedule is changed. Momentum is staggered between original
diagnostics, so the new per-step scan uses positions only.
"""
from pathlib import Path
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

P = Path(__file__).resolve().parent
OLD = P.parent / 'grut_3d_empirical_20261020'
BUILD = Path(os.environ.get('GRUT_WALL_BUILD', '/tmp/grut_issue13_wall_build'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

SCAN = r'''
// Read-only field scan. No momenta, forces, normalization or state writes.
void wall_scan(const Grid& g,const State& b,std::ofstream& out,double t){
 double lo=1e100,hi=-1e100,fneg=0;
 unsigned long long neg=0,n10=0,n25=0,n50=0,n75=0,n90=0,nvac=0,pvac=0;
 #pragma omp parallel for collapse(2) reduction(min:lo) reduction(max:hi) reduction(+:fneg,neg,n10,n25,n50,n75,n90,nvac,pvac) schedule(static)
 for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){
  size_t z=g.id(i,j,k);double S=b.q[2][z],F2=b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z];
  lo=std::min(lo,S);hi=std::max(hi,S);
  if(S<0){neg++;fneg+=F2;}
  n10+=S<-.1;n25+=S<-.25;n50+=S<-.5;n75+=S<-.75;n90+=S<-.9;
  nvac+=(S<-.5&&F2<.1);pvac+=(S>.5&&F2<.1);
 }
 out<<std::setprecision(17)<<t<<','<<lo<<','<<hi<<','<<neg<<','<<n10<<','<<n25<<','<<n50<<','<<n75<<','<<n90<<','<<nvac<<','<<pvac<<','<<(neg?fneg/neg:0)<<'\n';
}
'''

def replace_once(src, old, new):
    if src.count(old) != 1:
        raise RuntimeError('Frozen source hook no longer unique: '+old[:70])
    return src.replace(old, new, 1)

def main():
    protocol = json.loads((P/'PROTOCOL.json').read_text())
    source = OLD/'cartesian.cpp'
    if sha(source) != '3e8c04f2f92c4dc5881344e6f0d52fe16845f551ef07899f3dde3d7da1c632f0':
        raise RuntimeError('Original source identity changed')
    if not (OLD/'raw/production_coarse.execution.json').exists():
        subprocess.run([sys.executable, str(OLD/'archive_io.py'), 'unpack'], check=True)
    BUILD.mkdir(parents=True, exist_ok=True)
    (P/'runs').mkdir(exist_ok=True)
    src = replace_once(source.read_text(), 'int main(int argc,char** argv){try{', SCAN+'\nint main(int argc,char** argv){try{')
    src = replace_once(src, 'auto begin=std::chrono::steady_clock::now();diag(',
                       'std::ofstream walls(path+".walls.csv");walls<<"t,minS,maxS,negative,n_lt_m10,n_lt_m25,n_lt_m50,n_lt_m75,n_lt_m90,negative_lowX,positive_lowX,meanF2_negative\\n";wall_scan(g,b,walls,0);\n auto begin=std::chrono::steady_clock::now();diag(')
    src = replace_once(src, 'g.drift(b,dt);for(auto& t:tt)g.drift(t,dt);',
                       'g.drift(b,dt);for(auto& t:tt)g.drift(t,dt);wall_scan(g,b,walls,n*dt);')
    generated = BUILD/'cartesian_wall_scan.cpp'
    generated.write_text(src)
    exe = BUILD/'cartesian_wall_scan'
    compiler = ['g++','-O3','-std=c++17','-fopenmp','-ffp-contract=off',str(generated),'-o',str(exe)]
    subprocess.run(compiler, check=True)
    environment = {'python':sys.version,'platform':platform.platform(),
                   'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
                   'compiler_command':compiler, 'original_source_sha256':sha(source),
                   'instrument_script_sha256':sha(Path(__file__)),
                   'generated_source_sha256':sha(generated), 'executable_sha256':sha(exe),
                   'protocol_sha256':sha(P/'PROTOCOL.json')}
    (P/'RUN_ENVIRONMENT.json').write_text(json.dumps(environment,indent=2)+'\n')
    for name in protocol['diagnostic_runs']:
        old = json.loads((OLD/'raw'/f'{name}.execution.json').read_text())
        opts = dict(zip(old['command'][1::2],old['command'][2::2]))
        opts.update({'--profile':str(OLD/'profile.dat'),'--modes':str(OLD/'modes.dat'),
                     '--out':str(P/'runs'/f'{name}.csv')})
        cmd = [str(exe)] + [s for kv in opts.items() for s in kv]
        start = time.time()
        result = subprocess.run(cmd,capture_output=True,text=True)
        record = {'case':old['case'],'command':cmd,'returncode':result.returncode,
                  'stdout':result.stdout,'stderr':result.stderr,'wall_seconds':time.time()-start,
                  'instrumented_source_sha256':sha(generated),
                  'original_execution_sha256':sha(OLD/'raw'/f'{name}.execution.json')}
        (P/'runs'/f'{name}.execution.json').write_text(json.dumps(record,indent=2)+'\n')
        print(name,result.returncode,result.stdout.strip(),flush=True)
        if result.returncode:
            raise RuntimeError(result.stderr)

if __name__ == '__main__':
    main()
