#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mode=${1:-checks}
if [[ "$mode" == full ]]; then
  python code/model.py
  python code/trajectory.py --tag holder
  python code/trajectory.py --tag reference --original --ks none
  python code/trajectory.py --tag refined --rtol 2e-10
  python code/trajectory.py --tag finite_k --ks .001,.0001
  python code/perturbations.py
  python code/elastic.py --nodes 128 --angles 16 --source 96
  python code/elastic.py --nodes 256 --angles 16 --source 96
  python code/elastic.py --nodes 256 --angles 32 --source 96
  python code/elastic.py --nodes 512 --angles 16 --source 96
  python code/elastic.py --nodes 256 --angles 16 --source 64
  python code/linear_response.py --tag linear
  python code/linear_response.py --tag linear_refined --nq 48 --L 220 --rtol 2e-8 --switch 300
  python code/linear_response.py --tag linear_more --ks 100000,300000,1000000 --nq 96 --L 260 --rtol 1e-8 --switch 450
  python code/linear_response.py --tag linear_collisionless --collision 0 --nq 48 --L 220 --rtol 2e-8 --switch 300
  python code/linear_response.py --tag linear_isothermal --seed 0 --nq 48 --L 220 --rtol 2e-8 --switch 300
  python code/linear_response.py --tag linear_fit --shape fit --nq 48 --L 220 --rtol 2e-8 --switch 300
  python code/linear_response.py --tag linear_fit_refined --shape fit --ks 10000,30000,100000,300000,1000000 --nq 96 --L 260 --rtol 1e-8 --switch 450
  python code/linear_response.py --tag linear_rate_half --shape fit --ks 30000,100000,300000,1000000 --collision .5 --nq 48 --L 220 --rtol 2e-8 --switch 300
  python code/linear_response.py --tag linear_rate_double --shape fit --ks 30000,100000,300000,1000000 --collision 2 --nq 48 --L 220 --rtol 2e-8 --switch 300
  python code/linear_response.py --tag linear_replay_refined --shape fit --ks 10000 --nq 96 --L 260 --rtol 1e-8 --switch 450
elif [[ "$mode" == class ]]; then
  : "${GRUT_CLASS_BINARY:?Set GRUT_CLASS_BINARY to the pinned CLASS executable}"
  python code/class_control.py --binary "$GRUT_CLASS_BINARY" --tag coldmicro --kmax 100000
  python code/class_control.py --binary "$GRUT_CLASS_BINARY" --tag warmmicro --warm --kmax 100000
  python code/class_control.py --binary "$GRUT_CLASS_BINARY" --tag coldrk20 --kmax 20 --evolver 0
  python code/class_control.py --binary "$GRUT_CLASS_BINARY" --tag warmrk20 --warm --kmax 20 --nq 60 --evolver 0
  python code/class_control.py --binary "$GRUT_CLASS_BINARY" --tag coldrk_refined --kmax 20 --tol 1e-6 --evolver 0
  python code/class_control.py --binary "$GRUT_CLASS_BINARY" --tag warmrk120 --warm --kmax 20 --nq 120 --tol 1e-6 --evolver 0
elif [[ "$mode" != checks ]]; then
  printf 'Usage: bash reproduce.sh [checks|full|class]\n' >&2
  exit 2
fi
python code/verify.py
python code/notes_check.py
python code/summarize.py
python code/package.py
