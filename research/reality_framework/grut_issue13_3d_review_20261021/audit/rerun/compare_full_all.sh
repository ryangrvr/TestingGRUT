#!/usr/bin/env bash
# Compare every completed full-size rerun (runs/full_<row>) with the archived row raw/<row>.
SCR=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/rerun; A=/home/user/TestingGRUT/research/reality_framework/grut_issue13_3d_review_20261021/audit/rerun
for f in $SCR/runs/full_*.cmd.json; do n=$(basename $f .cmd.json); row=${n#full_};
 python3 -I $A/compare.py full $SCR/runs/$n $SCR/raw/$row $A/results/rerun_$n.json > /dev/null
 python3 -I -c "
import json;r=json.load(open('$A/results/rerun_$n.json'))
cnt=['fplus_04','fminus_04','fplus_05','fminus_05','fplus_06','fminus_06','genuine_plus','genuine_minus','min_S','max_S']
ok=all(c in r['bitwise_identical_columns'] for c in cnt)
mx=max([v['max_rel'] for k,v in r['differing_columns'].items() if k not in ('mean_S','linear_chi_projection')]+[0])
s=r['rerun_summary']; a=r['archive_summary']
print('$row','counts_bitwise',ok,'topo_identical',r['topology_records_identical'],'max_rel_reduction_diff %.1e'%mx,'a_form',s['a_form'],a['a_form'],'f+- %.6g %.6g'%(s['end_fplus_05'],s['end_fminus_05']),'wind',s.get('end_winding_plus_any'),s.get('end_winding_minus_any'),'dQ re %.2e ar %.2e'%(s['charge_rel_drift_max'],a['charge_rel_drift_max']),'bud re %.4g ar %.4g'%(s['budget_pointwise'],a['budget_pointwise']))
"; done
