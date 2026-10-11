#!/bin/bash
# ##CHRIS 2026-10-11 (decision 27; 261012 sec. 4.7.56): SUPPORTING COUNTS for Report G2, INFORMATION. Written AFTER the first run of
# the task B script (testG_koa_posthoc_261011.py); not one of task B's tables. Reads the fetched Test G logs only; writes nothing.
#   1. (B5) the contact audit of every run (HD_CONTACT_AUDIT=1): executed physical events and the largest |contact distance| [px]
#      per class, per engine and mass;
#   2. (C1) gen3's event mix per trajectory from its run record ([EDMD3-HEALTH]): the basis of the audit cost estimate;
#   3. (C3) the executed physical events per trajectory, gen3 against gen2 per mass (a static equilibrium property: the collision
#      rates follow from the contact values of the equilibrium distribution).
# usage (from hspist3/): bash experiments_gen3_gate_koa_261009/REPORT_G2_261011/supporting_counts.sh
set -eu
T=experiments_gen3_gate_koa_261009/testG
echo "## 1. contact audit (B5): executed physical events and max |contact distance| [px] per class (24 px = 1 sigma)"
for e in gen2 gen3; do for M in 50 300 1500 2000; do
  grep -h "^\[EDMD-CONTACT\]" "$T/$e/epi8_H_H10_L10/m_$M/run.log" | awk -v e=$e -v M=$M '{gsub(/,/,"");
    for(i=1;i<=NF;i++){if($i=="events")ev+=$(i+1); if($i=="disk-disk"&&$(i+1)>dd)dd=$(i+1); if($i=="divider"&&$(i+1)>dv)dv=$(i+1); if($i=="walls"&&$(i+1)>ow)ow=$(i+1)}; n++}
    END{printf "%s M=%s runs=%d events=%.6g max|c| disk-disk=%s outer-walls=%s divider=%s\n", e, M, n, ev, dd, ow, dv}'
done; done
echo; echo "## 2. gen3 event mix per trajectory (C1 cost basis), means over the runs of [EDMD3-HEALTH]"
for M in 50 300 1500 2000; do
  grep -h "^\[EDMD3-HEALTH\]" "$T/gen3/epi8_H_H10_L10/m_$M/run.log" | awk -v M=$M '{for(i=1;i<=NF;i++){split($i,a,"="); v[a[1]]+=a[2]}; n++}
    END{printf "M=%s runs=%d ev_pair %.0f ev_wall %.0f ev_cross %.0f ev_div %.0f ev_stale %.0f heap_max %.0f run_s %.2f engine_s %.2f (KOA E5-2670 v3)\n",
        M, n, v["ev_pair"]/n, v["ev_wall"]/n, v["ev_cross"]/n, v["ev_div"]/n, v["ev_stale"]/n, v["heap_max"]/n, v["run_s"]/n, v["engine_s"]/n}'
done
echo; echo "## 3. executed physical events per trajectory, gen3 against gen2 (C3)"
for M in 50 300 1500 2000; do
  python3 - "$T" "$M" <<'EOF'
import re, sys, numpy as np
T, M = sys.argv[1], sys.argv[2]
ev = {e: np.array([float(m.group(1)) for m in re.finditer(r"^\[EDMD-CONTACT\] executed events (\d+);",
      open(f"{T}/{e}/epi8_H_H10_L10/m_{M}/run.log", errors="replace").read(), re.M)]) for e in ("gen2", "gen3")}
a, b = ev["gen2"], ev["gen3"]; r = b.mean() / a.mean() - 1
se = np.hypot(a.std(ddof=1) / np.sqrt(len(a)), b.std(ddof=1) / np.sqrt(len(b))) / a.mean()
print(f"M={M}: per trajectory gen2 {a.mean():.1f} (n {len(a)}), gen3 {b.mean():.1f} (n {len(b)}); gen3/gen2 - 1 = {1e4 * r:+.2f} +- {1e4 * se:.2f} x 1e-4 "
      f"(z {r / se:+.2f}); per-run SD {100 * a.std(ddof=1) / a.mean():.3f} %")
EOF
done
