# ##CHRIS 2026-10-09 (261012 sec. 4.7.10, decision 5): summary of trace_paper_inputs.sh's guard trace (the per-script path lists stay
# in scratch: 60k paths). Writes paper_inputs_by_campaign.txt (campaign, files, directories, scripts), the campaign list that
# validation/gen2_defect_scan_261009.py marks as "read by the paper scripts".
# usage (from hspist3/): python3 experiments_gen2_defect_scan_261009/summarize_trace.py <trace dir>
import collections, glob, os, sys
T = sys.argv[1]; HS = os.getcwd(); OUT = "experiments_gen2_defect_scan_261009/paper_inputs_by_campaign.txt"
def campaign(rel):
    p = rel.split("/")
    if p[0] == "experiments_speed_of_sound" and len(p) > 4: return "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/" + p[4]
    if p[0] == "experiments_energy_transfer": return "experiments_energy_transfer/" + p[1]
    return "/".join(p[:2])
c = collections.defaultdict(lambda: [set(), set(), set()]); outside = set()
for f in sorted(glob.glob(os.path.join(T, "paths_*.txt"))):
    s = os.path.basename(f)[6:-4]
    for line in open(f):
        p = line.strip()
        if not p: continue
        if not p.startswith(HS + "/"): outside.add(p); continue
        rel = os.path.relpath(p, HS); k = campaign(rel)
        c[k][0].add(rel); c[k][1].add(rel if os.path.isdir(p) else os.path.dirname(rel)); c[k][2].add(s)
with open(OUT, "w") as o:
    o.write("# campaign | distinct files | distinct directories | paper scripts that read them (run_set.sh names)\n")
    for k in sorted(c):
        o.write(f"{k} | {len(c[k][0])} | {len(c[k][1])} | {','.join(sorted(c[k][2]))}\n")
    o.write(f"# distinct guarded paths outside hspist3 (temporary copies made by the scripts): {len(outside)}\n")
print(open(OUT).read(), end="")
