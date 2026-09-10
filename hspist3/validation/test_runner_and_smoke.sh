#!/usr/bin/env bash
# ##CHRIS: test 7 (runner inert-guard) + test 8 (dense smoke tests).
# Test 7 locks out the failure where a trajectory with no collisions at all
# reported Z_pair = 1.000, Z_wall = 0.000 and valid = 1.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
fails=0
ck(){ if [ "$1" = "1" ]; then printf "  [PASS] %s\n" "$2"; else printf "  [FAIL] %s -- %s\n" "$2" "${3:-}"; fails=$((fails+1)); fi; }

echo "7. runner inert-guard"
# A measurement window so short that no collision can occur must NOT be valid,
# even though every health counter is legitimately zero.
out=$(./validation/pressure_validation 0.30 400 20260907 1 0.0000000001 20 \
        /tmp/inert_t.csv /tmp/inert_b.csv 0.5 random 2>&1 | tail -2)
v=$(echo "$out" | grep -o 'valid=[0-9]*' | cut -d= -f2)
inert=$(echo "$out" | grep -c INERT || true)
ck "$([ "${v:-1}" = "0" ] && echo 1 || echo 0)" "zero-collision trajectory is NOT valid" "valid=$v"
ck "$([ "${inert:-0}" -ge 1 ] && echo 1 || echo 0)" "INERT diagnostic is emitted" "no INERT line"

echo
echo "8. dense smoke tests (initializer, geometry, liveness, health)"
printf "  %-22s %-9s %-11s %-10s %-9s %-9s %-8s %-8s %s\n" \
       "config" "init(s)" "min_pair_sep" "min_wall_gap" "psi6_glob" "psi6_loc" "pair_ev" "wall_ev" "health"
for spec in "0.60 400" "0.69 900" "0.72 1600"; do
  set -- $spec; eta=$1; N=$2
  line=$(./validation/lattice_smoke "$eta" "$N" 2>/dev/null)
  [ -z "$line" ] && { ck 0 "eta=$eta N=$N smoke" "no output"; continue; }
  printf "  %s\n" "$line"
  h=$(echo "$line" | awk '{print $NF}')
  pe=$(echo "$line" | awk '{print $(NF-2)}')
  ck "$([ "$h" = "0" ] && [ "$pe" -gt 0 ] && echo 1 || echo 0)" \
     "eta=$eta N=$N: clean and dynamically active" "health=$h pair_ev=$pe"
done

echo
[ "$fails" = "0" ] && echo "ALL TESTS PASSED  (0 failures)" || echo "TESTS FAILED  ($fails failures)"
exit $([ "$fails" = "0" ] && echo 0 || echo 1)
