squeue -u charing
cd ~/harddisks_gen3_k2/hspist3
ls logs/gen3-preview_*.out | wc -l
grep -lx 'chunk [0-9]* done' logs/gen3-preview_*.out | wc -l
cat logs/gen3-preview_*.out | grep -c -e FAILED -e REFUSED
grep -h '^chunk [0-9]* on ' logs/gen3-preview_*.out | sed 's/: lines.*//; s/^chunk [0-9]* on [^ ]* //' | sort | uniq -c
bash cluster/gen3_preview_261011/preview_stopped_chunks.sh logs gen3-preview 2
