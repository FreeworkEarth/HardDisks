squeue -u charing
sacct -X -u charing -S 2026-10-10 --name=gen3-testG --format=JobID%20,State,Elapsed,ExitCode | grep -v ' COMPLETED '
ls ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | wc -l
grep -lx 'chunk [0-9]* done' ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | wc -l
cat ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | grep -c -e FAILED -e REFUSED -e STOP
grep -h '^chunk [0-9]* on ' ~/harddisks_gen3/hspist3/logs/gen3-testG_*.out | sed 's/: lines.*//; s/^chunk [0-9]* on [^ ]* //' | sort | uniq -c
