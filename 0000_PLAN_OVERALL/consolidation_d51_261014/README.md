# Consolidation of engine-gen3 into main (plan-author decision 51; notes 261012 sec. 4.7.114, Report D51 sec. 4.7.115)

<!-- ##CHRIS 2026-10-14 -->
Prepared on 2026-10-14, read-only: nothing merged, moved or deleted. Executed only on the plan author's go, after B1600 is merged and fetched. The runsheet is in Report D51.

- `d51_dryrun.sh`: the dry run in a scratch clone (git clone --local without a checkout, a sparse checkout of two files, git merge --no-ff --no-commit engine-gen3), with the checks of Report D51 and the runsheet's merge lines rehearsed on a branch of the clone. `d51_dryrun_261014.txt` is its output at main c64d482e / engine-gen3 37748898. On the day it runs again, into a new folder.
- `d51_move_untracked.py`: runsheet step (b). It moves every file the engine-gen3 worktree holds outside git to the same relative path in Repo/HardDisks, and what cannot go there to Repo/HardDisks_d51_hold. It never deletes or overwrites, and by default it only lists. `d51_move_list_261014.txt` is today's list: NOT READY, because 00ALLINONE.c is uncommitted and five processes work in the worktree.
- `d51_worktree_remove_test.sh`: what `git worktree remove` does on this Mac's git (2.50.1) with a populated gitlink, an emptied gitlink folder and an ignored file. `d51_worktree_remove_test_261014.txt` is its output.
- `d51_mover_rehearsal.sh`: the mover end to end (the list, --go, git worktree remove, the checks) on a synthetic layout with every kind of file the worktree holds. `d51_mover_rehearsal_261014.txt` is its output.
