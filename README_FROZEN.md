# FROZEN: the engine-gen3 worktree, kept read-only

<!-- ##CHRIS 2026-10-15 (plan-author decisions 51 and 55, addenda 1 and 2) -->
This folder, Repo/HardDisks_gen3, is the git worktree of the branch engine-gen3. It was frozen on 2026-10-15
after plan-author decisions 51 and 55 (addenda 1 and 2): engine-gen3 was merged into main (merge commit 0544b7e8),
every untracked and ignored file was moved to Repo/HardDisks (what could not go there went to
Repo/HardDisks_d51_hold), and the worktree was locked:

    git worktree lock --reason "frozen after decision 51/55; work on main in Repo/HardDisks" ../HardDisks_gen3

Work on main in Repo/HardDisks: all fetches, builds and reports go there. Nothing in this folder is changed
or deleted; this file is the one tracked change made after the freeze.

It stays locked for one week. On 2026-10-22, and only on Chris's word, from Repo/HardDisks:

    git worktree unlock ../HardDisks_gen3
    git worktree remove ../HardDisks_gen3

The record: Repo/HardDisks/0000_PLAN_OVERALL/ALL_MARKDOWNS/261012_paper1_confinement.md, sec. 4.7.114
(decision 51), 4.7.115 (Report D51, the runsheet), 4.7.127 (decision 55), 4.7.129 and 4.7.131 (its addenda 1
and 2), 4.7.130 (step b with the folder kept) and 4.7.132 (the merge, the checks and this file).
