# HardDisks
Investigating Complex System Interactions with Hard Disks

HSPIST3


New preset: --experiment=wall_mid

One internal wall exactly at the mid-point (1/2 width).
Default split is 50/50 unless overridden by the particle-box options below.
Per-compartment particle allocation via CLI

For 2 boxes (one wall): use either
--particles-box-left=COUNT and/or --particles-box-right=COUNT
or the generic list: --particles-boxes=left,right
For N boxes (N = walls + 1): use the list form left→right
--particles-boxes=c1,c2,...,cN
You can give counts or weights; they’re normalized and leftover is resolved so the total always matches (no manual math needed).
Wall position for the “two-wall preset, variant 1” remains one-third from the left (leaving two-thirds to the right). The mid-wall preset gives you the single wall in the center explicitly.

How to use











# Test sigmoid protocol with energy measurement
./00ALLINONE --no-experiments --protocol=sigmoidal --energy-measurement

# Test multi-wall energy transfer
./00ALLINONE --no-experiments --num-walls=3 --wall-positions="0.0,5.0,10.0" --energy-measurement --protocol=sigmoidal

# Compare protocols
./00ALLINONE --no-experiments --protocol=step --energy-measurement
./00ALLINONE --no-experiments --protocol=sigmoidal --energy-measurement
./00ALLINONE --no-experiments --protocol=linear --energy-measurement

./00ALLINONE --no-experiments --protocol=sigmoidal --energy-measurement





Single mid wall, 50/50 split:
./00ALLINONE --experiment=wall_mid --no-experiments

Single mid wall, 25/75 split (100 total particles example):
./00ALLINONE --experiment=wall_mid --particles-box-left=25 --particles-box-right=75 --no-experiments

Generic list form (also works with >2 boxes):
./00ALLINONE --experiment=wall_mid --particles-boxes=25,75 --no-experiments

Your original example with 2/3 space on the right and custom 100/200 counts (for 300 particles):
./00ALLINONE --experiment=2wall_exp_with_2_wall:1 --particles-boxes=100,200 --no-experiments

Notes

The allocator is applied in both placement paths (grid/honeycomb and random), so your chosen split appears at t=0.
If you pass only one of left/right, the other side is inferred from the remainder.
Counts/weights are normalized; you don’t need the exact total.
If the single-wall variant still appears on the wrong side in your current preset run, switch to wall_mid for an explicit center wall, or keep 2wall_exp_with_2_wall:1 and set --particles-boxes as above; the divider placement logic is now consistent, and the allocation will respect your counts.