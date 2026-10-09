#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.15, Test G, the A-fixed information row): the registered Method A reduction
(cluster/confinement_20261013/reduce_AF.py, unchanged, called as a subprocess) on an event log that never reaches the disk.
The binary writes its event log (HD_PISTON_EVENTS) into a FIFO; this process reads the FIFO to its end into memory, waits until
the worker has created the marker file (the worker creates it after the binary has EXITED, so the trace tr_<seed>.csv is complete
and closed), and then runs reduce_AF.py with the log on its standard input (/dev/stdin), byte for byte as written.
Why: one A-fixed event log is 17 MB (gzip 6.4 MB); 400 of them do not fit the programme's disk limit (rule 6, 1.5 GB of new data).
The first 10 seeds per engine keep their event log on disk (testG_worker.sh); on those, testG_worker.sh --af-check compares this
path with the file path (red_ rows byte-identical).
usage: af_stream_reduce.py <fifo> <marker> <tr csv> <red csv> <t0> <t1>      exit code = reduce_AF.py's (3 if the marker never came)
"""
import os, subprocess, sys, time
fifo, marker, tr, red, t0, t1 = sys.argv[1:7]
HERE = os.path.dirname(os.path.abspath(__file__))
RAF = os.path.join(os.path.dirname(HERE), "confinement_20261013", "reduce_AF.py")
with open(fifo, "rb") as fh:
    data = fh.read()                                   # until the binary closes the log (at its exit)
t = time.time()
while not os.path.exists(marker):                       # the worker: the binary has exited
    if time.time() - t > 3600: sys.exit(3)
    time.sleep(0.05)
p = subprocess.run([sys.executable, RAF, "/dev/stdin", tr, red, t0, t1], input=data)
print(f"af_stream_reduce: {len(data)} bytes of event log reduced, exit {p.returncode}")
sys.exit(p.returncode)
