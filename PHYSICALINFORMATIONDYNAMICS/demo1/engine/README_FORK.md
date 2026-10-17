# demo1/engine: the PID fork of gen3 (Demo 1, brief `02_DEMO1_BRIEF.md`)

- **Source.** `edmd_gen3.c`, `edmd_gen3.h` and `edmd.h` were copied from `hspist3/edmd_core/` at origin/main `6b29aa5b`. Each file's 4-line header gives that file's last change (`a71f394d` for gen3, `942c3b63` for `edmd.h`) and its source blob: `edmd_gen3.c` `e7d9ce3e`, `edmd_gen3.h` `6616f999`, `edmd.h` `4ccd5da0`.
- **1a.** Each file is its original plus the 4-line header, so original line n = fork line n + 4. `make test` (test 2) checks this byte for byte, and git confirms the recorded commits and blobs.
- **gen2 is not forked.** The M1/M2 harnesses' run D needs `hspist3/edmd_core/edmd.c` (blob `b9921b0f`), so both sides of the box reference link the same unmodified object, `bin/gen2.o`. Brief §2 rule 8 is read as: fork gen3 and original gen3 are never in one binary.
- **Compiler exemptions** (`Makefile`; pre-check `02_DEMO1_BRIEF_PRECHECK.md` A1, A2), identical on both sides. gcc only: `-Wno-format-truncation` for gen3 (original line 842, fork line 846, which 1a may not edit). gen2: `-D_DEFAULT_SOURCE`, plus `-Wno-misleading-indentation` under gcc.
- **Test 1 coverage.** It compares the fork with the original byte for byte: the M1/M2 harnesses, plus `tests/box_final_state.c`. The second covers final states, records, the collision hook, checkpoint restart, the event log and the impulse/virial accessors. It cannot see a detector that has been switched off, so 1b needs positive controls (a white-box test that plants an error the audit must report).
- **1b** changes `edmd_gen3.c` and `.h` only (brief §4). Every change site carries a comment that opens with `/* PID: `.
  - New public-struct fields go at the end, and `EDMD3_NCLS` stays 5 (pre-check B8), so test 1's record prefixes stay comparable.
  - The hook and the checkpoints follow brief §4 "Also" and §4.8, with pre-check B7 decided there.
  - Test 3 is test 1's full byte comparison with an empty primitive list, not only its hashes.
- **`edmd.h` stays byte-identical below its header,** because the unmodified `edmd.c` is compiled against the original. Dividers, pistons and bands are kept: `pid_demo1` will not use them, but the box reference exercises them.
