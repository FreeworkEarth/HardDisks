/* PID Demo 1: the window of Demo 0 replaced by a stub, so that demo1/Makefile can build Demo 0 headless from its unchanged
   sources without SDL2 and without make in demo0/ (pre-check 02_DEMO1_BRIEF_PRECHECK.md B12; brief sec. 8). Test 7 runs it
   with --headless, which never calls view_run. */
#include <stdio.h>
#include "pid.h"
int view_run(Run* r){ (void)r; fprintf(stderr, "pid_demo0 (demo1's headless build): no window; run with --headless\n"); return 1; }
