# PID Demo 1 -- the one Makefile (brief 02_DEMO1_BRIEF.md sec. 2 and 6; pre-check 02_DEMO1_BRIEF_PRECHECK.md part A).
# Run make in this folder. It writes only into bin/ and runs/ here (both git-ignored); hspist3/edmd_core/ is only read.
#   make                 the box-reference binaries for both sides (fork and unmodified original) and the final-state driver
#   make test            delivery 1a, tests 1-2 (tests/run_tests.sh): a full rebuild, then about 5 min on 4 cores; outputs in
#                        runs/box_reference_<cc>_<UTC time>/. One make test at a time per checkout: bin/ is shared.
#   make test CC=clang   the same with another compiler (one run always compares the two sides under ONE compiler and flag set)
CC       ?= cc
CSTD     := -std=c11 -O2 -Wall -Wextra -ffp-contract=off -Werror
HSP      := ../../hspist3/edmd_core
# 1 under clang (Apple clang included), 0 under gcc: only clang's preprocessor replaces __clang__ by 1.
IS_CLANG := $(shell echo __clang__ | $(CC) -E -P - 2>/dev/null | grep -c '^1$$')

# Flags for code that 1a may not edit, the same on both sides (pre-check A1, A2):
#   gen2, hspist3/edmd_core/edmd.c (read-only; linked unchanged for the harnesses' run D): under -std=c11 glibc declares M_PI
#   only with a feature macro (M_PI is used on 9 lines from edmd.c:556; gcc reports 8), hence -D_DEFAULT_SOURCE (macOS ignores
#   it); gcc also warns at edmd.c:719, 765, 1591, 1592 (-Wmisleading-indentation).
#   gen3, fork and original alike: gcc warns at line 842 of the original (fork line 846: snprintf "D%d" into char kb[8],
#   -Wformat-truncation); the line cannot change in 1a (test 2), so gcc's warning is off for that file on both sides.
GEN2_X   := -D_DEFAULT_SOURCE $(if $(filter 1,$(IS_CLANG)),,-Wno-misleading-indentation)
GEN3_X   := $(if $(filter 1,$(IS_CLANG)),,-Wno-format-truncation)

# The include trap (pre-check A3): the harnesses include "../edmd_gen3.h", "../edmd.h" and "gen3_state_io.h", which the
# compiler looks up next to the harness file before any -I, so a harness compiled in place always sees hspist3's headers.
# The fork side therefore compiles unchanged COPIES of the harness files in bin/hfork/tests/, whose parent bin/hfork/ holds
# copies of the fork's headers; each fork-side link stops if a copy differs from its source or if its dependency list names
# a header of the other side.
MIR      := bin/hfork
HARN     := gen3_m1_harness gen3_m2_harness
FORK_H   := engine/edmd_gen3.h engine/edmd.h
ORIG_H   := $(HSP)/edmd_gen3.h $(HSP)/edmd.h
BINS     := $(addprefix bin/orig/,$(HARN) box_final_state) $(addprefix bin/fork/,$(HARN) box_final_state)
# Rebuild when the compiler, the flags or this Makefile change (make test CC=clang after make test): every object and
# program depends on the stamp (rewritten only when its content changes) and on the Makefile. Quotes are dropped from the
# compiler's version line so that the shell line below cannot break.
STAMP    := bin/.compiler
CC_LINE  := $(subst ',,$(CC) | $(shell $(CC) --version 2>/dev/null | head -1) | $(CSTD) | $(GEN2_X) | $(GEN3_X))
$(shell mkdir -p bin; echo '$(CC_LINE)' | cmp -s - $(STAMP) || echo '$(CC_LINE)' > $(STAMP))
DEPS     := $(STAMP) Makefile

.PHONY: all test
all: $(BINS)

bin/gen2.o: $(HSP)/edmd.c $(HSP)/edmd.h $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(GEN2_X) -c -o $@ $<
bin/orig/edmd_gen3.o: $(HSP)/edmd_gen3.c $(ORIG_H) $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(GEN3_X) -c -o $@ $<
bin/fork/edmd_gen3.o: engine/edmd_gen3.c $(FORK_H) $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) $(GEN3_X) -c -o $@ $<

$(MIR)/edmd_gen3.h: engine/edmd_gen3.h
	@mkdir -p $(@D) && cp $< $@
$(MIR)/edmd.h: engine/edmd.h
	@mkdir -p $(@D) && cp $< $@
$(MIR)/tests/%: $(HSP)/tests/%
	@mkdir -p $(@D) && cp $< $@

$(addprefix bin/orig/,$(HARN)): bin/orig/%: $(HSP)/tests/%.c $(HSP)/tests/gen3_state_io.h $(ORIG_H) bin/orig/edmd_gen3.o bin/gen2.o $(DEPS)
	$(CC) $(CSTD) -MM $< > $@.deps
	@if grep -q -e engine/ -e hfork $@.deps || ! grep -q 'edmd_core/tests/\.\./edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against hspist3's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -o $@ $< bin/orig/edmd_gen3.o bin/gen2.o -lm
$(addprefix bin/fork/,$(HARN)): bin/fork/%: $(MIR)/tests/%.c $(MIR)/tests/gen3_state_io.h $(MIR)/edmd_gen3.h $(MIR)/edmd.h bin/fork/edmd_gen3.o bin/gen2.o $(DEPS)
	@for f in edmd_gen3.h edmd.h; do cmp -s engine/$$f $(MIR)/$$f || { echo "stale copy: $(MIR)/$$f differs from engine/$$f"; exit 1; }; done; \
	 for f in $*.c gen3_state_io.h; do cmp -s $(HSP)/tests/$$f $(MIR)/tests/$$f || { echo "stale copy: $(MIR)/tests/$$f differs from its source"; exit 1; }; done
	$(CC) $(CSTD) -MM $< > $@.deps
	@if grep -q hspist3 $@.deps || ! grep -q 'hfork/tests/\.\./edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -o $@ $< bin/fork/edmd_gen3.o bin/gen2.o -lm

# The final-state driver (tests/box_final_state.c, pre-check A4): one source, built once per side; its "edmd_gen3.h" is found
# through -I only (tests/ holds no such file), and the same dependency check applies.
bin/orig/box_final_state: tests/box_final_state.c $(ORIG_H) bin/orig/edmd_gen3.o $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -I$(HSP) -MM $< > $@.deps
	@if grep -q -e engine/ -e hfork $@.deps || ! grep -q 'edmd_core/edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against hspist3's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -I$(HSP) -o $@ $< bin/orig/edmd_gen3.o -lm
bin/fork/box_final_state: tests/box_final_state.c $(FORK_H) bin/fork/edmd_gen3.o $(DEPS)
	@mkdir -p $(@D)
	$(CC) $(CSTD) -Iengine -MM $< > $@.deps
	@if grep -q hspist3 $@.deps || ! grep -q 'engine/edmd_gen3\.h' $@.deps; then \
	  echo "include check: $@ is not built against the fork's headers (see $@.deps)"; exit 1; fi
	$(CC) $(CSTD) -Iengine -o $@ $< bin/fork/edmd_gen3.o -lm

# make test rebuilds everything first (a few seconds), so the tested programs are always built from the files as they are;
# every command of that build, as make runs it, goes into the run's manifest (pre-check A10).
test:
	@$(MAKE) -s --no-print-directory -B all CC="$(CC)"
	@$(MAKE) -s --no-print-directory -n -B all CC="$(CC)" > bin/build_lines.txt
	@CC="$(CC)" IS_CLANG="$(IS_CLANG)" CSTD="$(CSTD)" sh tests/run_tests.sh
