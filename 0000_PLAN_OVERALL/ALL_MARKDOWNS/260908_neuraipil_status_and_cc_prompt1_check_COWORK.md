# 260908 — CC Prompt 1 check, and NeurAIpil: what is on disk, what is missing for a Cowork ↔ GPT loop (COWORK)

Written 2026-09-08 by Cowork. Sources: your pasted Claude Code output for Prompt 1, and the `neuraipil` repo read from disk (README, CHARTER, ARCHITECTURE, ROADMAP, WORKFLOW, ADR-0003/0004/0005/0009/0011, the A3 and A4 plans, `examples/*.yaml`, `live_tests/*`, `adapters/anthropic/claude_code_cli.py`, `adapters/openai/codex_cli.py`, `adapters/cli_process.py`, `cli/main.py`, `prompting/templates.py`, `config/loader.py`, `core/artifacts/context_bundle.py`, `core/context/model.py`, `core/events/model.py`, `.git/logs/HEAD`). Nothing below is from memory of earlier chats.

---

## 1. CC's Prompt 1 execution, checked against the plan

What matches (and the numbers agree with my 260907 pressure note):

- **TASK 1.** "17 DISCARD lines + PROD line 273" is the same 18 I counted: 17 with health during equilibration, plus the η = 0.60 / N = 900 trajectory that failed in production. "0 accepted with nonzero health" is what the contract promises; good that it was checked rather than assumed.
- **TASK 2.** Calibrate on a *disposable, equilibrated* instance (100 time units at a pre-chunk of min(0.5·320/N, 0.05) — which is 0.05 for every N ≥ 400, so it is just a fixed conservative pre-chunk), then verify the candidate over 20 chunks, halve on failure; production chunk = min(0.6·min(A, B), 320/N) for η ≥ 0.6. That is options (a) + (b) from my note combined, and it is the right combination. But be clear in the methods text about which of the two is the guard: the 320/N cap is the protection, the calibration is advisory (see point 1 below).
- **TASK 4.** Ladder η = 0.69, N = 900, equilibration 400 / 1600 / 6400, 3 seeds, 30 × 20 blocks at chunk 0.3 (0.3 ≤ 320/900 = 0.356, fine), equilibration logging on. The smoke test "Z_pair bit-identical 2.0955759 with and without equil_log" is the right regression: the logging must not touch the RNG stream or the trajectory, and bit identity proves it did not.
- **TASK 3** gated on `.launch_ok` (21 rows, no `CALIBRATION_FAILED`, no chunk above the cap). Correct gating. The disposable-instance seeds 911000001/2 do not collide with any production seed I have seen in the logs (they must not, or the calibration would be a partial replay of a production trajectory).
- Telemetry (16) and seeder (27) regression suites passing before relaunch: as required by the "no core changes during a campaign" rule — the runner changed, the core did not.

Things to fix or watch:

1. **CC's full scan of the avalanche warnings corrects my note.** I wrote (from the first discards I looked at) that the forced advances "sit at t ≈ 7–45 in equilibration". CC's scan says `n=345, min=7.15, median=100.39, max=423.20, 344/345 in equilibration`. So the bursts are spread over the whole 400-unit equilibration with the median at t ≈ 100; they are a property of the dense fluid at those chunk sizes, not a lattice-start transient. Two consequences: the chunk cap ∝ 1/N is the essential fix (a burst can come at any time, and the 250k-event budget is per call), and the new calibration's verification window (t ≈ 100 to ≈ 116 for 20 chunks of ≤ 0.8) still samples only a slice, so a calibration that passes proves little on its own. The 1/345 outside equilibration is the production discard at 0.60/900. I will correct the 260907 note in place.
2. **Mixed-method calibration table.** The pruned `chunk_calibration.csv` keeps 15 rows from the old method (fresh lattice, 5-unit verify) and adds 6 rows from the new one. Production is protected by the cap either way, but the table will end up in the methods section, and a table with two methods and no method column is a provenance hole. Either regenerate all 21 rows with the new method (cheap: one disposable instance of ~120 time units per row) or add a `method` column. Keeping the audit copy `chunk_calibration_20260907_freshseed_method.csv` was the right move; make sure the 6 bogus `CALIBRATION_FAILED` rows from the relative-`cd` bug are in neither file.
3. **Fail closed, not open.** `lookup_chunk` skipping `CALIBRATION_FAILED` must never fall through to a default chunk (0.8) for that cell. For a cell with no valid calibration row the campaign should stop and say so. Please have CC quote the line that guarantees that.
4. **Include the 0.60/900 seed 27388007 rerun in TASK 3.** That cell has 3 accepted seeds at chunk 0.8; a fourth at the capped chunk also tells us whether the chunk size changes Z at all. It should not: the chunk only sets how often the runner interrupts free flight, the collision sequence is the same physics, so Z must agree within block error across chunk sizes. If it moves by more than that, something else is wrong and we want to know before the paper, not after.
5. **What I need from the ladder** to judge the "under-equilibrated" hypothesis: per equilibration length, Z_pair, Z_wall, T, ψ6 global and local (mean over the 3 seeds ± sem), plus one plot of Z(t) with t on a log axis from the equilibration blocks (negative indices) of the 6400 runs, and ψ6 local on the same axis. Decision rule: if Z after 6400 is within the block error of Z after 400 (currently 9.926, −2.66 % vs KR 10.197), then longer equilibration is dead as an explanation, and what remains is the hard-wall state itself or the reference. The next discriminating test would then be a periodic-boundary control, which the pressure runner does not have — that is a Paper 1 scope decision, not a quick run.
6. Keep the machine for the campaign until TASK 3/4 land; Prompt 2 (speed of sound) and Prompt 3 (divider) are queued behind it, not parallel to it.

---

## 2. NeurAIpil: what is actually on disk (verified, not from memory)

- **Git history:** 14 commits, 2026-08-04 → 2026-08-19. Last commit `db72e73` "docs: define A4 execution review boundary" on 2026-08-19 13:15. The A4 contract files (`core/execution/models.py`, `core/review/models.py`, `ports/executor.py`, `ports/reviewer.py`, `adapters/mock/execution.py`, `adapters/mock/review.py`, `core/roles/a4_capabilities.py`, `tests/test_a4_contracts.py`) are dated 13:23–13:30 the same day, i.e. **after** the last commit. So A4 has been started and sits uncommitted, and nothing has happened in the repo for 20 days.
- **A3 is complete and was proven live.** README and ROADMAP: the live acceptance run used Codex `gpt-5.6-sol` and Claude `opus`; 7 provider calls (2 proposals + 2 critiques + 2 revisions + 1 synthesis — the live test asserts exactly 7 receipts); independent proposals, bidirectional critiques, revisions, designated synthesis; it reached `awaiting_plan_approval` with `manual_transfer_count == 0`; 208 offline tests pass; recovery reused a durable Codex receipt instead of repeating the call.
- **What A3 does, concretely:** one task text plus a frozen, content-addressed evidence bundle (UTF-8 text files under one root, limits set in YAML) → each participant proposes without seeing the other → each critiques the other (one target per call) → each revises → one designated participant synthesizes → `runs/<run_id>/run.md` on disk → human plan gate: approve (checksum-bound), bounded edit, or reject → successor run in which the rejection feedback is binding context (ADR-0008). Transports: `codex exec --ephemeral --sandbox read-only --json` with `model_reasoning_effort="high"` hardcoded, and `claude --print --output-format json --tools "" --no-session-persistence --safe-mode --disable-slash-commands --no-chrome --strict-mcp-config --permission-mode dontAsk --model <m>`. API-key variables are stripped from the child processes, so the subscriptions are what gets used. One timeout applies per call: `retry.timeout_seconds` (180 s in the examples).
- **So "let Claude and ChatGPT discuss and reach a solution without me pasting" already exists** — for planning tasks, at the plan gate. It has never been run on a HardDisks question.
- Not implemented, per README/ROADMAP: execution, reviewer, IdeaLabs packets (A5; ADR-0005 defines the schema, no code), NeoDepends/ArchAgent (C-track, deferred), UI, daemon, local models.

---

## 3. What is missing for a physics-review loop (small), and what is not going to happen (Cowork as a participant)

- **Cowork cannot be a participant.** Nothing on your Mac can call this session. Claude the model is in the room through the Claude Code CLI transport. My place is on the human side of the gate, which is where I already work: I read `runs/<id>/run.md` from disk the way I read everything else, and my critique goes back in through the gate — as your bounded plan edit (`neuraipil plan edit --from <file>`) or as the rejection feedback that becomes binding context of the successor run. No paste in either direction. A third participant that is "Cowork's model" via the same CLI is possible (two anthropic participants are allowed by the loader) but adds little independence; I would not.
- **Prompt wording.** Template ids are fixed (`proposal-v1`, `one-target-critique-v1`, `revision-v1`, `synthesis-v1`) and `PromptRenderer` refuses anything else. The proposal request line says "Produce an independent implementation plan." For a physics review that is off, but the task text dominates. Zero-code path: write the task so that "plan" means "decision memo". Clean path, ~50 lines: an `assessment-v1` template family in `prompting/templates.py` plus loader acceptance and tests; no core change.
- **Output and time budgets.** `max_output_tokens: 2000` in the examples is too small for a review; use 8000. 180 s is too short for high reasoning on a 50–100k-token prompt; use 900.
- **Evidence location.** The bundle has one root, and paths cannot traverse `..`. The evidence (`run.log`, `chunk_calibration.csv`, summary tables) lives under `hspist3/...` and the notes under `ALL_MARKDOWNS`. Copy the evidence files into a `context/` folder next to the config. That is a feature: the bundle is frozen and hashed, so `run.md` records exactly which bytes were reviewed.
- **CLI drift.** The Claude Code flags above and `codex exec --json` worked on Aug 19; both CLIs change often. Always start with the existing live smoke test.
- **Usage.** A two-round run is 11 provider calls at high reasoning, each carrying the whole bundle. That is the same ChatGPT subscription whose limit gave the chat export its title. Budget for it, and keep the bundle small (the 26k-line chat export is not evidence).

---

## 4. Recipe: the first real HardDisks deliberation (zero code changes)

Folder (keep it this shallow; the bridge cannot stage files deeper than 7 folders below `HardDisks`):

```
HardDisks/0000_PLAN_OVERALL/neuraipil_reviews/pressure_claim_range_01/
    config.yaml
    task.md
    context/
        evidence/
            260907_pressure_results_and_chat_export_review_COWORK.md
            260908_pressure_final_analysis.md      (CC's, when it exists)
            run.log                                (copied after TASK 3/4 finish)
            chunk_calibration.csv
            chunk_calibration_20260907_freshseed_method.csv
            ladder_0p69_N900_summary.csv           (CC's TASK 4 table)
    runs/                                          (created by neuraipil)
```

The `evidence/` subfolder is not decoration: `_logical_path` refuses a selection that resolves to the context root itself, so `root: context` + `path: .` would fail. A directory selection one level down is the supported shape.

`config.yaml` (checked against `config/loader.py`; every key is required, unknown keys are rejected):

```yaml
schema_version: a3-run-1
run:
  run_id: pressure-claim-range-01
  condition_id: harddisks-pressure-review
  task:
    file: task.md
  rounds: 2
  completion_policy: full-deliberation
storage:
  runs_root: runs
context:
  root: context
  sources:
    - path: evidence
      kind: directory
      recursive: false
      extensions: [.md, .log, .csv]
  limits:
    max_files: 20
    max_file_bytes: 1048576
    max_total_bytes: 4194304
retry:
  max_attempts: 2
  timeout_seconds: 900
  initial_backoff_seconds: 2
  maximum_backoff_seconds: 30
prompts:
  proposal: proposal-v1
  critique: one-target-critique-v1
  revision: revision-v1
  synthesis: synthesis-v1
participants:
  - participant_id: codex-reviewer
    role: solution-architect
    adapter: openai
    transport: codex-cli
    model: gpt-5.6-sol            # the model the Aug 19 acceptance run used; confirm it is still the right id
    capabilities: [planning]
    input_limit_tokens: 150000
    generation:
      max_output_tokens: 8000
  - participant_id: claude-reviewer
    role: critic
    adapter: anthropic
    transport: claude-code-cli
    model: opus                   # as in the acceptance run; use whatever `claude --model` id you want in the room
    capabilities: [planning, synthesis]
    input_limit_tokens: 150000
    generation:
      max_output_tokens: 8000
  - participant_id: executor-intent
    role: executor
    adapter: unbound
    model: unassigned
    capabilities: []
synthesis:
  policy: designated-participant
  participant_id: claude-reviewer
telemetry:
  retain_provider_request_ids: true
  pricing_catalog: null
```

`task.md` (the question; this is the trusted instruction, the evidence blocks are marked untrusted by the renderer):

```
You are reviewing a hard-disk EDMD pressure-validation campaign for a paper.
"Plan" in this workflow means a DECISION MEMO, not an implementation plan.
Nothing is to be executed.

Question: what claim range should Paper 1 make for the equation-of-state
validation (Z vs eta, hard-wall box, N = 400/900/1600, 1/sqrt(N) extrapolation),
and what is the minimal set of additional runs needed to defend it?

Return, in this order:
1. Claim range per eta (validated / onset of departure / exploratory), each
   with the number from the evidence that justifies it, quoted verbatim with
   its file name.
2. Ranked explanations for the eta = 0.69 deficit (Z_inf about 5 % below
   Kolafa-Rottner, deviation growing with N, blocks stationary), each with
   ONE discriminating test and its expected outcome under each explanation.
3. Required vs optional follow-up runs, with estimated cost in wall-clock or
   trajectory count, and what claim each one protects.
4. Explicit disagreements you expect a careful referee to raise.

Rules: quote numbers only from the evidence; if a number you need is not in
the evidence, say so instead of estimating it. Do not reopen settled points
(the accelerated core is excluded; the health contract is not negotiable).
The deliberation is complete when the claim range and the required-run list
are agreed; remaining disagreements are listed, not re-litigated.
```

Commands (from the `neuraipil` repo, in its `.venv`; `neuraipil` is the console script from `pyproject.toml`, or `python -m neuraipil.cli.main`):

```
R=/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL/neuraipil_reviews/pressure_claim_range_01
neuraipil create  --config $R/config.yaml
neuraipil advance pressure-claim-range-01 --config $R/config.yaml      # runs all phases, pauses at the plan gate
neuraipil document pressure-claim-range-01 --runs-root $R/runs --path  # prints the run.md path — that is what I read
neuraipil plan show pressure-claim-range-01 --config $R/config.yaml
# then one of:
neuraipil plan edit    pressure-claim-range-01 --from $R/cowork_edit.md --config $R/config.yaml
neuraipil plan approve pressure-claim-range-01 --expected-document-sha256 <sha> --config $R/config.yaml
neuraipil plan reject  pressure-claim-range-01 --feedback "$(cat $R/cowork_feedback.md)" --config $R/config.yaml
neuraipil replan successor pressure-claim-range-01 --new-run-id pressure-claim-range-02 --operation-id op-02 --config $R/config.yaml
```

`replan successor` re-freezes the context from the same sources, so you can drop new evidence (the ladder table) into `context/` before replanning and the successor run sees it, with the rejection feedback as binding context. That is the loop.

Timing: run this **after** TASK 3/4 land so the bundle contains the ladder result. Before that, only the smoke test (section 6, part A).

---

## 5. IdeaLabs and ArchAgent: sequence, don't combine yet

The roadmap order is A4 (executor + reviewer) → A5 (IdeaLabs work/result packets; ADR-0005 defines them, no code) → C-track (NeoDepends/ArchAgent adapters, explicitly deferred). None of A5 or C exists as code. The order that pays off soonest:

1. Use A3 on one real HardDisks question (section 4). That tests the tool where it hurts and produces the physics decision you need anyway.
2. Finish A4. The contracts are written and uncommitted since Aug 19. A4 is what makes "start building things" real: approved plan → CC implements in an isolated worktree → Codex reviews the diff against the plan → you accept. The pressure-runner changes CC made yesterday (calibration rewrite, `equil_log`) are exactly the kind of change to a measurement instrument during a campaign that should go through that gate.
3. A5 packets when there is an idea that needs deliberation; ArchAgent last — different job, different repo, and it would use the same packet pattern.

One more thing from the repo: `docs/potentialpaper/AI_AUDITLOOP.md` (Aug 5) describes precisely the coordination overhead you are complaining about — "the human repeatedly transported prompts and results between planning chat, the coding environment and independent review" — and proposes a bounded review protocol: classify findings BLOCKER / REQUIRED / OPTIONAL, repair only the first two, verify the named repairs, stop when nothing required remains, expand scope only by explicit human decision. The physics deliberation should carry the same stop rule (it is in the task text above). Without it, two models will happily re-litigate η = 0.69 forever.

---

## 6. Prompt 4 for Claude Code (NeurAIpil; run part A now, part B after Prompt 1 finishes)

```
Context: repo /Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/neuraipil.
Before touching anything, read README.md, CHARTER.md, ARCHITECTURE.md, docs/WORKFLOW.md,
docs/ROADMAP.md and all ADRs under docs/decisions/ (AGENTS.md requires it). Rules from
AGENTS.md apply: narrow scope, no commits, report contradictions instead of inventing
architecture, show git status --short and git diff --stat at the end. Do not touch the
uncommitted A4 files. Do not modify anything under hspist3/.

PART A — plumbing check (no new code)
1. Report `git status --short` (I expect uncommitted A4 contract files from 2026-08-19)
   and `git log --oneline -3`.
2. Report `claude --version` and `codex --version`, and whether both are logged in
   (do not print tokens or credential files).
3. Run the offline suite: `python -m unittest discover -s tests` from the .venv.
   Report the count (README says 208 for A3).
4. I authorize ONE live subscription-backed smoke run now:
   NEURAIPIL_LIVE_TESTS=1 NEURAIPIL_LIVE_ACK_SUBSCRIPTIONS=yes
   NEURAIPIL_CODEX_MODEL=<codex model id> NEURAIPIL_CLAUDE_MODEL=<claude model id>
   python -m unittest live_tests.test_two_subscription_deliberation
   Use the model ids I give you in the reply, not guesses. Report the JSON line it
   prints (status, run_document_path, provider_call_count, manual_transfer_count)
   and, if any CLI flag has been removed since Aug 19, the exact stderr; do not
   patch the adapter without asking.

PART B — first HardDisks deliberation folder (run only after Prompt 1 TASK 3/4 are done)
5. Create HardDisks/0000_PLAN_OVERALL/neuraipil_reviews/pressure_claim_range_01/ with
   config.yaml and task.md exactly as given in
   0000_PLAN_OVERALL/ALL_MARKDOWNS/260908_neuraipil_status_and_cc_prompt1_check_COWORK.md
   section 4, and a context/evidence/ folder containing COPIES of:
   260907_pressure_results_and_chat_export_review_COWORK.md,
   260908_pressure_final_analysis.md, the campaign run.log,
   chunk_calibration.csv, chunk_calibration_20260907_freshseed_method.csv,
   and the TASK 4 ladder summary table. Copies, not links (symlinks are rejected).
6. Run `neuraipil create --config .../config.yaml` and `neuraipil context inspect
   pressure-claim-range-01 --runs-root .../runs`; report the artifact list with byte
   counts and the bundle sha256. If the loader or the bundle freezer rejects anything,
   quote the error and stop; do not change core or adapter code to get past it.
7. STOP there. Do not run `advance` — that spends both subscriptions and I want to
   look at the frozen bundle first.

Report everything with the actual commands and their output; quote numbers from
files, not from memory.
```

Deliberately not in this prompt: the `assessment-v1` template family. Decide after the first run whether "implementation plan" wording actually hurt the output; if it did, that is a 50-line scoped A3.1 change with its own test, and it goes through the same AGENTS.md rules.

---

## 7. For GPT, if you still paste today (5 lines)

```
CC executed Prompt 1: 18 discards confirmed (17 equil + 1 prod), calibration now
runs on an equilibrated disposable instance with a hard production cap
chunk <= 320/N, eta=0.69 N=900 equilibration ladder 400/1600/6400 running (9/9),
TASK 3 reruns gated on the recalibration table. CC's scan of the forced-advance
warnings: n=345, t from 7 to 423, median ~100 -> bursts are a property of the
dense fluid at those chunks, not a start-up transient; the cap is the real fix.
```
