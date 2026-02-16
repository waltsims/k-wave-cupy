## Bug workflow (mandatory)

When the user reports a bug, do **not** start by proposing fixes.

### 1) Reproduce-first gate
Before any fix attempt:
- Create or identify a **minimal reproducible test** that fails and demonstrates the bug.
- Prefer automated tests (unit/integration/regression). If automation is not feasible, write a deterministic repro script and expected output.
- Record: environment assumptions, exact command(s) to run, and expected vs actual behavior.

**Hard rule:** No code changes are allowed until a failing repro exists.

### 2) Subagent fix attempts (after repro exists)
Once a failing test exists, spawn subagents to attempt fixes in parallel:
- Each subagent must propose a fix **and** explain how it will make the repro test pass.
- Each subagent must provide either:
  - a patch + evidence the test would pass, or
  - a reasoning-based proof tied directly to the failing test and change.

### 3) Merge criteria
Accept a fix only if:
- The original repro test now passes.
- No related tests regress (or at least, the agent explains any new failures).
- The repro test is kept as a **regression test** (do not delete it).

### 4) Output format when a bug is reported
Always respond in this order:
1. “Repro/Test Plan” (what test you’ll add/run)
2. “Current hypothesis” (optional, but no fixes yet)
3. “Subagent tasks” (what each subagent will try after repro is in place)
4. “Fix + Proof” (only after the failing test exists)