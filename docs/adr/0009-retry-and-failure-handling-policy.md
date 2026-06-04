# 0009. Retry and failure-handling policy

- **Status:** Accepted (error-action choice superseded by [0010](0010-error-strategy-finish.md); retry policy still in force)
- **Date:** 2026-06-04
- **Deciders:** Sean Davis

## Context

Tasks fail for several distinct reasons on HPC schedulers, and they are not
equivalent:

- **Resource kills** — the scheduler kills a task that exceeds its memory
  request (OOM, exit `137` = 128+SIGKILL) or, less commonly, a related signal
  in the `137..140` range.
- **Wall-time exhaustion** — the scheduler kills a task that exceeds its time
  limit (`TIMEOUT`; surfaces as SIGTERM/`143`, or a null exit code).
- **Cancellation / node failure** — `CANCELLED`, `NODE_FAIL`, preemption; the
  task vanishes with no clean exit code.
- **Application errors** — the tool itself exits non-zero (e.g. `1`) on bad
  input or a genuine bug.

Two further constraints shape the policy:

- **Alpine (and similar ACCESS clusters) cap per-job wall time** — Alpine's QOS
  rejects jobs without an explicit `time`, and the usable ceiling is 24 h.
  Escalating *time* across retries is therefore not available to us there.
- **The run driver is itself a scheduler job.** With `executor = 'slurm'`, the
  Nextflow process runs in its own sbatch allocation and submits each task as a
  separate job. If the *driver* job dies (its own wall-time/memory, a `scancel`,
  or a node failure), every in-flight task is cancelled as collateral and
  appears `ABORTED` — with zero `FAILED` tasks. This failure mode lives *above*
  any per-task `errorStrategy` and cannot be recovered by task-level retries.

Historically the only documented retry behaviour was a config comment. A run
that died from a driver-level failure presented as "many ABORTED, no FAILED,
no .nextflow.log", which is not diagnosable from the telemetry alone and forced
manual `sacct` investigation.

## Decision

We will adopt an explicit, memory-first retry policy in `conf/base.config`,
with these rules:

1. **Escalate memory, not time.** Per-process memory is
   `{ baseline * task.attempt }`; `time` stays flat (24 h on Alpine via
   `conf/profiles/alpine.config`). Most retryable failures here are memory
   pressure, and time escalation is unavailable on the capping clusters anyway.
2. **Retry the resource-kill exit family, terminate otherwise.** The default is
   `errorStrategy = { (task.exitStatus in 137..140 && task.attempt <= 3) ? 'retry' : 'terminate' }`
   with `maxRetries = 3` (downloads get `maxRetries = 4` via the
   `download_retry` label). Walltime/cancel/node-fail/app-error exits do not
   retry.
3. **Resources are declared per process.** Tool processes that are imported
   under aliases for the full/rarefied branches (kraken2, bracken, resistome,
   fastqc, rarefy_fastq) set resources *in-body*; the rest use `withName`
   blocks in `conf/base.config`. (See the comment in
   `modules/processes/kraken.nf`.)
4. **The driver is sized independently of tasks.** Driver memory is set on the
   nf-client wrapper job (outside this repo), because driver death is not a task
   failure and is not covered by `errorStrategy`.

## Alternatives considered

- **Escalate both memory and time on retry** — rejected: Alpine's QOS caps
  wall time at 24 h, so a time bump is either a no-op or a submission rejection.
  Memory is the lever that actually exists.
- **Broaden the retry predicate to all non-zero exits** — rejected for now:
  blindly retrying genuine application errors (exit `1` on bad data) wastes a
  full allocation three times before failing, and hides real bugs. We retry only
  the kill family, where a bigger allocation plausibly fixes the next attempt.
- **`errorStrategy = 'finish'` (or `'ignore'`) instead of `'terminate'`** —
  attractive because `terminate` aborts the entire run (and thus every other
  sample in the batch) on a single unrecoverable task. Deferred, not rejected —
  see Consequences. Changing it warrants its own ADR.

## Consequences

- Retries recover the common case (a task that needs more RAM on the second
  attempt) without operator intervention.
- **Known cost — blast radius.** `terminate` means one task that exhausts its
  retries, or hits a non-`137..140` exit, takes down the whole run; its
  siblings become `ABORTED` collateral. With per-sample batching this can fail a
  batch for one bad sample. Switching to `'finish'` (let running tasks complete,
  stop launching new ones) would cap this. This is the main open question and
  the most likely subject of a follow-up ADR.
- **Walltime and node failures do not retry.** On a cluster where these are
  common, work will terminate rather than recover. Acceptable on Alpine/Anvil
  today; revisit per-site if it bites.
- Driver-level deaths remain outside this policy by construction. They are
  addressed operationally (driver memory sizing) and surfaced via the telemetry
  run-classification (`wrapper-failed` / `ended-no-log`).

## References

- `conf/base.config` (`errorStrategy`, `maxRetries`, per-process `memory`).
- `conf/profiles/alpine.config` (`process.time = '24h'`, executor = slurm).
- `modules/processes/kraken.nf` (rationale for in-body resource directives).
- Driver-death diagnosis and the `wrapper-failed` / `ended-no-log` run
  classification: seandavi/nextflow_telemetry#105.
