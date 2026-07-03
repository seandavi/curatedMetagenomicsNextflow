# 0011. Isolate per-sample failures with errorStrategy 'ignore'; keep shared processes fail-hard

- **Status:** Accepted (supersedes [0010](0010-error-strategy-finish.md))
- **Date:** 2026-07-03
- **Deciders:** Sean Davis

## Context

[0010](0010-error-strategy-finish.md) switched the non-retryable branch from
`'terminate'` to `'finish'` to stop one bad sample from aborting the whole
batch, and **explicitly rejected `'ignore'`** on the grounds that it "would let
the run report success while silently dropping failed samples."

Production experience (2026-07-02/03) showed the `'finish'` premise was wrong.
`'finish'` lets *already-running* tasks complete but **stops launching new
ones**. With per-sample batching (up to 25 samples per Nextflow run), a failure
in an *early* per-sample process — e.g. `fasterq_dump` exit 3 on a bad/unavailable
SRA accession — strands every batch-mate whose downstream tasks had not launched
yet. Those samples never reach the `MARK_COMPLETE` sentinel and are
dead-lettered. Incident: **4 bad downloads failed 4 batches and dead-lettered 56
samples, 52 of them healthy collateral** (14× amplification). `sacct` confirmed
the driver jobs exited `1:0` — not OOM, not walltime — and the telemetry showed
`MARK_COMPLETE = 0` for entire batches. So `'finish'` isolates only the tasks
that happened to be running at the moment of failure, which at batch size 25 is
a minority of the batch.

Crucially, **0010's objection to `'ignore'` does not hold in this system.**
Per-sample completion is tracked by the `MARK_COMPLETE` sentinel, not by the run
exit code. A sample that is `ignore`d is dropped from its run, never emits
`MARK_COMPLETE`, and is therefore still caught by the orchestrator's close-run
sweep and dead-lettered — visibly. The run exiting `0` hides nothing, because
completion was never inferred from the exit code.

## Decision

Per-sample compute processes use **retry-then-`'ignore'`**:

```groovy
errorStrategy = { (task.exitStatus in 137..140 && task.attempt <= 4) ? 'retry'
                                                                     : (task.attempt <= 2 ? 'retry' : 'ignore') }
```

Retry the OOM / scheduler-kill family up to 4× (memory escalates with
`task.attempt` for those processes), retry any other failure once, then
`'ignore'` — drop *that one sample* and let its batch-mates keep flowing to
`MARK_COMPLETE`.

Processes shared across the batch, or that finalize a sample, stay **fail-hard**
(never `'ignore'`) via `withLabel` overrides that retry transient failures then
`'finish'`:

- **`db_setup`** — every reference/database process in
  `modules/processes/databases.nf`. A broken shared DB must stop the run loudly,
  not silently ruin every sample.
- **`finalize`** — `MARK_COMPLETE` (completion sentinel) and `sample_manifest`
  (publish/manifest). `ignore`ing these would mark a sample complete without its
  outputs, or drop its manifest.

## Alternatives considered

- **`'finish'` (0010).** Rejected: it does not actually isolate — it strands
  not-yet-launched batch-mates, which is most of a 25-sample batch. The problem
  0010 set out to solve is unmet by `'finish'`.
- **`'terminate'` (0009).** Rejected as in 0010 — kills running work immediately.
- **Blanket `'ignore'` on *all* processes (including DB/finalize).** Rejected: a
  systemic failure (broken reference DB, bad container, missing input) would
  silently drop *every* sample while the run reported success. Keeping
  `db_setup`/`finalize` fail-hard makes systemic problems fail loudly; per-sample
  `MARK_COMPLETE` tracking makes per-sample drops visible. The combination is the
  safety net that makes `'ignore'` safe for per-sample processes.
- **One sample per Nextflow run (no batching).** Would eliminate batch-poisoning
  outright but multiplies driver/scheduler overhead ~25×. Not warranted when
  failure isolation solves the same problem.

## Consequences

- A single bad sample drops itself; its batch-mates complete. The 14×
  amplification (4 bad → 56 dead) is eliminated.
- A dropped sample is still caught and dead-lettered (no `MARK_COMPLETE`), so
  failures stay visible — no silent loss. This resolves 0010's concern via the
  sentinel rather than the exit code.
- Systemic failures still fail the run loudly (`db_setup`/`finalize` fail-hard),
  so a broken DB or container surfaces immediately instead of quietly producing
  an empty corpus.
- Genuinely unrecoverable samples (e.g. an unfetchable accession) dead-letter
  *individually* after their retries, instead of taking down batch-mates.
- **Depends on correct labelling.** A per-sample process mislabelled `db_setup`
  or `finalize` would wrongly fail-hard and re-introduce batch poisoning; a
  shared/finalize process left unlabelled would wrongly `ignore`. The semantic
  labels are load-bearing.

## References

- `conf/base.config` — `errorStrategy` default plus `withLabel: db_setup` /
  `withLabel: finalize` overrides.
- Supersedes [0010](0010-error-strategy-finish.md); the retry-family behaviour
  from [0009](0009-retry-and-failure-handling-policy.md) is preserved.
- Diagnosis + rollout: `curatedMetagenomicsNextflow` PR #68, release 2.0.7; the
  batch-poisoning triage in the `nextflow_telemetry` backend.
