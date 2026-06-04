# 0010. Use errorStrategy 'finish' instead of 'terminate'

- **Status:** Accepted (supersedes the error-action choice in [0009](0009-retry-and-failure-handling-policy.md))
- **Date:** 2026-06-04
- **Deciders:** Sean Davis

## Context

[0009](0009-retry-and-failure-handling-policy.md) set the retry policy and used
`'terminate'` for non-retryable failures, while explicitly flagging
`'finish'` as deferred-not-rejected. Production experience settled it: with
`executor = 'slurm'` and per-sample batching, a single unrecoverable task (e.g.
the RGI image that failed to pull) made Nextflow **terminate the whole run** —
cancelling every in-flight sibling task, which then surfaced as a wave of
`ABORTED` tasks and failed/dead-lettered jobs for samples that had nothing wrong
with them. One bad sample failed the entire batch.

## Decision

We will use **`errorStrategy = 'finish'`** for the non-retryable branch (the
retry-the-137..140-kill-family behaviour from 0009 is unchanged):

```groovy
errorStrategy = { (task.exitStatus in 137..140 && task.attempt <= 3) ? 'retry' : 'finish' }
```

`'finish'` lets tasks already running complete and stops launching new ones,
then ends the run with a non-zero status — instead of `'terminate'`, which kills
running tasks immediately.

## Alternatives considered

- **`'terminate'`** (the 0009 status quo) — rejected: one unrecoverable task
  wastes all concurrent work in the batch and inflates the failure/DLQ counts
  with collateral, obscuring the real cause.
- **`'ignore'`** — rejected: it would let the run report success while silently
  dropping failed samples. We want the run to end non-zero so the orchestrator
  sweeps the genuinely-incomplete jobs, just without nuking healthy siblings.

## Consequences

- A single bad sample no longer aborts the batch; sibling tasks finish and their
  outputs/telemetry are preserved.
- Runs with a failed task still end non-zero, so the telemetry sweep
  (close-run → retry/DLQ) still routes the incomplete jobs correctly — now
  scoped to the actually-failed samples rather than the whole batch.
- Slightly longer wind-down on failure (running tasks drain rather than being
  killed). Accepted — the wasted-work cost of `terminate` was larger.

## References

- `conf/base.config` (`errorStrategy`).
- Supersedes the error-action choice in [0009](0009-retry-and-failure-handling-policy.md).
- Trigger case: RGI container tag fix (manifest unknown → terminate cascade).
