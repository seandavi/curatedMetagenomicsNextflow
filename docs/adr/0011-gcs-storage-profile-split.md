# 0011. Split the GCS storage profile (publish vs work)

- **Status:** Accepted
- **Date:** 2026-06-04
- **Deciders:** Sean Davis

## Context

The `gcs` profile bundled two independent choices: publishing pipeline
**outputs** to GCS *and* routing the Nextflow **work** directory to GCS. That
bundling made it unusable on HPC: composing `-profile alpine,gcs` to get
durable GCS output also forced `workDir = gs://cmgd-data/work`, which the SLURM
executor cannot use without a GCS-backed work filesystem (fusion/gcsfuse). The
practical effect was that nobody enabled `gcs` on Alpine, so outputs went to the
default `${launchDir}/results` on scratch — and were deleted with the work dir
at the end of each job. Nothing durable was produced, and the per-task
`.command*` files (a useful debugging/progress view) were lost.

The two concerns are separable: GCS *publish* only needs the driver to copy
results to `gs://` (it has `GOOGLE_APPLICATION_CREDENTIALS`); GCS *workDir* only
makes sense for cloud compute that reads/writes work natively from `gs://`.

## Decision

We will split the profile in two:

- **`gcs`** — publish outputs to GCS (`publish_base_dir = gs://…`,
  `publish_mode = 'copy'`, `google.project`). No `workDir` override. Safe to
  compose with HPC: `-profile alpine,gcs` = local scratch workDir, GCS output.
- **`gcswork`** — only sets `workDir = gs://cmgd-data/work`. For cloud compute:
  `-profile google,gcs,gcswork`.

HPC production uses `-profile alpine,gcs`.

## Alternatives considered

- **Keep the combined profile** — rejected: couples output durability to a
  cloud workDir that doesn't work on SLURM, so it stayed off and outputs were
  ephemeral.
- **Set GCS workDir on HPC via gcsfuse/fusion** — deferred: more moving parts on
  the compute nodes than warranted; local scratch workDir is fine, only the
  results need to persist.

## Consequences

- `-profile alpine,gcs` now publishes durable outputs (including each process's
  `.command*`) to `gs://cmgd-data/results/cMDv<cmgd_version>/<name>/<version>/…`
  while keeping a fast local scratch workDir.
- Validation gating: the driver's publishDir copy to `gs://` depends on
  Nextflow's GCS filesystem support + credentials being present in the job — the
  first real `alpine,gcs` run is the true test (cf. ADR-0007's note).
- Cloud-compute runs must now add `gcswork` explicitly to get GCS work staging.

## References

- `conf/profiles/gcs.config`.
- Output path derivation: `nextflow.config` `publish_base_dir`.
