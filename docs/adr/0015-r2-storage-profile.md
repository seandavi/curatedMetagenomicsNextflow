# 0015. Publish production outputs to Cloudflare R2 (`r2` profile)

- **Status:** Accepted
- **Date:** 2026-10-03
- **Deciders:** Sean Davis

## Context

Production runs `-profile alpine,gcs`, which publishes to
`gs://cmgd-data/results/cMDv4/` ([ADR-0011](0011-gcs-storage-profile-split.md)
split that profile so it only moves outputs). nextflow_telemetry ADR 0008
moved all cloud object storage for these projects to Cloudflare R2: downloads
to onclappc02 and to outside users cost nothing in egress, and the control
plane and DuckLake data already live there. Every sample is being re-run
under readset ids (nextflow_telemetry ADR 0007), so existing GCS outputs are
not migrated.

Two constraints showed up when testing on Alpine:

- Alpine's `nextflow` module is 24.04.4. Its `nf-amazon` (AWS SDK v1) reads the
  bucket's location hint from R2 (`WNAM`) and sends the write to
  `s3.WNAM.amazonaws.com`, ignoring `aws.client.endpoint`, so every publish
  fails with "Failed to create publish directory"
  ([nextflow-io/nextflow#4873](https://github.com/nextflow-io/nextflow/issues/4873),
  fixed by #5779 in 25.04). No setting of `aws.region` avoids it.
- Nextflow 26.04 enables the strict syntax parser by default, and this
  pipeline fails `nextflow lint` under it (180 errors).

## Decision

We will add an `r2` storage profile and make `alpine,r2` the production
profile.

- `r2` sets `publish_base_dir = 's3://cmgd-raw'` and `publish_mode = 'copy'`,
  and points `aws.client.endpoint` at
  `https://$R2_ACCOUNT_ID.r2.cloudflarestorage.com` with path-style access and
  region `us-east-1` (R2 treats it as `auto`). Outputs land at
  `s3://cmgd-raw/cmgd_nextflow/<version>/<sample_id>/<step>/`.
- Keys come from `R2_ACCOUNT_ID`, `R2_ACCESS_KEY_ID`, `R2_SECRET_ACCESS_KEY`
  in the driver's environment, never from the repo. Non-`AWS_*` names keep the
  awscli calls inside task containers from picking them up.
- The `r2` profile needs Nextflow >= 25.04. Alpine runs 25.10.8 from a
  standalone launcher on Java 18 (nextflow_telemetry `config/nf_tel.env.alpine`),
  not the 24.04 module.
- `gcs` stays for now and is removed once the v2 re-run has replaced the GCS
  outputs.

## Alternatives considered

- **Upgrade to Nextflow 26.04.** Rejected for now: the pipeline does not pass
  the strict parser. Running 26.04 with `NXF_SYNTAX_PARSER=v1` would work but
  adds a second override for no benefit over 25.10.
- **Stay on Nextflow 24.04 and work around the endpoint bug.** Rejected: no
  `aws.region` value avoids the regional redirect (tried `auto`, `us-east-1`,
  `WNAM`).
- **Publish locally and sync to R2 with rclone after the run.** Rejected: a
  second copy step outside Nextflow, with its own failure modes and no
  `publishDir` semantics.

## Consequences

- A test run of sample `05f281407e15a03e65eba0dd74f30fae` (ERR866581) with
  `alpine,r2` on Nextflow 25.10.8 published the same 128 files as its 2.2.1 GCS
  output. Every MetaPhlAn and Kraken/Bracken profile is identical after
  decompression. The other differences are read order, KMA's date header, and
  FastQC's sampled duplicate estimate.
- Anvil cannot use `r2` yet: it only offers Java 8/11, which caps Nextflow at
  23.10.1. Anvil stays on `anvil,gcs` until it has Java 17.
- The manifest's `nextflow_version` changes from 24.04.4 to 25.10.8 for new
  Alpine runs.

## References

- nextflow_telemetry ADR 0008 (object storage on R2), ADR 0007 (readsets)
- [ADR-0011](0011-gcs-storage-profile-split.md)
- Issue #79
