# Architecture Decision Records

This directory holds **Architecture Decision Records (ADRs)** — short documents
that capture a single significant decision: the context that forced it, the
choice we made, the alternatives we rejected, and the consequences we accept.

We use ADRs because decisions were previously scattered across GitHub issues,
PR threads, and config comments, where they are hard to find and easy to lose.
An ADR is the durable, greppable home for *why* the pipeline is shaped the way
it is.

## When to write one

Write an ADR when a decision:

- is hard or expensive to reverse (e.g. output directory contract, container
  strategy, taxonomy choice), or
- has non-obvious rationale that a future maintainer would otherwise have to
  reconstruct, or
- was contested — there were real alternatives with real trade-offs.

Do **not** write one for routine, easily-reversible changes (a bug fix, a
dependency bump, a rename). Those belong in the commit message.

## How to write one

1. Copy [`template.md`](template.md) to `NNNN-short-title.md`, where `NNNN` is
   the next zero-padded number in sequence.
2. Fill it in. Keep it short — one decision, one page.
3. Set the status (see below).
4. Add a row to the index below.
5. Commit it alongside the change it documents where possible.

ADRs are immutable once `Accepted`. To change a decision, write a **new** ADR
that supersedes the old one, and update both `Status` lines to point at each
other. This preserves the historical record rather than rewriting it.

## Status values

- **Proposed** — under discussion, not yet acted on.
- **Accepted** — the decision is in force.
- **Accepted (not yet implemented)** — agreed, but the code does not exist yet.
- **Superseded by [NNNN](NNNN-...md)** — replaced by a later decision.
- **Deprecated** — no longer relevant, but kept for the record.

## Index

| ADR | Title | Status |
|-----|-------|--------|
| [0000](0000-record-architecture-decisions.md) | Record architecture decisions | Accepted |
| [0001](0001-container-strategy.md) | Single base image, per-process biocontainers for new tools | Accepted (HUMAnN-in-base-image amended by [0016](0016-humann-bundles.md); self-built images by [0017](0017-self-built-images-on-ghcr.md); base-image build by [0019](0019-base-image-on-ghcr.md)) |
| [0002](0002-defer-humann.md) | Defer HUMAnN functional profiling pending version alignment | Superseded by [0016](0016-humann-bundles.md) |
| [0003](0003-dual-branch-rarefied-profiling.md) | Dual-branch profiling: full depth + rarefied | Accepted |
| [0004](0004-gtdb-conversion-vendored-mapping.md) | GTDB conversion via a vendored, store_dir-backed mapping table | Superseded by [0013](0013-remove-gtdb-conversion.md) |
| [0005](0005-per-sample-manifest.md) | Per-sample provenance & read-accounting manifest | Accepted |
| [0006](0006-kraken2-bracken-complementary-profiler.md) | Complementary read-based profiling with Kraken2 + Bracken | Accepted |
| [0007](0007-resistome-rgi-card.md) | Resistome profiling with RGI against CARD | Superseded by [0012](0012-resistome-kma-card.md) |
| [0008](0008-per-sample-qc-reporting.md) | Per-sample post-decontamination FastQC | Accepted |
| [0009](0009-retry-and-failure-handling-policy.md) | Retry and failure-handling policy | Accepted (error-action superseded by 0010) |
| [0010](0010-error-strategy-finish.md) | Use errorStrategy 'finish' instead of 'terminate' | Accepted |
| [0011](0011-gcs-storage-profile-split.md) | Split the GCS storage profile (publish vs work) | Accepted |
| [0012](0012-resistome-kma-card.md) | Resistome profiling with KMA against CARD | Accepted |
| [0013](0013-remove-gtdb-conversion.md) | Remove in-pipeline GTDB conversion (do it as post-processing) | Accepted |
| [0014](0014-ena-first-read-acquisition.md) | Acquire reads ENA-first with SRA fallback | Accepted |
| [0015](0015-r2-storage-profile.md) | Publish production outputs to Cloudflare R2 (`r2` profile) | Accepted |
| [0016](0016-humann-bundles.md) | HUMAnN as a version-pinned bundle with its own MetaPhlAn pass | Accepted (bundle field list amended by [0018](0018-metaphlan-profiles.md)) |
| [0017](0017-self-built-images-on-ghcr.md) | Publish self-built tool images to GHCR from GitHub Actions | Accepted (extended to the base image by [0019](0019-base-image-on-ghcr.md)) |
| [0018](0018-metaphlan-profiles.md) | Named MetaPhlAn profiles; HUMAnN bundles reference them and reuse the main pass when they match | Accepted (amends [0016](0016-humann-bundles.md); default profile and base image moved by [0019](0019-base-image-on-ghcr.md)) |
| [0019](0019-base-image-on-ghcr.md) | Build the base image in GitHub Actions, publish to GHCR, drop HUMAnN from it | Accepted (amends [0001](0001-container-strategy.md), [0017](0017-self-built-images-on-ghcr.md)) |
