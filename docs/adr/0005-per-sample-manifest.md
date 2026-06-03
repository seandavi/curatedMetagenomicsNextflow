# 0005. Per-sample provenance & read-accounting manifest

- **Status:** Accepted
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

The pipeline emits many per-step files (profiles, marker tables, logs) but
nothing machine-readable that captures, per sample: how many reads/bases went
in and survived QC, what was run, and with which versions. Downstream consumers
need read depth and read-length statistics to filter samples and as covariates —
depth is the dominant confounder in nearly every analysis. Today this
information is scattered across the KneadData log and line-count files, and tool
versions live in per-process `versions.yml` files that are never consolidated.

## Decision

Emit a single `manifest.json` at each sample's published root that compiles
everything about the sample: provenance (pipeline/Nextflow version, container,
command line, parameters, input accessions, timestamps), read accounting (raw
and post-decontamination read counts, base counts, and read-length / GC
statistics), the rarefaction parameters when that branch is active, and the
consolidated tool versions.

Specific choices:

- **One file, sample root** — not per-branch. Read accounting describes the real
  sample; the rarefied branch is a downsample and is represented by its
  parameters, not a parallel accounting (see
  [0003](0003-dual-branch-rarefied-profiling.md)).
- **No strict schema alignment** to the consumer's column names. We emit clear,
  stable field names; mapping to the curated resource's schema happens post-hoc.
- **No in-pipeline checksums.** Production output lands in cloud object storage,
  which records per-object MD5/CRC32C for free; computing our own would add
  complexity for no benefit. (Trade-off: local/non-cloud runs get no checksums.)
- **Read statistics computed in pure Python** over the FASTQs (single streaming
  pass), so the manifest step needs no new tool/container — Python is already in
  the base image (see [0001](0001-container-strategy.md)).

## Alternatives considered

- **MultiQC report** — great for human inspection, but produces HTML/aggregate
  artifacts rather than a per-sample machine-readable record keyed to our output
  layout. Complementary, not a substitute; may be added later.
- **Add `seqkit`/`fastp` for stats** — faster, but pulls in another container
  for something a short Python pass does adequately. Rejected per
  [0001](0001-container-strategy.md).
- **Match the consumer schema exactly in-pipeline** — brittle coupling to an
  external schema that changes independently; the consumer maps post-hoc instead.

## Consequences

- Every sample gets a durable, greppable provenance + QC record that downstream
  filtering and normalization can key on directly.
- Per-process `versions.yml` files are consolidated into the manifest, giving
  one place to read the full toolchain version set per run.
- The manifest is a late/finalize-style step that depends on read inputs and the
  collected version files; it does not gate the profiling branches.
- Read-statistics cost is one extra streaming pass over each FASTQ per sample.

## References

- `modules/processes/manifest.nf`, `bin/build_manifest.py`, `main.nf` (wiring).
