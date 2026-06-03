# 0008. Per-sample QC reporting: post-decontamination FastQC + MultiQC

- **Status:** Accepted (implemented)
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

Before a full-scale run we want a quick, per-sample quality-control view. Raw-read
FastQC already runs in **both** input modes (embedded in `fasterq_dump` and
`local_fastqc`), so that is not the gap. What is missing:

1. No FastQC on the reads that actually feed profiling — i.e. after
   KneadData host decontamination and trimming. Without it there is no
   before/after view of what QC did to the reads.
2. No consolidated, human-readable per-sample QC report; the FastQC outputs and
   logs are scattered across subdirectories.

## Decision

Add two steps, both **per-sample only** (no cross-sample aggregation — that is a
separate downstream concern and explicitly out of scope here):

1. A **FastQC pass on the host-decontaminated reads** (`fastqc` in
   `modules/processes/qc.nf`), published under `<sample>/fastqc/`.
2. A **per-sample MultiQC** (`multiqc`) that aggregates the raw and
   post-decontamination FastQC into one `multiqc_report.html` (+ data) under
   `<sample>/multiqc/`. The two FastQC reports are staged into `raw/` and
   `clean/` subdirectories and MultiQC is run with `--dirs` so the otherwise
   identical FastQC sample names are disambiguated.

Both are gated by `skip_multiqc` (default false, opt-out, consistent with the
other optional steps). FastQC runs in the base image; MultiQC uses a pinned
biocontainer (per [0001](0001-container-strategy.md)).

## Alternatives considered

- **Cross-sample MultiQC** (the usual way MultiQC is used) — explicitly out of
  scope: this is a per-sample pipeline; any across-sample reporting is handled
  later, downstream.
- **Re-run raw FastQC inside the QC module** for symmetry — wasteful; the raw
  FastQC already exists, so MultiQC consumes it directly rather than recomputing.
- **Feed Kraken2/Bracken/MetaPhlAn into MultiQC now** — deferred. v1 keeps the
  report to the FastQC before/after to stay low-risk; the MultiQC input set can
  be widened later (MultiQC has Kraken modules).

## Consequences

- Each sample gets a before/after read-QC report in one HTML, plus a FastQC
  profile of the decontaminated reads.
- New MultiQC container; its exact biocontainer tag is **not exercised by stub
  tests** (no containers run there), so confirm it on the first real run — same
  caveat class as the Kraken/RGI containers.
- The MultiQC input set is intentionally minimal for now; widening it (Kraken,
  etc.) is future work, not a new architectural decision.

## References

- ADR [0001](0001-container-strategy.md) (per-process container).
- Implemented in: `modules/processes/qc.nf` (`fastqc`, `multiqc`), `main.nf`
  wiring, `nextflow.config` (`skip_multiqc`).
