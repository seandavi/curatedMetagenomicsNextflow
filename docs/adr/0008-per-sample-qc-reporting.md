# 0008. Per-sample post-decontamination FastQC

- **Status:** Accepted (implemented)
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

Before a full-scale run we want a quick, per-sample quality-control view. Raw-read
FastQC already runs in **both** input modes (embedded in `fasterq_dump` and
`local_fastqc`), so that is not the gap. What is missing is FastQC on the reads
that actually feed profiling — i.e. after KneadData host decontamination and
trimming — so there is no before/after view of what QC did to the reads.

## Decision

Add a **FastQC pass on the host-decontaminated reads** (`fastqc` in
`modules/processes/qc.nf`), published under `<sample>/fastqc/`. It runs in the
base image (FastQC is already present), so **no additional container** is
introduced. Gated by `skip_fastqc` (default false, opt-out). Per-sample only.

## Alternatives considered

- **Add a per-sample MultiQC report** aggregating the raw + post-decontamination
  FastQC. Considered and **rejected for now**: MultiQC is not in the base image,
  and adding another container is not worth it for a per-sample report over just
  two FastQC outputs. The raw and post-decontamination FastQC reports remain
  available individually under the sample directory.
- **Cross-sample MultiQC** — out of scope; this is a per-sample pipeline.
- **Re-run raw FastQC inside this module for symmetry** — wasteful; raw FastQC
  already exists upstream.

## Consequences

- Each sample gets a FastQC profile of the decontaminated reads alongside the
  existing raw-read FastQC, giving a before/after without any new dependency.
- No consolidated single-file QC report; consumers read the two FastQC outputs
  (or aggregate them downstream). Revisiting MultiQC later would mean accepting a
  new container.

## References

- ADR [0001](0001-container-strategy.md) (per-process container).
- Implemented in: `modules/processes/qc.nf` (`fastqc`), `main.nf` wiring,
  `nextflow.config` (`skip_fastqc`).
