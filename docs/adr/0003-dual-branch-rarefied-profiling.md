# 0003. Dual-branch profiling: full depth + rarefied

- **Status:** Accepted
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

Sequencing depth is a dominant, study-specific confounder in metagenomic
profiling: richness and the detectability of low-abundance taxa scale with read
count, which complicates comparison across samples and studies sequenced at
different depths. A curated resource that mixes many studies needs a way to
support depth-controlled comparisons.

## Decision

By default (`skip_rarefied = false`) run profiling **twice**: once on the full
host-decontaminated reads, and once on reads downsampled with `seqtk` to a fixed
depth (`rarefy_reads`, default 1,000,000) using a fixed seed (`rarefy_seed`,
default 42) for reproducibility. The two runs publish into branch-specific
subdirectories (`full_data/` and `rarefied_data/`) so results can be compared
side by side. Setting `skip_rarefied = true` restores the original single-branch
layout with no `full_data/` prefix, preserving backward compatibility.

## Alternatives considered

- **Full depth only** — simplest, but pushes all depth normalization onto
  downstream users, who then need the raw reads to redo it. Rejected for a
  resource meant to be analysis-ready.
- **Rarefied only** — discards information from deeply sequenced samples.
  Rejected.
- **Rarefy after profiling (subsample the profile)** — statistically inferior to
  rarefying reads before alignment for marker-based profilers. Rejected.

## Consequences

- Roughly doubles the profiling compute per sample when enabled (two MetaPhlAn /
  marker passes), accepted as the cost of depth-controlled comparability.
- Output layout gains a branch dimension (`full_data/` vs `rarefied_data/`) that
  every downstream profiling step (MetaPhlAn lists/markers, StrainPhlAn, GTDB)
  must respect via `meta.branch`.
- Per-sample **read accounting describes the real sample**, not the rarefied
  copy; the rarefaction target/seed are recorded as parameters rather than as a
  parallel accounting (see [0005](0005-per-sample-manifest.md)).

## References

- `modules/processes/rarefaction.nf`, `main.nf` (branch wiring), `nextflow.config`
  (`skip_rarefied`, `rarefy_reads`, `rarefy_seed`).
- Commit `82a6d68` (dual profiling branches with seqtk rarefaction).
