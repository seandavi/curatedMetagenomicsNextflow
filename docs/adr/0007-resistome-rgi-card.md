# 0007. Resistome profiling with RGI against CARD

- **Status:** Accepted (implemented)
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

Antimicrobial-resistance (AMR) gene content is one of the most-cited reuses of
shotgun metagenomes, and the pipeline currently produces none. We want a
per-sample resistome profile without requiring assembly, to keep it cheap and
applicable to shallow/historical samples.

## Decision

Add a read-based resistome step using **RGI (Resistance Gene Identifier) against
the CARD database**, emitting a per-sample AMR gene abundance table (gene plus
AMR family / drug class), assembly-free. Like the other reference databases,
CARD is downloaded/cached under `store_dir`, and RGI runs in its own pinned
container per [0001](0001-container-strategy.md).

## Alternatives considered

- **AMR++ / MEGARes** — purpose-built for read-based resistome *quantification*
  with a hierarchical ontology; a reasonable alternative. Chosen RGI/CARD instead
  for CARD's curation and broad recognition in the gut/clinical literature, which
  favors interpretability and comparability for a curated resource.
- **Assembly + contig-based ARG calling (e.g. AMRFinderPlus)** — higher
  precision and gives genomic context (co-localization, mobility), but requires
  per-sample assembly, which is a separate, heavier scope decision we have not
  taken. Rejected for the read-based first pass.

## Consequences

- Adds a widely-requested AMR data product at modest, assembly-free cost.
- Read-based calls lack genomic context (no contig/mobility information); that is
  accepted for this first pass and could be revisited if assembly is later added.
- Virulence-factor profiling (e.g. VFDB) is a natural sibling and could reuse the
  same read-based pattern, but is out of scope here.

## Update (implementation, 2026-06-03)

- **Full-depth branch only.** Unlike MetaPhlAn/GTDB/Kraken (which run on both
  the full and rarefied branches), the resistome runs only on the full-depth
  reads: ARGs are rare and the rarefied 1M-read depth yields a sparse,
  low-value resistome. It still uses the conventional **branch-aware**
  publishDir (`<sample>/<branch>/resistome`), so with the rarefied branch active
  it publishes under `full_data/resistome/`.
- **Opt-out (default on)** via `skip_resistome`, consistent with `skip_gtdb` /
  `skip_kraken`.
- **RGI `bwt` read-based workflow.** `rgi load --local` + `rgi card_annotation`
  + `rgi bwt --read_one ... --local`. The reads are single-end (the pipeline is
  unpaired throughout), so only `--read_one` is supplied. `rgi_aligner`
  (default `bowtie2`) selects the bwt aligner. Container:
  `quay.io/biocontainers/rgi:6.0.5` (set in the process body, alias-irrelevant
  here but consistent with the Kraken processes).
- **Caveat:** the RGI command sequence follows CARD documentation but is **not
  exercised by stub tests** (which do not run containers) — the first real run
  is the true validation, as with the Kraken DB/container specifics. The exact
  biocontainer build suffix should be confirmed against current quay.io tags.

## References

- CARD / RGI: https://card.mcmaster.ca/
- ADR [0001](0001-container-strategy.md) (per-process container), `store_dir`
  database pattern.
- Implemented in: `modules/processes/resistome.nf`,
  `modules/processes/databases.nf` (`card_db`), `main.nf` wiring,
  `nextflow.config` parameters.
