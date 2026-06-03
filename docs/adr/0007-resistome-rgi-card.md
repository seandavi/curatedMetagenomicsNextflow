# 0007. Resistome profiling with RGI against CARD

- **Status:** Accepted (not yet implemented)
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

## References

- CARD / RGI: https://card.mcmaster.ca/
- ADR [0001](0001-container-strategy.md) (per-process container), `store_dir`
  database pattern.
- To be implemented: `modules/processes/resistome.nf`, `databases.nf` (CARD
  download), `main.nf` wiring.
