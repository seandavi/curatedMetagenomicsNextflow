# 0012. Resistome profiling with KMA against CARD

- **Status:** Accepted (supersedes [0007](0007-resistome-rgi-card.md))
- **Date:** 2026-07-04
- **Deciders:** Sean Davis

## Context

[ADR-0007](0007-resistome-rgi-card.md) added read-based resistome profiling with
RGI's `bwt` workflow against CARD. In practice RGI is heavyweight for this use:
it requires a multi-step `load` / `card_annotation` / `load` dance to build a
run-local database before every alignment, it is comparatively slow, and its
cost meant we restricted it to the full-depth branch only — leaving the rarefied
branch without a resistome for side-by-side comparison.

We still want an assembly-free, per-sample AMR gene profile against CARD, but
with a simpler, faster tool that we can afford to run on both branches.

## Decision

We will replace RGI with **KMA** (k-mer alignment,
https://bitbucket.org/genomicepidemiology/kma) mapping the host-decontaminated
reads directly against a **KMA-indexed CARD reference**.

- The CARD "broadstreet" release is downloaded once into `store_dir`; its
  homolog-model nucleotide FASTA (`nucleotide_fasta_protein_homolog_model.fasta`)
  is indexed with `kma index` into a separate `store_dir`-backed asset
  (`card_kma_db`) so the index is built once and reused across runs.
- Alignment runs on **both the full and rarefied branches** (KMA is cheap
  enough), matching MetaPhlAn and Kraken, and keeps the branch-aware publishDir
  (`<sample>/<branch>/resistome`).
- KMA runs in its own pinned biocontainer per [0001](0001-container-strategy.md)
  (`quay.io/biocontainers/kma`); indexing runs in the same image (the base image
  has no KMA).
- KMA is invoked with `-ef` (extended features / `.mapstat`) for ARG
  quantification; all plain-text outputs are gzipped before publishing.

## Alternatives considered

- **Keep RGI/CARD (ADR-0007)** — well-recognized, but the multi-step
  load/annotate setup and cost pushed us to full-branch-only. Rejected for
  simplicity and to enable both branches.
- **Build a custom KMA container** — the `kma-adr.md` spec suggested building
  one, but a maintained biocontainer already exists and matches ADR-0001's
  per-process-biocontainer policy, so we avoid owning an image lifecycle.
- **AMR++ / MEGARes, assembly-based AMRFinderPlus** — the same trade-offs noted
  in ADR-0007 still apply; out of scope for this read-based pass.

## Consequences

- **Output contract changes.** The `resistome/` directory now contains KMA
  outputs (`card_kma.res.gz`, `card_kma.mapstat.gz`, `card_kma.aln.gz`,
  `card_kma.fsa.gz`, `card_kma.frag.gz`) instead of the RGI
  `rgi_bwt.*_mapping_data.txt.gz` tables. The pipeline version is bumped to
  2.1.0 accordingly.
- A resistome is now produced for the rarefied branch as well; ARGs are sparse
  at 1M reads, so treat the rarefied resistome as low-sensitivity.
- `rgi_aligner` is removed; `card_db_url` now points at the broadstreet tarball.
- A second `store_dir` asset (`card_kma_db`, the KMA index) is introduced
  alongside the extracted CARD data.
- **Caveat:** the KMA index/align commands are not exercised by stub tests
  (which do not run containers); the first real run is the true validation, as
  with the Kraken and (former) RGI steps.

## References

- KMA: https://bitbucket.org/genomicepidemiology/kma/src/master/
- CARD: https://card.mcmaster.ca/
- Supersedes [0007](0007-resistome-rgi-card.md).
- Implemented in: `modules/processes/resistome.nf` (`resistome_kma`),
  `modules/processes/databases.nf` (`card_db`, `card_kma_db`), `main.nf` wiring,
  `nextflow.config` parameters.
