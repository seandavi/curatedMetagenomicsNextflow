# 0006. Complementary read-based profiling with Kraken2 + Bracken

- **Status:** Accepted (not yet implemented)
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

MetaPhlAn is marker-based: precise, but it only sees taxa represented in its
marker database, and the unclassified ("unknown") fraction can be large. A
k-mer-based profiler (Kraken2 + Bracken) classifies the whole read set, giving a
read-count-based, methodologically independent view — useful for triangulation
and for sample types where the marker view is incomplete.

The cohort spans multiple body sites (gut, airway, skin, vaginal), so the
reference database must be broad rather than gut-specific. The deployment target
includes HPC where the database sits on a **shared filesystem** and **every
process runs on a fresh node** (no opportunity to reuse a node-local staged
copy). Read lengths skew short — many samples are historical (~100 bp).

## Decision

Add Kraken2 + Bracken as an **opt-out** complementary profiler
(`skip_kraken`-style flag), running on the host-decontaminated reads on the same
branches as MetaPhlAn (full and rarefied).

- **Database:** PlusPF (RefSeq bacteria/archaea/viral + protozoa + **fungi** +
  human decoy), **capped to 32 GB** by default, switchable to 64 GB via a
  `kraken_db_url`-style parameter. PlusPF (not gut-specific UHGG) because the
  cohort is multi-body-site; fungi matter for skin (*Malassezia*).
- **Loading:** run `kraken2` in **default mode (load DB into RAM)**, reading the
  DB directly from its staged `store_dir` location — *not* memory-mapping and
  *not* staging to node-local scratch. We have the RAM, and this avoids
  cluster-specific scratch-path logic; it is also the same idiom the existing DB
  inputs use. Provision ~48 GB for the 32 GB DB.
- **Throttle:** cap concurrency with `maxForks` on the Kraken process to prevent
  a batch of fresh nodes from each reading the DB off the shared filesystem
  simultaneously (a "thundering herd" on shared storage).
- **Bracken:** use the prebuilt distribution matching the chosen DB at **100 bp**.

## Alternatives considered

- **UHGG / GTDB gut-specific DB** — highest gut completeness and taxonomy
  coherent with our GTDB output, but prokaryote-only and gut-tuned; it would
  misclassify airway/skin/vaginal samples and miss the mycobiome. Rejected
  because the cohort is multi-site.
- **Full / larger DB (PlusPF full, or 64 GB default)** — better rare-tail
  sensitivity, but doubles per-task shared-FS read and RAM. The gain is in
  low-abundance taxa we trust least (especially in low-biomass skin/airway
  samples), and MetaPhlAn remains the authoritative profiler. Kept as a
  parameter-flip option rather than the default.
- **8 GB capped DB** — too lossy for an authoritative resource. Rejected.
- **Stage-to-scratch + `--memory-mapping`** — decouples RAM from DB size, but
  adds per-cluster scratch handling and a redundant copy for no benefit when RAM
  is available. Rejected.

## Consequences

- Adds a read-count-based, whole-community profile complementing MetaPhlAn's
  relative abundances; helps explain the MetaPhlAn unknown fraction.
- In CPU-hours Kraken2+Bracken is comparable-to-cheaper than MetaPhlAn; the real
  cost is a higher **memory tier** (~48 GB) and **shared-FS read I/O** paid per
  task (no node reuse), mitigated by the cap and the `maxForks` throttle.
- Kraken/Bracken report a different taxonomy than the GTDB output; cross-profiler
  comparison must account for that.
- New tool ⇒ per-process container per [0001](0001-container-strategy.md); DB
  download follows the `store_dir` pattern.

## References

- ADR [0001](0001-container-strategy.md) (per-process container).
- Bracken read-length rationale: predominantly ~100 bp historical reads.
- To be implemented: `modules/processes/kraken.nf`, `databases.nf`
  (`kraken_db`), `main.nf` wiring, `conf/base.config` resource tier + `maxForks`.
