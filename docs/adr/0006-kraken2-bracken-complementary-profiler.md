# 0006. Complementary read-based profiling with Kraken2 + Bracken

- **Status:** Accepted (implemented)
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

## Update (implementation, 2026-06-03)

Reconciling the decision with what is actually available and practical:

- **DB cap is 16 GB, not 32 GB.** The prebuilt PlusPF indexes
  (https://benlangmead.github.io/aws-indexes/k2) are offered only at 8 GB,
  16 GB, or full — there is no 32 GB build, contrary to the original Decision.
  The default `kraken_db_url` therefore points at the **16 GB** PlusPF cap
  (release `20260226`); the parameter lets a user select the full (better
  rare-tail sensitivity) or 8 GB build. The 16 GB cap is the closest available
  match to the intent and loads comfortably into the provisioned RAM.
- **Two processes, two StaPH-B containers.** Implemented as `kraken2`
  (`staphb/kraken2:2.1.3`) feeding `bracken` (`staphb/bracken:2.9`) rather than
  one combined container, because a single image confirmed to carry both tools
  could not be verified. At abundance-estimation time Bracken needs only the DB
  (which ships the kmer distributions) and the Kraken2 report, so it does not
  need Kraken2 in its container.
- **Directives in the process body, not `conf/base.config`.** Because the
  processes are imported under aliases (`*_full`/`*_rarefied`), container, cpus,
  memory (with `task.attempt` scaling), and `maxForks` are set in the process
  bodies so they apply regardless of how selectors match aliases.
- Memory tier is ~32 GB for the 16 GB DB (loaded into RAM); bump it if pointing
  at the full DB.

## References

- ADR [0001](0001-container-strategy.md) (per-process container).
- Bracken read-length rationale: predominantly ~100 bp historical reads.
- Implemented in: `modules/processes/kraken.nf`, `modules/processes/databases.nf`
  (`kraken_db`), `main.nf` wiring, `nextflow.config` parameters.
