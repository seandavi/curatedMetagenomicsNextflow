# 0001. Single base image, per-process biocontainers for new tools

- **Status:** Accepted
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

Historically the pipeline runs every process inside one container image,
`seandavi/curatedmetagenomics:metaphlan4.2.2`, pinned in `conf/base.config` and
built manually via `docker/cloudbuild.yaml`. This was simple and kept the
toolchain (KneadData, MetaPhlAn, HUMAnN, bowtie2, …) versioned together.

As we add tools that are unrelated to MetaPhlAn — GTDB conversion, a per-sample
manifest, Kraken2/Bracken, and a resistome caller (RGI/CARD) — the single-image
approach starts to hurt: each new tool forces a rebuild of the one image,
couples unrelated version upgrades, and grows the image toward a monolith whose
dependency graph is hard to keep consistent.

## Decision

Keep the single base image for the **existing** MetaPhlAn/KneadData/HUMAnN
toolchain, but introduce **per-process containers** (pinned biocontainers) for
**newly added, independent tools** via Nextflow's `withName` container override.
Prefer pure-stdlib implementations where a new tool would otherwise be pulled in
only for something trivial (e.g. the manifest's read statistics are computed in
Python, which the base image already provides, rather than adding a `seqkit`
container).

## Alternatives considered

- **Bake every new tool into the base image** — keeps "one container," but
  produces a monolithic image, couples unrelated upgrades, and makes every
  addition a full rebuild-and-repush of a large image.
- **Move fully to per-process biocontainers immediately** — cleanest long-term
  nf-core-style model, but a large, risky migration of the stable existing
  toolchain for no immediate benefit.

## Consequences

- New tools get independently pinned, smaller images and decoupled version
  management; the stable MetaPhlAn image stops being a bottleneck.
- The pipeline now mixes container provenance models. Each process must be
  explicit about its container, and the base image remains the default.
- A future ADR may migrate the existing toolchain to per-process containers if
  the mixed model proves awkward; this decision deliberately defers that.

## References

- `conf/base.config` (`process.container` default and `withName` overrides).
- `docker/cloudbuild.yaml` (base image build).
- ADRs [0004](0004-gtdb-conversion-vendored-mapping.md),
  [0005](0005-per-sample-manifest.md),
  [0006](0006-kraken2-bracken-complementary-profiler.md),
  [0007](0007-resistome-rgi-card.md) all rely on this strategy.
