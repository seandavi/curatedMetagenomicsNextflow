# 0002. Defer HUMAnN functional profiling pending version alignment

- **Status:** Accepted
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

The pipeline contains a complete HUMAnN functional-profiling branch (gene
families, pathways, EC abundances) wired through `modules/processes/humann.nf`
and the ChocoPhlAn/UniRef/utility-mapping database downloads. Functional
potential is one of the most-requested data products for a curated metagenomic
resource.

However, the branch is disabled by default (`skip_humann = true`). The active
MetaPhlAn index (`mpa_vJan25_CHOCOPhlAnSGB_202503`) is newer than the
MetaPhlAn/HUMAnN version combination this HUMAnN path can consume. Running the
branch as-is produces incompatible or incorrect results, so it is not safe to
enable without coordinated version alignment and validation. Today this rationale
lives only as a comment in `nextflow.config`.

## Decision

Keep the HUMAnN branch in the codebase but **disabled by default**
(`skip_humann = true`). Treat it as preserved-but-dormant. Re-enabling requires
deliberately realigning the MetaPhlAn index with a HUMAnN-compatible version
combination (the older `mpa_vOct22_CHOCOPhlAnSGB_202403` index is noted in config
as the HUMAnN4-compatible option) and validating outputs before flipping the
default.

## Alternatives considered

- **Remove the HUMAnN branch** — sheds dead code, but discards working plumbing
  for a high-value data product we intend to restore. Rejected.
- **Enable it now** — produces wrong results given the version mismatch.
  Rejected.
- **Pin MetaPhlAn back to the HUMAnN-compatible index globally** — would
  downgrade the taxonomic profiling everyone else depends on just to satisfy
  HUMAnN. Rejected; the taxonomic index is the higher-priority contract.

## Consequences

- No functional profiles are published by default; downstream users get
  taxonomy only until this is resolved.
- The branch must be kept compiling and stub-testable so it does not rot while
  dormant.
- Re-enabling is a scoped, validated task (version alignment), not new
  architecture — tracked as future work.

## References

- `nextflow.config` (`skip_humann`, `metaphlan_index`, and the
  HUMAnN-compatible index comment).
- `modules/processes/humann.nf`.
