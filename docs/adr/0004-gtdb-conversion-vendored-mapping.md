# 0004. GTDB conversion via a vendored, store_dir-backed mapping table

- **Status:** Superseded by [0013](0013-remove-gtdb-conversion.md)
- **Date:** 2026-06-03
- **Deciders:** Sean Davis

## Context

MetaPhlAn 4 reports taxonomy using its own SGB clade names. Many downstream
consumers want GTDB taxonomy. MetaPhlAn ships an official converter,
`metaphlan/utils/sgb_to_gtdb_profile.py`, plus a versioned SGB→GTDB assignment
table (`mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv`). The mapping is not purely
1:1 — some SGBs collapse n:1 onto a shared GTDB taxon, which the official script
handles by binning and lineage aggregation.

Two friction points with using the upstream script as-is:

1. It hard-codes the assignment table path relative to the installed MetaPhlAn
   package. We want the pipeline to own and pin that table (download once, cache,
   reproduce), not depend on whatever copy happens to be baked into the container
   image — which may be older than, or drift from, our pinned `metaphlan_index`.
2. It imports MetaPhlAn package internals (`util_fun`), so it is awkward to run
   standalone.

## Decision

**Vendor** a self-contained adaptation of the official script as
`bin/sgb_to_gtdb_profile.py` (Nextflow stages `bin/` onto `PATH`), with two
deliberate, minimal changes: (a) the assignment table is supplied via
`-m/--mapping` instead of being hard-coded, and (b) the `info`/`error` helpers
are inlined to drop the package import. The **conversion logic is unchanged**
from upstream (1:1 substitution, n:1 binning, lineage aggregation). The table
itself is downloaded once into `store_dir` by a dedicated `sgb_to_gtdb_db`
process (URL is the `sgb2gtdb_url` parameter) and passed in explicitly.

## Alternatives considered

- **Call the script bundled in the container** — couples us to whatever version
  the image ships and to its hard-coded table path; defeats pinning the mapping
  to our `metaphlan_index`. Rejected.
- **Reimplement the conversion from scratch** — risks diverging from the
  validated upstream binning/aggregation behavior. Rejected in favor of a
  faithful, minimally-patched copy.

## Consequences

- The SGB→GTDB table is a first-class, pinned, cached pipeline input; conversion
  output is reproducible and decoupled from the container image's bundled copy.
- We carry a vendored script that must be re-synced if upstream changes its
  conversion logic. The file header records its upstream provenance to make that
  sync deliberate.
- `sgb2gtdb_url` must be kept in lockstep with `metaphlan_index`; this coupling
  is documented in `nextflow.config` and the README.

## References

- PR #52 (implementation); closes issue #18.
- `bin/sgb_to_gtdb_profile.py`, `modules/processes/gtdb.nf`,
  `modules/processes/databases.nf` (`sgb_to_gtdb_db`).
- Upstream: https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/sgb_to_gtdb_profile.py
- Relies on the container strategy in [0001](0001-container-strategy.md).
