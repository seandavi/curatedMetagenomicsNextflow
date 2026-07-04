# 0013. Remove in-pipeline GTDB conversion (do it as post-processing)

- **Status:** Accepted (supersedes [0004](0004-gtdb-conversion-vendored-mapping.md))
- **Date:** 2026-07-04
- **Deciders:** Sean Davis

## Context

[ADR-0004](0004-gtdb-conversion-vendored-mapping.md) added an in-pipeline step
(`metaphlan_to_gtdb`) that translated each MetaPhlAn SGB profile into a
GTDB-taxonomy profile, backed by a `store_dir`-cached SGB→GTDB assignment table
(`sgb_to_gtdb_db`) and a vendored conversion script. It ran on both branches and
published a `gtdb/` subdirectory per branch.

In practice the conversion is just a **static relational mapping** applied to the
already-published MetaPhlAn profiles: substitute 1:1 SGB→GTDB mappings, sum n:1
mappings into the shared taxon, aggregate up the lineage. Nothing about it needs
the pipeline's compute context, the container, or per-run execution — it can be
done once, downstream, over the published profiles.

Keeping it in the pipeline carried ongoing cost: a vendored script that collided
with the container's bundled copy (see CHANGELOG 2.0.6), an `sgb2gtdb_url` that
had to be kept in lockstep with `metaphlan_index`, an extra DB-setup process, and
two aliased compute steps plus their completion-gating wiring.

## Decision

We will **remove the GTDB conversion from the pipeline entirely** and perform it
as a post-processing step over the published MetaPhlAn profiles instead. Removed:
the `gtdb.nf` module (`metaphlan_to_gtdb`), the `sgb_to_gtdb_db` database
process, the vendored `bin/cmgd_sgb_to_gtdb.py`, the `skip_gtdb` / `sgb2gtdb_url`
parameters, the `skip_gtdb` field in the manifest, and all related config and
docs.

## Alternatives considered

- **Keep it, gated by `skip_gtdb`** — status quo. Rejected: it still carries the
  script-collision, version-lockstep, and wiring maintenance for something with
  no per-run compute dependency.
- **Move it to a single post-run process at the end of the pipeline** — would
  still live in-repo and re-introduce the container/version coupling. Rejected in
  favor of true out-of-band post-processing over the published `mpa` profiles.

## Consequences

- **Output-contract change:** the per-branch `gtdb/` subdirectory and
  `gtdb_profile.tsv.gz` are no longer produced. Pipeline version bumped to
  **2.2.0**.
- **Manifest change:** the `parameters.skip_gtdb` field is dropped from
  `manifest.json`.
- The SGB→GTDB translation now lives downstream; the mapping table
  (MetaPhlAn's `*_SGB2GTDB.tsv`, matched to the active index) must be applied
  there. The published MetaPhlAn profiles carry everything needed to do so.
- Fewer moving parts in the pipeline: one less DB download, one less vendored
  script, no `metaphlan_index`/`sgb2gtdb_url` lockstep to maintain.

## References

- Supersedes [0004](0004-gtdb-conversion-vendored-mapping.md).
- GTDB: https://gtdb.ecogenomic.org/
- Removed from: `modules/processes/gtdb.nf`, `modules/processes/databases.nf`
  (`sgb_to_gtdb_db`), `bin/cmgd_sgb_to_gtdb.py`, `main.nf` wiring,
  `nextflow.config`, `conf/base.config`, `modules/processes/manifest.nf`,
  `bin/build_manifest.py`.
