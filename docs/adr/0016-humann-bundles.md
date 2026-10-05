# 0016. HUMAnN as a version-pinned bundle with its own MetaPhlAn pass

- **Status:** Accepted (supersedes [0002](0002-defer-humann.md); amends [0001](0001-container-strategy.md))
- **Date:** 2026-10-05
- **Deciders:** Sean Davis

## Context

[ADR-0002](0002-defer-humann.md) kept HUMAnN dormant (`skip_humann = true`)
because the pipeline's MetaPhlAn index (`mpa_vJan25_CHOCOPhlAnSGB_202503`,
MetaPhlAn 4.2.2) is newer than anything HUMAnN can consume, and it left
"realign the versions" as unscoped future work. The facts that scope it
(verified 2026-10-05, issue #82):

- **HUMAnN 3.9 requires the `vJun23` MetaPhlAn database.**
  `humann/config.py@v3.9` sets `metaphlan_v4_db_version="vJun23"`, and
  `humann/search/prescreen.py` (lines ~151, 207-212) exits when the supplied
  taxonomic profile's header lacks it. HUMAnN 3.9 also cannot parse MetaPhlAn
  4.2.x output at all, so the HUMAnN pass needs MetaPhlAn 4.1.x.
- **HUMAnN 4.0.0a1 needs `vOct22_CHOCOPhlAnSGB_202403`** and names its outputs
  differently (`_2_genefamilies`, `_3_reactions`, `_4_pathabundance`,
  `_5_pathcoverage`; unmapped row `READS_UNMAPPED` instead of `UNMAPPED`).
- The databases `humann_databases` fetches are addressed by short names
  (`full`, `uniref90_ec_filtered_diamond`) whose content depends on the HUMAnN
  release, so a cache keyed by name alone is ambiguous across releases (see the
  cache convention in `modules/processes/databases.nf`, #83).
- [ADR-0001](0001-container-strategy.md) kept HUMAnN in the single base image.
  Version alignment between MetaPhlAn and HUMAnN cannot be achieved there
  without changing the main taxonomy pass.

The published MetaPhlAn 4.2.2 / vJan25 taxonomy must not change to satisfy
HUMAnN.

## Decision

HUMAnN runs as a self-contained, **bundle**-driven subworkflow on the
full-depth branch only, beside (not instead of) the main taxonomy pass.

- **Bundle.** A named, tested set of pins defined in one registry,
  `conf/humann_bundles.config` (`params.humann_bundles`): HUMAnN container,
  MetaPhlAn container, MetaPhlAn version and index, and the ChocoPhlAn / UniRef /
  utility-mapping database names. `params.humann_bundle` selects one (default
  `humann3.9`). There are deliberately no free-form per-tool version parameters:
  compatibility is a single upstream pin per HUMAnN release, and a new release is
  one more registry entry. An unknown name fails at workflow start with the list
  of valid names, but only when `!skip_humann`, so a stale value never blocks
  runs that skip HUMAnN.
- **`humann3.9` bundle.** HUMAnN
  `quay.io/biocontainers/humann:3.9--py312hdfd78af_0`; MetaPhlAn
  `quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0` with the
  `mpa_vJun23_CHOCOPhlAnSGB_202307` index (present on the MetaPhlAn database
  server and recognised by MetaPhlAn 4.1.1's `--index`); ChocoPhlAn `full`,
  UniRef `uniref90_ec_filtered_diamond`, utility mapping `full` (all three
  confirmed by `humann_databases --available` in the 3.9 container, so the
  previous defaults carry over).
- **Why a separate MetaPhlAn container.** The HUMAnN 3.9 biocontainer also ships
  MetaPhlAn 4.1.1, but its `bowtie2` is 2.2.3 (the biocontainer's MetaPhlAn
  image has 2.5.4). The MetaPhlAn pass is the step that maps whole samples
  against the vJun23 index, so it runs in the dedicated MetaPhlAn 4.1.1 image;
  HUMAnN keeps its own image.
- **Processes.** `metaphlan_for_humann` (bundle MetaPhlAn container; input the
  host-decontaminated full-depth reads; `-t rel_ab_w_read_stats`) produces the
  `--taxonomic-profile` for `humann` (bundle HUMAnN container). Containers,
  resources and `maxForks` are set in the process bodies as for
  Kraken2 ([ADR-0006](0006-kraken2-bracken-complementary-profiler.md)); the
  container is a closure resolved from the bundle when a task launches. `humann`
  renormalizes (cpm, relab), splits stratified tables and gzips; its outputs are
  globs of HUMAnN-native filenames, and renormalization/splitting loops over
  whatever tables HUMAnN wrote, so a bundle with a different file set (HUMAnN 4)
  needs no process change.
- **Databases.** A bundle-pinned MetaPhlAn index process and the three HUMAnN
  database processes are invoked from `DATABASES` when `!skip_humann`
  (so `--databases_only --skip_humann false` stages them) and run in the bundle's
  containers. Cache keys: MetaPhlAn index → `metaphlan/<index>`; ChocoPhlAn,
  UniRef, utility mapping → `<db>/<bundle>`. Bundle, not database name, is the
  version, because the same name can denote different content across releases.
- **Publish layout.** `<sample>/humann/<bundle>/` for the HUMAnN tables, with
  the profile that drove the stratification under
  `<sample>/humann/<bundle>/metaphlan/` as provenance: functional profiles are
  stratified by this bundle's taxonomy, not by the published MetaPhlAn 4.2.2
  profile.
- **Manifest.** `manifest.json` keeps the main pass under `metaphlan` / `bowtie2`
  and records the HUMAnN pass under distinct keys (`metaphlan_humann`,
  `bowtie2_humann`, `humann`), because `merge_versions()` in
  `bin/build_manifest.py` is a flat last-file-wins map. When HUMAnN ran, the
  manifest `parameters` also carry `humann_bundle` and `humann_metaphlan_index`.
- **Default stays off.** `skip_humann` remains `true`. Enabling it by default is
  a separate decision (#90) gated on real-sample validation and cost data (#88)
  and on the choice of default bundle (#86).

## Alternatives considered

- **Downgrade the main MetaPhlAn index to a HUMAnN-compatible one** — would
  degrade the taxonomy every consumer depends on (and HUMAnN 3.9 still could not
  read 4.2.x output). Rejected.
- **Derive HUMAnN's profile from the main pass** — profiles from MetaPhlAn 4.2.2
  / vJan25 lack the `vJun23` header and species set HUMAnN needs. Rejected.
- **Free-form `humann_version` / `metaphlan_version` / index parameters** —
  admits untested combinations; compatibility is one pin per release. Rejected
  for a registry of named bundles.
- **Use the HUMAnN container's bundled MetaPhlAn** — fewer images, but its
  bowtie2 (2.2.3) is much older than the MetaPhlAn image's (2.5.4) for the most
  expensive mapping step. Rejected for the dedicated MetaPhlAn image; revisit
  if validation (#88) shows no difference.
- **Keep HUMAnN in the base image (ADR-0001's position)** — ties the HUMAnN
  version to the base image rebuild cycle and to MetaPhlAn 4.2.2. Rejected;
  removing it from the base image is tracked separately (#89).

## Consequences

- Functional profiles can be produced without touching the published taxonomy;
  the cost is a second, HUMAnN-version-specific MetaPhlAn mapping per sample
  and an extra MetaPhlAn database (vJun23) in `store_dir`.
- HUMAnN-side databases live in new cache paths (`chocophlan/humann3.9/`, …);
  previously cached `chocophlan/full/` etc. are not reused and a one-time
  re-download is needed per `store_dir`.
- Adding a HUMAnN release (e.g. 4.0.0a1, #87) is a registry entry plus a
  container image; the publish layout and manifest keys already accommodate it.
- Stub tests do not exercise the containers, database URLs or tool command
  lines; the first real-sample run (#88) is the validation of the commands and
  of whether a HUMAnN 3.9 table set from a vJun23 profile is acceptable.
- Enabling HUMAnN by default is explicitly not decided here.

## References

- Issues #82 (design), #85 (this implementation), #83 (cache keys),
  #84 (`DATABASES`), #86, #87, #88, #89, #90.
- `conf/humann_bundles.config`, `modules/processes/humann.nf`,
  `modules/subworkflows/humann.nf`, `modules/processes/databases.nf`.
- HUMAnN 3.9 `humann/config.py`, `humann/search/prescreen.py`.
- bioBakery forum thread on HUMAnN 3.9 and MetaPhlAn 4.2.2 output.
