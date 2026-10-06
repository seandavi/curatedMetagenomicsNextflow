# 0018. Named MetaPhlAn profiles; HUMAnN bundles reference them and reuse the main pass when they match

- **Status:** Accepted (amends [0016](0016-humann-bundles.md): the bundle field list)
- **Date:** 2026-10-05
- **Deciders:** Sean Davis

## Context

[ADR-0016](0016-humann-bundles.md) made each HUMAnN bundle carry its own
MetaPhlAn pins (`metaphlan_container`, `metaphlan_version`, `metaphlan_index`,
`metaphlan_db_option`), while the main taxonomy pass took its index from
`params.metaphlan_index` and its container from the base image. MetaPhlAn
version facts therefore lived in three places, and the Alpine pilot (#88)
tripped over the version-specific CLI twice (#95):

- MetaPhlAn 4.1.x calls the database directory option `--bowtie2db` and writes
  the alignment map with `--bowtie2out`; 4.2.x uses `--db_dir` and `--mapout`
  and rejects `--bowtie2out` (verified in the containers).
- 4.1.2's `--version` prints an extra line that breaks HUMAnN 4.0.0a1's version
  parse, so the 4.1.x pins are 4.1.1.

The pipeline also installed the MetaPhlAn index with two near-duplicate
processes (`install_metaphlan_db`, `metaphlan_db_humann`), and always ran a
second MetaPhlAn pass for HUMAnN even when a HUMAnN release can consume the main
pass's profile (HUMAnN 4 final with vJan25).

## Decision

A **MetaPhlAn profile** is a named, pinned unit in one registry,
`conf/metaphlan_profiles.config` (`params.metaphlan_profiles`): `container`,
`version`, `index`, `db_option`, `map_option` and `map_input_type` (the
`--input_type` value that reads a map file back: `bowtie2out` in 4.1.x,
`mapout` in 4.2). Registered profiles: `mpa4.2.2_vJan25` (the base image, index
`mpa_vJan25_CHOCOPhlAnSGB_202503`), `mpa4.1.1_vJun23` and `mpa4.1.1_vOct22`
(`quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0`).

- **`params.metaphlan_profile`** (default `mpa4.2.2_vJan25`) replaces
  `params.metaphlan_index` with no alias. It selects the container, index and CLI
  options of every main-pass MetaPhlAn step and of the index install. An unknown
  name fails at workflow start, listing the valid names, whether or not HUMAnN
  runs (the main pass always needs a profile).
- **Bundles name a profile** (`metaphlan_profile`) instead of carrying
  `metaphlan_container` / `metaphlan_version` / `metaphlan_index` /
  `metaphlan_db_option`. The other bundle fields (`humann_container`, DB names,
  `humann_utility_db_option`) are unchanged. A bundle naming an unregistered
  profile fails at start-up the same way.
- **Reuse.** When the bundle's profile equals `params.metaphlan_profile`,
  `HUMANN` takes the main full-branch `rel_ab_w_read_stats` profile; neither
  `metaphlan_for_humann` nor a second index install runs. A copy is published
  under `humann/<bundle>/metaphlan/` (process `humann_reuse_metaphlan`), so the
  provenance layout is the same either way. When the profiles differ, behaviour
  is as in ADR-0016.
- **One install process.** `install_metaphlan_db` takes a profile name and runs
  in that profile's container with its `db_option`. `DATABASES` calls it for the
  main profile and, only when a bundle needs a different one, again (aliased as
  `metaphlan_db_humann`).
- **The main pass's container does not change.** `mpa4.2.2_vJan25` names the
  base image's container string, equal to `process.container` in
  `conf/base.config`, so the published taxonomy is unchanged by this decision.
  Moving the main pass to a biocontainer could change the bowtie2/MetaPhlAn
  builds and so the taxonomy, and needs a file-for-file comparison; it is a
  separate decision.
- **Cache keys stay index-based** (`<store_dir>/metaphlan/<index>/`), so existing
  stores are reused and two profiles with the same index share one database.
- **Manifest.** `parameters.metaphlan_index` becomes `parameters.metaphlan_profile`,
  and `parameters.humann_metaphlan_index` becomes `humann_metaphlan_profile`
  (equal to the main profile on reuse).

## Alternatives considered

- **Keep per-bundle MetaPhlAn fields plus `params.metaphlan_index`** — version
  facts stay in three places and each new CLI difference needs a new field.
  Rejected.
- **Keep an alias for `metaphlan_index`** — two ways to say the same thing, and
  the index alone cannot say which container or CLI to use. Rejected.
- **Always run the HUMAnN-side MetaPhlAn pass** — a redundant mapping of every
  sample against the same index once a HUMAnN release accepts the main profile.
  Rejected.
- **Move the main pass to a biocontainer in the same change** — out of scope; see
  above.

## Consequences

- A new MetaPhlAn version (or CLI spelling) is one registry entry; a new HUMAnN
  release that accepts the main profile is a bundle whose `metaphlan_profile` is
  the main one, with no further code.
- The registry's `mpa4.2.2_vJan25` container and `process.container` are two
  spellings of one pin and must be kept equal by hand.
- Downstream readers of the manifest must read `metaphlan_profile` where they
  read `metaphlan_index`.
- No registry bundle uses the reuse path yet (it needs a HUMAnN release that
  accepts vJan25); it is exercised by a stub test with a test-only bundle.

## References

- #99 (this change), #82 (plan), #95 (CLI difference), #87 (HUMAnN 4 alpha
  bundle), #88 (pilot).
- `conf/metaphlan_profiles.config`, `conf/humann_bundles.config`,
  `modules/lib/metaphlan_profiles.nf`, `modules/processes/databases.nf`,
  `modules/subworkflows/humann.nf`.
- ADR-0016 (HUMAnN bundles), ADR-0017 (self-built images).
