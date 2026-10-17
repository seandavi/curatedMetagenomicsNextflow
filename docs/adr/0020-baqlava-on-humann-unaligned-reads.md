# 0020. Run BAQLaVa on the reads the bundle's HUMAnN left unaligned

- **Status:** Proposed
- **Date:** 2026-10-17
- **Deciders:** Sean Davis

## Context

BAQLaVa 1.2.0 (viral profiling, #119) starts with a bacterial-depletion step:
`utility_scripts/run_bacterial_depletion.py` runs
`humann --input <reads> --bypass-translated-search [--taxonomic-profile P]`
and keeps HUMAnN's `<name>_bowtie2_unaligned.fa`. That call passes no
`--nucleotide-database` or `--utility-database`, so HUMAnN uses its configured
defaults. In our image (`docker/humann4a-baqlava`) those are the packaged
`chocophlan_DEMO` and `utility_DEMO`, and the config cannot be changed at run
time because the container is read-only under Apptainer. Depletion therefore
"succeeds" against a three-species demo database. Without `--taxonomic-profile`
HUMAnN runs MetaPhlAn, which has no index in the image, and BAQLaVa falls back
to the undepleted input. Both failures are silent: the step exits 0 either way.

The pipeline already runs that exact HUMAnN nucleotide search in the bundle's
`humann` task, against the real ChocoPhlAn with the bundle's MetaPhlAn profile.
That task keeps `out_humann_temp/out_bowtie2_unaligned.fa`. BAQLaVa accepts
such a file directly with `--humann-depleted-fasta`, which is how bioBakery
workflows run it.

## Decision

`baqlava` takes the `humann` task's `out_bowtie2_unaligned.fa` (a new,
unpublished `unaligned` output of `humann`) and runs
`baqlava --humann-depleted-fasta`. BAQLaVa's own depletion step does not run.
BAQLaVa therefore needs HUMAnN. `skip_baqlava=false` fails at start-up unless
`skip_humann=false` and the bundle sets `baqlava_container`, a new bundle field
(today only `humann4.0.0a1`).

When HUMAnN selects no species from the profile it skips the nucleotide search
and writes no unaligned file. `humann` then writes the input reads in that
file's format, the same thing BAQLaVa does when it falls back.

The step gates `MARK_COMPLETE`, like HUMAnN. A sample without its viral profile
is incomplete for a BAQLaVa registration.

## Alternatives considered

- **Kneaddata reads plus `--taxonomic-profile`**, as #119 first described.
  Depletion would run against the demo ChocoPhlAn (see Context) and repeat
  HUMAnN's most expensive nucleotide step for every sample.
- **A `humann` wrapper on PATH that injects `--nucleotide-database`
  `--utility-database`.** This would fix the database problem, but it hides
  behaviour in a shim and still pays for a second nucleotide search.
- **Not gating `MARK_COMPLETE` on `baqlava`.** That would let a sample count
  as complete without the output its registration exists to produce.
  Short-read samples are not affected either way: they fail MetaPhlAn, and so
  HUMAnN, before BAQLaVa runs (#115).

## Consequences

- No extra MetaPhlAn or ChocoPhlAn pass. The depleted reads are identical to
  BAQLaVa's own depletion with the same profile and databases (checked on the
  BAQLaVa demo: same profile values).
- `baqlava` runs after `humann`, not in parallel with it.
- `humann` output declarations now also require
  `out_humann_temp/out_bowtie2_unaligned.fa`, for every bundle.
- A deterministic BAQLaVa failure costs the sample its retries and then
  dead-letters it, as a HUMAnN failure does today.

## References

- curatedMetagenomicsNextflow #117, #118, #119; seandavi/nextflow_telemetry#303
- BAQLaVa `baqlava/baqlava.py` (`--humann-depleted-fasta`) and
  `utility_scripts/run_bacterial_depletion.py` at 9a458b1
