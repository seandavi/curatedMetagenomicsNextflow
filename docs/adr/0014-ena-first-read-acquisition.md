# 0014. Acquire reads ENA-first with SRA fallback

- **Status:** Accepted
- **Date:** 2026-07-03
- **Deciders:** Sean Davis

## Context

`fasterq_dump` acquires reads by `curl`-ing the `.sra` object from the AWS Open
Data ODP mirror and running `fasterq-dump --split-files` on the local file.

A class of SRA runs is archived **without a QUALITY column**. For those,
`fasterq-dump` refuses to emit FASTQ and exits 3:

```
fasterq-dump err: the input data is missing the QUALITY-column
```

Under the current retry policy (ADR-0011) a non-OOM failure retries once and is
then `ignore`d — so the sample is dropped, never reaches `MARK_COMPLETE`, and is
dead-lettered by the orchestrator. In the 2.0.7 batch **all 14 dead-lettered
jobs** failed this way (the download itself succeeds; the *dump* fails). The
failure is 100% deterministic — it is a property of the archived object, not the
transfer — so retrying, and equally prefetch-then-dump, cannot help.

The underlying reads exist. Probing the EBI ENA `filereport` API showed FASTQ is
published and downloadable for every one of these runs (and for a 50/50 sample
of already-completed runs). ENA generates FASTQ even for quality-less runs, so
fetching from ENA sidesteps the `fasterq-dump` wall entirely. ENA coverage is
very high but not universal — brand-new runs (ingestion lag), submitted-only
formats (BAM/CRAM without generated FASTQ), and controlled-access data may have
no ENA FASTQ.

## Decision

We will acquire reads **ENA-first, per run accession**: resolve `fastq_ftp` via
the ENA `filereport` API, download over HTTPS, and verify against the provided
`fastq_md5`. When ENA serves no FASTQ for a run, fall back to the existing SRA
path (`curl` the `.sra` + `fasterq-dump`). The downstream output contract
(`out.fastq.gz` + lightweight QC) is unchanged.

## Alternatives considered

- **prefetch → fasterq-dump / nf-core `sratools` modules** — rejected: the
  process already fetches the `.sra` then dumps locally, so it is already
  prefetch-equivalent. The failure is in the dump step; a fetch-method change
  fixes none of these runs.
- **Switch ingestion to nf-core/fetchngs** — rejected for now: fetchngs is a
  download-only pipeline that speaks none of our contracts (`md5(sorted SRRs)`
  sample id, weblog telemetry, `MARK_COMPLETE`, `manifest.json`, publishDir
  layout). Revisit as a *decoupled* ingestion stage only if download
  volume/failure at 100k+ scale justifies a second orchestration surface.
- **ENA-only (drop SRA)** — rejected: ENA coverage is not 100%; the SRA path
  remains the backstop for the ENA-missing tail.

## Consequences

- Removes the QUALITY-column dead-letter class; the common path is also simpler
  and faster (plain HTTPS download, no `.sra` intermediate or dump step).
- Adds EBI/ENA as a runtime network dependency (transatlantic from US clusters);
  per-run fallback to SRA covers ENA being unavailable.
- ENA FASTQ for quality-less runs carries **synthetic** quality scores, so the
  downstream quality-trim (kneaddata/trimmomatic) is effectively a no-op for
  those samples. This is a curation consideration, not a pipeline bug — decide
  separately whether such samples belong in the dataset.
- ENA and `fasterq-dump` FASTQ are not byte-identical (read headers, paired-end
  splitting); immaterial for taxonomic profiling but a reproducibility note if a
  cohort mixes sources.
- No new container: `curl`, `md5sum`, `cut`, `gunzip` are all in the base image.
- The local-FASTQ input path (`local_fastqc`, cloud-bucket inputs) is out of
  scope here and still needs its own design pass (identity keys, metadata,
  user-supplied md5).

## References

- 2.0.7 dead-letter incident: 14/14 DLQ jobs = `fasterq-dump` exit 3,
  `missing the QUALITY-column`; ENA `filereport` confirmed FASTQ live for all.
- ENA filereport API: `https://www.ebi.ac.uk/ena/portal/api/filereport`
- ADR-0011 (retry/ignore policy), `modules/processes/preprocessing.nf`
