# Pre-full-scale-run verification checklist

The nf-test/`-stub-run` suite validates the **workflow wiring** (channels,
branching, output contract) but deliberately runs no containers and downloads no
databases. So a green test suite does **not** confirm that the pinned container
images exist, that the database URLs still resolve, or that the real tool
invocations behave as written.

Run through this checklist once before a full-scale run (and again whenever a
pinned tag or dated DB release is bumped). Everything here is cheap relative to a
full run and catches the failure modes stub tests cannot.

> Production uses Singularity (`singularity.enabled = true`). The checks below
> use `singularity`; Docker equivalents are noted where useful.

---

## 1. Container images resolve and carry the expected tools

Each image is pinned in code at the path shown. Confirm it pulls and the tool
runs.

- [ ] **Base image** — `docker://seandavi/curatedmetagenomics:metaphlan4.2.2`
      (`conf/base.config`)
      ```sh
      singularity exec docker://seandavi/curatedmetagenomics:metaphlan4.2.2 metaphlan --version
      singularity exec docker://seandavi/curatedmetagenomics:metaphlan4.2.2 kneaddata --version
      ```
- [ ] **Kraken2** — `docker://staphb/kraken2:2.1.3` (`modules/processes/kraken.nf`)
      ```sh
      singularity exec docker://staphb/kraken2:2.1.3 kraken2 --version
      ```
- [ ] **Bracken** — `docker://staphb/bracken:2.9` (`modules/processes/kraken.nf`)
      ```sh
      singularity exec docker://staphb/bracken:2.9 bracken -v
      ```
- [ ] **RGI** — `docker://quay.io/biocontainers/rgi:6.0.5--pyha8f3691_0`
      (`modules/processes/resistome.nf`). The `--pyha8f3691_0` build suffix is
      the most likely failure point; if the pull 404s, find the current build at
      <https://quay.io/repository/biocontainers/rgi?tab=tags> and update the tag.
      ```sh
      singularity exec docker://quay.io/biocontainers/rgi:6.0.5--pyha8f3691_0 rgi main --version
      ```

## 2. Database URLs resolve

- [ ] **MetaPhlAn index** `mpa_vJan25_CHOCOPhlAnSGB_202503` (`nextflow.config`
      `metaphlan_index`) — installs via `metaphlan --install`; confirm the index
      name is still offered for the MetaPhlAn version in the base image.
- [ ] **SGB→GTDB table** (`nextflow.config` `sgb2gtdb_url`) — must match
      `metaphlan_index`.
      ```sh
      curl -fsIL "https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan25_CHOCOPhlAnSGB_202503_SGB2GTDB.tsv" | head -1
      ```
- [ ] **Kraken2 PlusPF DB** (`nextflow.config` `kraken_db_url`) — dated release;
      confirm it still resolves (Langmead index releases roll over time, see
      <https://benlangmead.github.io/aws-indexes/k2>).
      ```sh
      curl -fsIL "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20260226.tar.gz" | head -1
      ```
- [ ] **Kraken DB has the Bracken distribution for the read length.** After the
      DB is fetched, the extracted `kraken_db/` must contain
      `database100mers.kmer_distrib` (matching `bracken_read_length = 100`). The
      PlusPF tarballs ship 50/75/100/150/200/250; if you change
      `bracken_read_length`, confirm the matching file exists.
- [ ] **CARD data** (`nextflow.config` `card_db_url`) — should resolve and the
      tarball should contain `card.json`.
      ```sh
      curl -fsIL "https://card.mcmaster.ca/latest/data" | head -1
      ```

## 3. Tool invocations not exercised by stub

These command sequences are written to documentation but only run for real
outside the stub. Validate them on **one real sample** (see §4) and confirm:

- [ ] **Kraken2** loads the DB into RAM within the provisioned memory tier
      (the `kraken2` process requests `32.GB * attempt`; the 16 GB DB plus
      overhead should fit — bump the process memory if pointing `kraken_db_url`
      at the full PlusPF DB).
- [ ] **Bracken** produces `bracken.species.txt` / `bracken.genus.txt` (no
      "kmer distribution not found" error → confirms §2 read-length item).
- [ ] **RGI bwt** completes for **single-end** input (`--read_one` only) and the
      `card_annotation` step produces a `card_database_v*.fasta` that the glob in
      `resistome.nf` picks up. Confirm `rgi_bwt.gene_mapping_data.txt` and
      `rgi_bwt.allele_mapping_data.txt` are non-empty.

## 4. End-to-end smoke test on one real sample

Run the **non-stub** pipeline on a single small accession with every branch
enabled, then confirm each published subdirectory is populated.

```sh
nextflow run . -profile <your_profile> \
  --sample_id SMOKE_TEST --run_ids <small_SRA_run_accession>
```

- [ ] `<sample>/manifest.json` exists and `read_accounting` shows non-zero
      `number_reads` for both `raw` and `decontaminated`.
- [ ] `<sample>/fastqc/` (decontaminated-read FastQC) is populated.
- [ ] `<sample>/full_data/metaphlan_lists/`, `metaphlan_markers/`,
      `strainphlan_markers/` populated.
- [ ] `<sample>/full_data/gtdb/gtdb_profile.tsv.gz` is non-empty.
- [ ] `<sample>/full_data/kraken/` has `kraken2.report.txt.gz` and the Bracken
      tables.
- [ ] `<sample>/full_data/resistome/` has the RGI mapping tables.
- [ ] `<sample>/rarefied_data/...` mirrors the full branch (minus `resistome/`,
      which is full-depth only).
- [ ] `<sample>/MARK_COMPLETE` exists (the run gated on every enabled branch).

## If something fails

- **Container tag rolled** → update the `container` line in `conf/base.config`
  (base image) or the relevant `modules/processes/*.nf` (in-body `container`).
- **DB dated release rolled** → update `kraken_db_url` in `nextflow.config`.
- **A step is genuinely broken on real data** → fix in the relevant module; the
  fix is isolated to that process. Re-run §4.

See the ADRs in [`docs/adr/`](adr/) for why each tool/DB was chosen
(0004 GTDB, 0006 Kraken2/Bracken, 0007 RGI/CARD, 0008 FastQC).
