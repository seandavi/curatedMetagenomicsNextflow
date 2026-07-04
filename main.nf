#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Curated Metagenomic Data pipeline
 *
 * Refactor invariant:
 * - Published directory names and file names must remain unchanged.
 * - Process output file names must remain unchanged.
 * - Workflow branching semantics must remain unchanged.
 *
 * This entrypoint intentionally focuses on orchestration. Process definitions
 * live under modules/processes/ so users can navigate the workflow graph
 * without scrolling through every task implementation first.
 */

include { fasterq_dump; local_fastqc } from './modules/processes/preprocessing'
include {
    install_metaphlan_db
    chocophlan_db
    utility_mapping_db
    uniref_db
    kneaddata_human_database
    kneaddata_mouse_database
    sgb_to_gtdb_db
    kraken_db
    card_db
    card_kma_db
} from './modules/processes/databases'
include {
    kneaddata
    metaphlan_unknown_viruses_lists as metaphlan_unknown_viruses_lists_full
    metaphlan_unknown_viruses_lists as metaphlan_unknown_viruses_lists_rarefied
    metaphlan_unknown_list
    metaphlan_markers as metaphlan_markers_full
    metaphlan_markers as metaphlan_markers_rarefied
    sample_to_markers as sample_to_markers_full
    sample_to_markers as sample_to_markers_rarefied
} from './modules/processes/profiling'
include { rarefy_fastq } from './modules/processes/rarefaction'
include {
    metaphlan_to_gtdb as metaphlan_to_gtdb_full
    metaphlan_to_gtdb as metaphlan_to_gtdb_rarefied
} from './modules/processes/gtdb'
include {
    kraken2 as kraken2_full
    kraken2 as kraken2_rarefied
    bracken as bracken_full
    bracken as bracken_rarefied
} from './modules/processes/kraken'
include {
    resistome_kma as resistome_kma_full
    resistome_kma as resistome_kma_rarefied
} from './modules/processes/resistome'
include { fastqc } from './modules/processes/qc'
include { humann } from './modules/processes/humann'
include { sample_manifest } from './modules/processes/manifest'
include { MARK_COMPLETE } from './modules/processes/finalize'

/*
 * Metadata helpers normalize both TSV-driven and single-sample invocations
 * into a common metadata structure used throughout the pipeline:
 *
 *   [sample: <sample_id>, accessions|fpaths: [...], meta: <source row or null>]
 */
def generate_row_tuple(row) {
    accessions = row.NCBI_accession.split(';');
    sample_id = row.sample_id;
    return [sample:sample_id, accessions: accessions, meta: row]
}

def generate_row_tuple_local(row) {
    fpaths = row.file_paths.split(';');
    sample_id = row.sample_id;
    return [sample:sample_id, fpaths: fpaths, meta: row]
}

def generate_sample_metadata_single_sample(sample_id, run_ids) {
    accessions = run_ids.split(';')
    return [sample: sample_id, accessions: accessions, meta: null]
}

workflow {

    samples = null

    /*
     * Input contract:
     * - Either provide a metadata TSV
     * - Or provide both run_ids and sample_id for a single sample
     */
    if (params.metadata_tsv == null) {
        if (params.run_ids == null || params.sample_id == null) {
            error "Either metadata_tsv or run_ids/sample_id must be provided"
        } else {
            samples = generate_sample_metadata_single_sample(
                params.sample_id,
                params.run_ids
            )
            fasterq_dump(samples)
        }
    } else {
        if (params.local_input) {
            samples = Channel
                .fromPath(params.metadata_tsv)
                .splitCsv(header: true, quote: '"', sep:'\t')
                .map { row -> generate_row_tuple_local(row) }

            local_fastqc(samples)
        } else {
            samples = Channel
                .fromPath(params.metadata_tsv)
                .splitCsv(header: true, quote: '"', sep:'\t')
                .map { row -> generate_row_tuple(row) }

            fasterq_dump(samples)
        }
    }

    /*
     * Database setup section
     *
     * Current behavior note:
     * both kneaddata database setup processes are invoked because downstream
     * wiring currently expects both channels to exist. This is a refactor
     * target, but the behavior is preserved here intentionally.
     */
    install_metaphlan_db()
    uniref_db()
    chocophlan_db()
    utility_mapping_db()
    kneaddata_human_database()
    kneaddata_mouse_database()
    if (!params.skip_gtdb) {
        sgb_to_gtdb_db()
    }
    if (!params.skip_kraken) {
        kraken_db()
    }
    if (!params.skip_resistome) {
        card_db()
        card_kma_db(card_db.out.card_db)
    }

    /*
     * Core preprocessing and profiling section
     *
     * When skip_rarefied is false (the default), both a full and a rarefied
     * profiling branch run in parallel. The branch is identified by the
     * meta.branch key ('full_data' or 'rarefied_data'), which is used by the
     * profiling processes to write outputs into distinct subdirectories.
     *
     * When skip_rarefied is true, only the full branch runs and meta.branch is
     * not set, so publishDir paths remain identical to previous pipeline
     * versions (backward-compatible output layout).
     */
    if (params.local_input) {
        raw_fastq_ch = local_fastqc.out.fastq
        raw_versions_ch = local_fastqc.out.versions
        raw_meta_ch = local_fastqc.out.meta
        kneaddata(
            local_fastqc.out.meta,
            local_fastqc.out.fastq,
            kneaddata_human_database.out.kd_genome.collect(),
            kneaddata_mouse_database.out.kd_mouse.collect())
    } else {
        raw_fastq_ch = fasterq_dump.out.fastq
        raw_versions_ch = fasterq_dump.out.versions
        raw_meta_ch = fasterq_dump.out.meta
        kneaddata(
            fasterq_dump.out.meta,
            fasterq_dump.out.fastq,
            kneaddata_human_database.out.kd_genome.collect(),
            kneaddata_mouse_database.out.kd_mouse.collect())
    }

    /*
     * Determine the meta channel for the full branch.
     *
     * When the rarefied branch is also active we tag meta with branch='full_data'
     * so that the profiling processes publish outputs into a branch-specific
     * subdirectory.  When skip_rarefied is true the meta is left unmodified so
     * that output paths remain identical to the pre-dual-branch layout.
     */
    if (!params.skip_rarefied) {
        full_meta_ch = kneaddata.out.meta.map { meta -> meta + [branch: 'full_data'] }
    } else {
        full_meta_ch = kneaddata.out.meta
    }

    metaphlan_unknown_viruses_lists_full(
        full_meta_ch,
        kneaddata.out.fastq,
        install_metaphlan_db.out.metaphlan_db.collect())

    metaphlan_markers_full(
        metaphlan_unknown_viruses_lists_full.out.meta,
        metaphlan_unknown_viruses_lists_full.out.metaphlan_bt2,
        install_metaphlan_db.out.metaphlan_db.collect())

    sample_to_markers_full(
        metaphlan_unknown_viruses_lists_full.out.meta,
        metaphlan_unknown_viruses_lists_full.out.metaphlan_sam,
        install_metaphlan_db.out.metaphlan_db.collect())

    /*
     * GTDB-taxonomy conversion of the full-depth MetaPhlAn profile.
     */
    if (!params.skip_gtdb) {
        metaphlan_to_gtdb_full(
            metaphlan_unknown_viruses_lists_full.out.meta,
            metaphlan_unknown_viruses_lists_full.out.metaphlan_unknown_list,
            sgb_to_gtdb_db.out.sgb2gtdb_db.collect())
    }

    /*
     * Complementary read-based profiling (Kraken2 + Bracken) on the full-depth
     * host-decontaminated reads.
     */
    if (!params.skip_kraken) {
        kraken2_full(
            full_meta_ch,
            kneaddata.out.fastq,
            kraken_db.out.kraken_db.collect())

        bracken_full(
            kraken2_full.out.meta,
            kraken2_full.out.report,
            kraken_db.out.kraken_db.collect())
    }

    /*
     * Resistome profiling (KMA against CARD) on the full-depth
     * host-decontaminated reads. KMA is cheap enough to run on both branches
     * (see the rarefied section below). See ADR-0012.
     */
    if (!params.skip_resistome) {
        resistome_kma_full(
            full_meta_ch,
            kneaddata.out.fastq,
            card_kma_db.out.card_kma_db.collect())
    }

    /*
     * Per-sample QC: FastQC on the host-decontaminated reads (raw-read FastQC
     * already runs in fasterq_dump/local_fastqc).
     */
    if (!params.skip_fastqc) {
        fastqc(
            kneaddata.out.meta,
            kneaddata.out.fastq)
    }

    /*
     * Per-sample manifest
     *
     * Join (by sample) the raw reads + their versions, the host-decontaminated
     * reads + versions, and the full-branch MetaPhlAn versions, then compile a
     * single manifest.json per sample. Each output channel is merged from the
     * same originating process so positions stay task-aligned before the join.
     */
    raw_for_manifest = raw_meta_ch
        .merge(raw_fastq_ch)
        .merge(raw_versions_ch)
        .map { m, fq, ver -> tuple(m.sample, m, fq, ver) }

    kneaddata_for_manifest = kneaddata.out.meta
        .merge(kneaddata.out.fastq)
        .merge(kneaddata.out.versions)
        .map { m, fq, ver -> tuple(m.sample, fq, ver) }

    metaphlan_versions_for_manifest = metaphlan_unknown_viruses_lists_full.out.meta
        .merge(metaphlan_unknown_viruses_lists_full.out.versions)
        .map { m, ver -> tuple(m.sample, ver) }

    manifest_input = raw_for_manifest
        .join(kneaddata_for_manifest)
        .join(metaphlan_versions_for_manifest)
        .map { sample, meta, raw_fq, raw_ver, kd_fq, kd_ver, mp_ver ->
            tuple(meta, raw_fq, raw_ver, kd_fq, kd_ver, mp_ver) }

    sample_manifest(manifest_input)

    /*
     * Optional rarefied branch
     *
     * Downsample reads with seqtk then run the same profiling steps under the
     * 'rarefied_data' subdirectory so the outputs can be compared side-by-side
     * with the full-depth results.
     */
    if (!params.skip_rarefied) {
        rarefied_meta_ch = kneaddata.out.meta.map { meta -> meta + [branch: 'rarefied_data'] }

        rarefy_fastq(
            rarefied_meta_ch,
            kneaddata.out.fastq)

        metaphlan_unknown_viruses_lists_rarefied(
            rarefy_fastq.out.meta,
            rarefy_fastq.out.fastq,
            install_metaphlan_db.out.metaphlan_db.collect())

        metaphlan_markers_rarefied(
            metaphlan_unknown_viruses_lists_rarefied.out.meta,
            metaphlan_unknown_viruses_lists_rarefied.out.metaphlan_bt2,
            install_metaphlan_db.out.metaphlan_db.collect())

        sample_to_markers_rarefied(
            metaphlan_unknown_viruses_lists_rarefied.out.meta,
            metaphlan_unknown_viruses_lists_rarefied.out.metaphlan_sam,
            install_metaphlan_db.out.metaphlan_db.collect())

        if (!params.skip_gtdb) {
            metaphlan_to_gtdb_rarefied(
                metaphlan_unknown_viruses_lists_rarefied.out.meta,
                metaphlan_unknown_viruses_lists_rarefied.out.metaphlan_unknown_list,
                sgb_to_gtdb_db.out.sgb2gtdb_db.collect())
        }

        if (!params.skip_kraken) {
            kraken2_rarefied(
                rarefy_fastq.out.meta,
                rarefy_fastq.out.fastq,
                kraken_db.out.kraken_db.collect())

            bracken_rarefied(
                kraken2_rarefied.out.meta,
                kraken2_rarefied.out.report,
                kraken_db.out.kraken_db.collect())
        }

        if (!params.skip_resistome) {
            resistome_kma_rarefied(
                rarefy_fastq.out.meta,
                rarefy_fastq.out.fastq,
                card_kma_db.out.card_kma_db.collect())
        }
    }

    /*
     * Optional HUMAnN branch (uses full-depth data only)
     */
    if (!params.skip_humann) {
        humann(
           metaphlan_unknown_viruses_lists_full.out.meta,
           kneaddata.out.fastq,
           metaphlan_markers_full.out.marker_rel_ab_w_read_stats,
           chocophlan_db.out.chocophlan_db,
           uniref_db.out.uniref_db,
           utility_mapping_db.out.utility_mapping_db)
    }

    /*
     * Build a per-sample completion channel that fires only once all active
     * branches have finished for a given sample.
     *
     * Full branch: join metaphlan_markers_full and sample_to_markers_full.
     * Rarefied branch (when enabled): also join metaphlan_markers_rarefied and
     * sample_to_markers_rarefied.
     * HUMAnN branch (when enabled): additionally gate on humann output.
     * GTDB branch (when enabled): additionally gate on the full-branch GTDB
     * conversion.
     */
    finished_ch = metaphlan_markers_full.out.meta
        .map { meta -> tuple(meta.sample, meta) }
        .join(sample_to_markers_full.out.meta.map { meta -> tuple(meta.sample, meta) })
        .map { sample_id, meta1, meta2 -> meta1 }

    if (!params.skip_gtdb) {
        // Gate completion on the full-branch GTDB conversion when enabled.
        finished_ch = finished_ch
            .map { meta -> tuple(meta.sample, meta) }
            .join(metaphlan_to_gtdb_full.out.meta.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }
    }

    if (!params.skip_kraken) {
        // Gate completion on the full-branch Kraken2/Bracken profiling.
        finished_ch = finished_ch
            .map { meta -> tuple(meta.sample, meta) }
            .join(bracken_full.out.meta.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }
    }

    if (!params.skip_resistome) {
        // Gate completion on the full-branch resistome profiling.
        finished_ch = finished_ch
            .map { meta -> tuple(meta.sample, meta) }
            .join(resistome_kma_full.out.meta.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }
    }

    if (!params.skip_fastqc) {
        // Gate completion on the post-decontamination FastQC.
        finished_ch = finished_ch
            .map { meta -> tuple(meta.sample, meta) }
            .join(fastqc.out.meta.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }
    }

    if (!params.skip_rarefied) {
        finished_rarefied_ch = metaphlan_markers_rarefied.out.meta
            .map { meta -> tuple(meta.sample, meta) }
            .join(sample_to_markers_rarefied.out.meta.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }

        finished_ch = finished_ch
            .map { meta -> tuple(meta.sample, meta) }
            .join(finished_rarefied_ch.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }
    }

    if (!params.skip_humann) {
        // Also gate completion on HUMAnN when that branch is enabled.
        finished_ch = finished_ch
            .map { meta -> tuple(meta.sample, meta) }
            .join(humann.out.meta.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }
    }

    // Gate completion on the per-sample manifest so the published sample
    // directory is not marked complete before manifest.json exists.
    finished_ch = finished_ch
        .map { meta -> tuple(meta.sample, meta) }
        .join(sample_manifest.out.meta.map { meta -> tuple(meta.sample, meta) })
        .map { sample_id, meta1, meta2 -> meta1 }

    MARK_COMPLETE(finished_ch)
}
