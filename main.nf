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
} from './modules/processes/databases'
include {
    kneaddata
    metaphlan_unknown_viruses_lists
    metaphlan_unknown_list
    metaphlan_markers
    sample_to_markers
} from './modules/processes/profiling'
include { humann } from './modules/processes/humann'
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
        if (params.run_ids == null) or (params.sample_id == null) {
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

    /*
     * Core preprocessing and profiling section
     */
    if (params.local_input) {
        kneaddata(
            local_fastqc.out.meta,
            local_fastqc.out.fastq,
            kneaddata_human_database.out.kd_genome.collect(),
            kneaddata_mouse_database.out.kd_mouse.collect())
    } else {
        kneaddata(
            fasterq_dump.out.meta,
            fasterq_dump.out.fastq,
            kneaddata_human_database.out.kd_genome.collect(),
            kneaddata_mouse_database.out.kd_mouse.collect())
    }

    metaphlan_unknown_viruses_lists(
        kneaddata.out.meta,
        kneaddata.out.fastq,
        install_metaphlan_db.out.metaphlan_db.collect())

    metaphlan_markers(
        metaphlan_unknown_viruses_lists.out.meta,
        metaphlan_unknown_viruses_lists.out.metaphlan_bt2,
        install_metaphlan_db.out.metaphlan_db.collect())

    sample_to_markers(
        metaphlan_unknown_viruses_lists.out.meta,
        metaphlan_unknown_viruses_lists.out.metaphlan_sam,
        install_metaphlan_db.out.metaphlan_db.collect())

    /*
     * Optional HUMAnN branch
     */
    if (!params.skip_humann) {
        humann(
           kneaddata.out.meta,
           kneaddata.out.fastq,
           metaphlan_markers.out.marker_rel_ab_w_read_stats,
           chocophlan_db.out.chocophlan_db,
           uniref_db.out.uniref_db,
           utility_mapping_db.out.utility_mapping_db)
    }

    /*
     * Build a per-sample completion channel that only emits once both
     * MetaPhlAn-derived branches have completed for the same sample.
     */
    finished_ch = metaphlan_markers.out.meta
        .map { meta -> tuple(meta.sample, meta) }
        .join(sample_to_markers.out.meta.map { meta -> tuple(meta.sample, meta) })
        .map { sample_id, meta1, meta2 -> meta1 }

    if (!params.skip_humann) {
        // Also gate completion on HUMAnN when that branch is enabled.
        finished_ch = finished_ch
            .map { meta -> tuple(meta.sample, meta) }
            .join(humann.out.meta.map { meta -> tuple(meta.sample, meta) })
            .map { sample_id, meta1, meta2 -> meta1 }
    }

    MARK_COMPLETE(finished_ch)
}
