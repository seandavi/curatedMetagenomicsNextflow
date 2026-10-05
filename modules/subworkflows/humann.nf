/*
 * HUMANN subworkflow (ADR-0016)
 *
 * Full-depth branch only. Runs the bundle's version-matched MetaPhlAn pass and
 * feeds its profile to HUMAnN. The bundle (containers, MetaPhlAn index, DB
 * names) is selected by params.humann_bundle; the databases come from
 * DATABASES. Nothing here touches the main MetaPhlAn taxonomy outputs.
 *
 * Samples are paired by sample id (not by channel position), so a failed or
 * reordered MetaPhlAn task can never hand HUMAnN another sample's reads.
 */

include { metaphlan_for_humann; humann } from '../processes/humann'

workflow HUMANN {

    take:
    meta_ch             // val meta (one per sample)
    fastq_ch            // host-decontaminated full-depth reads, aligned with meta_ch
    metaphlan_humann_db // bundle MetaPhlAn index (DATABASES.out.metaphlan_humann_db)
    chocophlan_db
    uniref_db
    utility_mapping_db

    main:
    metaphlan_for_humann(
        meta_ch,
        fastq_ch,
        metaphlan_humann_db.collect())

    reads_by_sample = meta_ch
        .merge(fastq_ch)
        .map { m, fq -> tuple(m.sample, m, fq) }

    humann_input = reads_by_sample
        .join(metaphlan_for_humann.out.profile.map { m, prof -> tuple(m.sample, prof) })
        .map { sample, m, fq, prof -> tuple(m, fq, prof) }

    humann(
        humann_input,
        chocophlan_db.collect(),
        uniref_db.collect(),
        utility_mapping_db.collect())

    emit:
    meta     = humann.out.meta
    versions = metaphlan_for_humann.out.meta
        .merge(metaphlan_for_humann.out.versions)
        .map { m, ver -> tuple(m.sample, ver) }
        .mix(humann.out.meta
            .merge(humann.out.versions)
            .map { m, ver -> tuple(m.sample, ver) })
}
