/*
 * HUMANN subworkflow (ADR-0016, amended by ADR-0018)
 *
 * Full-depth branch only. Produces the MetaPhlAn profile HUMAnN stratifies by
 * and feeds it to HUMAnN. The bundle (HUMAnN container, MetaPhlAn profile, DB
 * names) is selected by params.humann_bundle; the databases come from
 * DATABASES.
 *
 *   - Bundle profile != params.metaphlan_profile: the bundle's own MetaPhlAn
 *     pass runs (metaphlan_for_humann) against its own index.
 *   - Bundle profile == params.metaphlan_profile: no second MetaPhlAn. The main
 *     full-branch rel_ab_w_read_stats profile is HUMAnN's input, and a copy is
 *     published under humann/<bundle>/metaphlan/ (humann_reuse_metaphlan).
 *
 * Nothing here changes the main MetaPhlAn taxonomy outputs.
 *
 * Samples are paired by sample id (not by channel position), so a failed or
 * reordered MetaPhlAn task can never hand HUMAnN another sample's reads.
 */

include { humann_reuses_main_metaphlan } from '../lib/metaphlan_profiles'
include { metaphlan_for_humann; humann_reuse_metaphlan; humann } from '../processes/humann'

workflow HUMANN {

    take:
    meta_ch             // val meta (one per sample)
    fastq_ch            // host-decontaminated full-depth reads, aligned with meta_ch
    main_profile_ch     // tuple(meta, main full-branch rel_ab_w_read_stats.tsv); used only on reuse
    metaphlan_humann_db // bundle MetaPhlAn index (DATABASES.out.metaphlan_humann_db); empty on reuse
    chocophlan_db
    uniref_db
    utility_mapping_db

    main:
    if (humann_reuses_main_metaphlan()) {
        humann_reuse_metaphlan(main_profile_ch)
        profile_ch = humann_reuse_metaphlan.out.profile
        metaphlan_versions_ch = Channel.empty()
    } else {
        metaphlan_for_humann(
            meta_ch,
            fastq_ch,
            metaphlan_humann_db.collect())
        profile_ch = metaphlan_for_humann.out.profile
        metaphlan_versions_ch = metaphlan_for_humann.out.meta
            .merge(metaphlan_for_humann.out.versions)
            .map { m, ver -> tuple(m.sample, ver) }
    }

    reads_by_sample = meta_ch
        .merge(fastq_ch)
        .map { m, fq -> tuple(m.sample, m, fq) }

    humann_input = reads_by_sample
        .join(profile_ch.map { m, prof -> tuple(m.sample, prof) })
        .map { sample, m, fq, prof -> tuple(m, fq, prof) }

    humann(
        humann_input,
        chocophlan_db.collect(),
        uniref_db.collect(),
        utility_mapping_db.collect())

    emit:
    meta     = humann.out.meta
    versions = metaphlan_versions_ch
        .mix(humann.out.meta
            .merge(humann.out.versions)
            .map { m, ver -> tuple(m.sample, ver) })
}
