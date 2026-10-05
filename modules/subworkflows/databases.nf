/*
 * DATABASES subworkflow
 *
 * The single place reference-database processes are invoked. Used both by the
 * normal pipeline (which consumes the emitted channels) and by
 * `--databases_only`, which runs nothing else so the storeDir caches can be
 * populated ahead of a batch.
 *
 * Databases are gated by the same feature flags as the profiling steps:
 *   - MetaPhlAn + KneadData (human, mouse): always
 *   - Kraken2:                              unless skip_kraken
 *   - CARD + KMA index:                     unless skip_resistome
 *   - ChocoPhlAn, UniRef, utility mapping:  only when !skip_humann
 *
 * Channels for gated-off databases are empty; the main workflow only consumes
 * them under the same flag, so they are never read.
 */

include {
    install_metaphlan_db
    chocophlan_db
    utility_mapping_db
    uniref_db
    kneaddata_human_database
    kneaddata_mouse_database
    kraken_db
    card_db
    card_kma_db
} from '../processes/databases'

workflow DATABASES {

    main:
    install_metaphlan_db()

    // Both kneaddata database setup processes are invoked because downstream
    // wiring expects both channels to exist. The selected reference remains
    // controlled by organism_database.
    kneaddata_human_database()
    kneaddata_mouse_database()

    kraken_ch = Channel.empty()
    if (!params.skip_kraken) {
        kraken_db()
        kraken_ch = kraken_db.out.kraken_db
    }

    card_ch = Channel.empty()
    card_kma_ch = Channel.empty()
    if (!params.skip_resistome) {
        card_db()
        card_kma_db(card_db.out.card_db)
        card_ch = card_db.out.card_db
        card_kma_ch = card_kma_db.out.card_kma_db
    }

    chocophlan_ch = Channel.empty()
    uniref_ch = Channel.empty()
    utility_mapping_ch = Channel.empty()
    if (!params.skip_humann) {
        chocophlan_db()
        uniref_db()
        utility_mapping_db()
        chocophlan_ch = chocophlan_db.out.chocophlan_db
        uniref_ch = uniref_db.out.uniref_db
        utility_mapping_ch = utility_mapping_db.out.utility_mapping_db
    }

    emit:
    metaphlan_db       = install_metaphlan_db.out.metaphlan_db
    kd_genome          = kneaddata_human_database.out.kd_genome
    kd_mouse           = kneaddata_mouse_database.out.kd_mouse
    kraken_db          = kraken_ch
    card_db            = card_ch
    card_kma_db        = card_kma_ch
    chocophlan_db      = chocophlan_ch
    uniref_db          = uniref_ch
    utility_mapping_db = utility_mapping_ch
}
