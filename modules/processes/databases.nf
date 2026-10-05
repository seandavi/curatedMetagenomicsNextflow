/*
 * Reference database setup processes
 *
 * These tasks populate storeDir-backed assets so expensive downloads and
 * indexing work can be reused across runs.
 *
 * Cache layout convention (issue #83)
 * -----------------------------------
 * storeDir only checks that a process's declared outputs exist, so a cache
 * path that does not encode the database version silently serves stale data
 * after a parameter change. Every parameter-dependent database therefore has
 * its own storeDir
 *
 *     ${params.store_dir}/<db_name>/<version_key>/
 *
 * holding the database directory (named <db_name>, unchanged from before, so
 * the staged input path is the same shape as ever) plus that task's
 * `.command*` and `versions.yml`. The version key is:
 *
 *     metaphlan        params.metaphlan_index
 *     metaphlan (HUMAnN-side, metaphlan_db_humann)
 *                      the selected bundle's metaphlan_index (same key space
 *                      as the main MetaPhlAn database)
 *     chocophlan, uniref, utility_mapping
 *                      params.humann_bundle (the database names are bundle
 *                      fields, and the same name can differ across HUMAnN
 *                      releases, so the bundle is the version)
 *     kraken_db        basename of params.kraken_db_url without archive extension
 *     card_db          basename of params.card_db_url without archive extension
 *     card_kma_db      the same key as card_db (the KMA index is built from it)
 *
 * Changing a parameter therefore creates a new directory beside the old one
 * instead of reusing it, and two versions can coexist (e.g. two MetaPhlAn
 * indexes). Downstream processes must use the staged input path, never a
 * hardcoded directory name.
 *
 * Not versioned: the KneadData human_genome / mouse_C57BL databases. They take
 * no pipeline parameter that selects a release (`kneaddata_database --download`
 * always fetches the current prebuilt bowtie2 index), so there is no version
 * value to key them by; their fixed paths stay as they were.
 */

// Version key for a database URL: the basename without archive extensions,
// e.g. .../k2_pluspf_16_GB_20260226.tar.gz -> k2_pluspf_16_GB_20260226.
def db_url_key(url) {
    def name = url.toString().split('\\?')[0].tokenize('/').last()
    return name.replaceAll(/(\.tar)?\.(gz|bz2|xz|zip)$|\.(tgz|tbz2|tar)$/, '')
}

process install_metaphlan_db {
    label 'db_setup'
    label 'download_retry'

    cpus 4
    memory "8g"

    storeDir "${params.store_dir}/metaphlan/${params.metaphlan_index}"

    output:
    path 'metaphlan', emit: metaphlan_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    stub:
    """
    mkdir -p metaphlan
    touch metaphlan/db.fake
    touch .command.run
    touch versions.yml
    """

    script:
    """
    echo ${PWD}
    metaphlan --install --index ${params.metaphlan_index} --db_dir ./metaphlan

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS

    """
}

/*
 * HUMAnN-side databases (ADR-0016)
 *
 * These four processes belong to the selected HUMAnN bundle
 * (params.humann_bundles[params.humann_bundle], conf/humann_bundles.config)
 * and run in the bundle's own containers, so the DBs are always built by the
 * tool version that will read them. They are invoked only from DATABASES,
 * only when !skip_humann.
 *
 * The MetaPhlAn index is keyed by index name (the same key as the main pass:
 * an identical index name is an identical database). ChocoPhlAn, UniRef and
 * utility mapping are keyed by *bundle*, not by database name, because the
 * same names ("full", "uniref90_ec_filtered_diamond") can mean different
 * content under different HUMAnN releases.
 */

process metaphlan_db_humann {
    label 'db_setup'
    label 'download_retry'

    container { params.humann_bundles[params.humann_bundle].metaphlan_container }

    cpus 4
    memory { 8.GB * task.attempt }

    storeDir "${params.store_dir}/metaphlan/${params.humann_bundles[params.humann_bundle].metaphlan_index}"

    output:
    path 'metaphlan', emit: metaphlan_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    stub:
    """
    mkdir -p metaphlan
    touch metaphlan/db.fake
    touch .command.run
    touch versions.yml
    """

    script:
    def bundle = params.humann_bundles[params.humann_bundle]
    """
    metaphlan --install --index ${bundle.metaphlan_index} --db_dir ./metaphlan

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan_humann: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2_humann: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}

process chocophlan_db {
    label 'db_setup'
    label 'download_retry'

    container { params.humann_bundles[params.humann_bundle].humann_container }

    cpus 1
    memory { 1.GB * task.attempt }

    storeDir "${params.store_dir}/chocophlan/${params.humann_bundle}"

    output:
    path "chocophlan", emit: chocophlan_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    stub:
    """
    mkdir -p chocophlan
    touch chocophlan/db.fake
    touch .command.run
    touch versions.yml
    """

    script:
    def bundle = params.humann_bundles[params.humann_bundle]
    """
    humann_databases --update-config no --download chocophlan ${bundle.chocophlan} .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}

process utility_mapping_db {
    label 'db_setup'
    label 'download_retry'

    container { params.humann_bundles[params.humann_bundle].humann_container }

    cpus 1
    memory { 1.GB * task.attempt }

    storeDir "${params.store_dir}/utility_mapping/${params.humann_bundle}"

    output:
    path "utility_mapping", emit: utility_mapping_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    stub:
    """
    mkdir -p utility_mapping
    touch utility_mapping/db.fake
    touch .command.run
    touch versions.yml
    """

    script:
    def bundle = params.humann_bundles[params.humann_bundle]
    """
    humann_databases --update-config no --download utility_mapping ${bundle.utility_mapping} .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}

process uniref_db {
    label 'db_setup'
    label 'download_retry'

    container { params.humann_bundles[params.humann_bundle].humann_container }

    cpus 1
    memory { 1.GB * task.attempt }

    storeDir "${params.store_dir}/uniref/${params.humann_bundle}"

    output:
    path "uniref", emit: uniref_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    stub:
    """
    mkdir -p uniref
    touch uniref/db.fake
    touch .command.run
    touch versions.yml
    """

    script:
    def bundle = params.humann_bundles[params.humann_bundle]
    """
    humann_databases --update-config no --download uniref ${bundle.uniref} .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}

process kraken_db {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "4g"

    storeDir "${params.store_dir}/kraken_db/${db_url_key(params.kraken_db_url)}"

    output:
    path "kraken_db", emit: kraken_db, type: 'dir'
    path ".command*"

    stub:
    """
    mkdir -p kraken_db
    touch kraken_db/hash.k2d
    touch kraken_db/opts.k2d
    touch kraken_db/taxo.k2d
    touch kraken_db/database${params.bracken_read_length}mers.kmer_distrib
    touch .command.run
    """

    script:
    """
    echo ${PWD}
    mkdir -p kraken_db
    # Prebuilt Kraken2 index tarballs bundle the Bracken kmer distributions and
    # extract their .k2d files directly (no top-level directory), so unpack into
    # kraken_db/. Stream the download to disk to avoid buffering a large file.
    curl -fsSL "${params.kraken_db_url}" -o kraken_db.tar.gz
    tar -xzf kraken_db.tar.gz -C kraken_db
    rm -f kraken_db.tar.gz
    """
}

process card_db {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "2g"

    storeDir "${params.store_dir}/card_db/${db_url_key(params.card_db_url)}"

    output:
    path "card_db", emit: card_db, type: 'dir'
    path ".command*"

    stub:
    """
    mkdir -p card_db
    touch card_db/nucleotide_fasta_protein_homolog_model.fasta
    touch .command.run
    """

    script:
    """
    echo ${PWD}
    mkdir -p card_db
    # The CARD "broadstreet" release is a bzip2 tarball; we index the homolog-
    # model nucleotide FASTA with KMA (see card_kma_db). Extract with Python's
    # tarfile so we do not depend on a bzip2 binary being present in the image.
    curl -fsSL "${params.card_db_url}" -o card_data.tar.bz2
    python -c "import tarfile; tarfile.open('card_data.tar.bz2','r:bz2').extractall('card_db')"
    rm -f card_data.tar.bz2
    """
}

process card_kma_db {
    label 'db_setup'

    // KMA is not in the base image; index CARD in its own pinned biocontainer
    // (ADR-0001). This is a shared, storeDir-backed asset reused across runs.
    container 'docker://quay.io/biocontainers/kma:1.6.13--h118bc1c_0'

    cpus 2
    memory "8g"

    storeDir "${params.store_dir}/card_kma_db/${db_url_key(params.card_db_url)}"

    input:
    path card_db

    output:
    path "card_kma_db", emit: card_kma_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    stub:
    """
    mkdir -p card_kma_db
    touch card_kma_db/card_kma_db.comp.b
    touch card_kma_db/card_kma_db.length.b
    touch card_kma_db/card_kma_db.name
    touch card_kma_db/card_kma_db.seq.b
    touch .command.run
    touch versions.yml
    """

    script:
    """
    echo ${PWD}
    mkdir -p card_kma_db
    kma index \
        -i ${card_db}/nucleotide_fasta_protein_homolog_model.fasta \
        -o card_kma_db/card_kma_db

    cat <<-END_VERSIONS > versions.yml
    versions:
        kma: \$( kma -v 2>&1 | head -n1 | sed 's/^KMA-//' )
    END_VERSIONS
    """
}

process kneaddata_human_database {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "4g"

    storeDir "${params.store_dir}"

    output:
    path "human_genome", emit: kd_genome, type: "dir"
    path ".command*"

    stub:
    """
    mkdir -p human_genome
    touch human_genome/hg37dec_v0.1.1.bt2
    touch .command.run
    """

    script:
    """
    echo ${PWD}
    mkdir -p human_genome
    kneaddata_database --download human_genome bowtie2 human_genome
    """
}

process kneaddata_mouse_database {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "4g"

    storeDir "${params.store_dir}"

    output:
    path "mouse_C57BL", emit: kd_mouse, type: "dir"
    path ".command*"

    stub:
    """
    mkdir -p mouse_C57BL
    touch mouse_C57BL/mouse_C57BL_6NJ_Bowtie2_v0.1.bt2
    touch .command.run
    """

    script:
    """
    echo ${PWD}
    mkdir -p mouse_C57BL
    kneaddata_database --download mouse_C57BL bowtie2 mouse_C57BL
    """
}
