/*
 * Reference database setup processes
 *
 * These tasks populate storeDir-backed assets so expensive downloads and
 * indexing work can be reused across runs.
 */

process install_metaphlan_db {
    label 'db_setup'
    label 'download_retry'

    cpus 4
    memory "8g"

    storeDir "${params.store_dir}"

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
    metaphlan --install --index ${params.metaphlan_index} --db_dir metaphlan

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS

    """
}

process chocophlan_db {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "1g"

    storeDir "${params.store_dir}"

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
    """
    echo ${PWD}
    humann_databases --update-config no --download chocophlan ${params.chocophlan} .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}

process utility_mapping_db {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "1g"

    storeDir "${params.store_dir}"

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
    """
    echo ${PWD}
    humann_databases --update-config no --download utility_mapping full .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}

process uniref_db {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "1g"

    storeDir "${params.store_dir}"

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
    """
    echo ${PWD}
    humann_databases --update-config no --download uniref ${params.uniref} .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}

process sgb_to_gtdb_db {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "1g"

    storeDir "${params.store_dir}"

    output:
    path "sgb_to_gtdb", emit: sgb2gtdb_db, type: 'dir'
    path ".command*"

    stub:
    """
    mkdir -p sgb_to_gtdb
    touch sgb_to_gtdb/${file(params.sgb2gtdb_url).name}
    touch .command.run
    """

    script:
    """
    echo ${PWD}
    mkdir -p sgb_to_gtdb
    # Download with the Python stdlib so we do not depend on curl/wget being
    # present in the container (the MetaPhlAn image always ships Python).
    python -c "import urllib.request; urllib.request.urlretrieve('${params.sgb2gtdb_url}', 'sgb_to_gtdb/${file(params.sgb2gtdb_url).name}')"
    """
}

process kraken_db {
    label 'db_setup'
    label 'download_retry'

    cpus 1
    memory "4g"

    storeDir "${params.store_dir}"

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

    storeDir "${params.store_dir}"

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

    storeDir "${params.store_dir}"

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
