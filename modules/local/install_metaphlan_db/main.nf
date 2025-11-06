process INSTALL_METAPHLAN_DB {
    label 'process_low'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'
    storeDir "${params.store_dir}"

    output:
    path 'metaphlan', emit: metaphlan_db, type: 'dir'
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo ${PWD}
    metaphlan --install --index latest --bowtie2db metaphlan

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(echo \$(metaphlan --version 2>&1) | awk '{print \$3}')
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p metaphlan
    touch metaphlan/db.fake
    touch versions.yml
    """
}
