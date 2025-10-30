process UNIREF_DB {
    label 'process_single'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'
    storeDir "${params.store_dir}"

    output:
    path "uniref", emit: uniref_db, type: 'dir'
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo ${PWD}
    humann_databases --update-config no --download uniref ${params.uniref} .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1) | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p uniref
    touch uniref/db.fake
    touch versions.yml
    """
}
