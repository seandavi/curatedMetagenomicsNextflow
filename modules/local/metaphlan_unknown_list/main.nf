process METAPHLAN_UNKNOWN_LIST {
    tag "${meta.sample}"
    label 'process_medium'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'

    input:
    tuple val(meta), path(metaphlan_bt2)
    path metaphlan_db

    output:
    tuple val(meta), path('metaphlan_unknown_list.tsv'), emit: metaphlan_unknown_list
    path 'metaphlan_unknown_list.tsv.gz', emit: metaphlan_unknown_list_gz
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    metaphlan \\
        --input_type bowtie2out \\
        --index ${params.metaphlan_index} \\
        --bowtie2db metaphlan \\
        --nproc ${task.cpus} \\
        --unclassified_estimation \\
        -o metaphlan_unknown_list.tsv \\
        <( gunzip -c ${metaphlan_bt2} )

    gzip -c metaphlan_unknown_list.tsv > metaphlan_unknown_list.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(echo \$(metaphlan --version 2>&1) | awk '{print \$3}')
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch metaphlan_unknown_list.tsv
    touch metaphlan_unknown_list.tsv.gz
    touch versions.yml
    """
}
