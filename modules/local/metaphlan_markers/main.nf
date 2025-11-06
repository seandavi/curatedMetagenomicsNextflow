process METAPHLAN_MARKERS {
    tag "${meta.sample}"
    label 'process_low'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'

    input:
    tuple val(meta), path(metaphlan_bt2)
    path metaphlan_db

    output:
    tuple val(meta), path("marker_abundance.tsv.gz"), emit: marker_abundance
    tuple val(meta), path("marker_presence.tsv.gz"), emit: marker_presence
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    metaphlan --input_type bowtie2out \\
        --index ${params.metaphlan_index} \\
        --bowtie2db metaphlan \\
        -t marker_pres_table \\
        -o marker_presence.tsv \\
        <( gunzip -c ${metaphlan_bt2} )
    
    metaphlan --input_type bowtie2out \\
        --index ${params.metaphlan_index} \\
        --bowtie2db metaphlan \\
        -t marker_ab_table \\
        -o marker_abundance.tsv \\
        <( gunzip -c ${metaphlan_bt2} )
    
    gzip *.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(echo \$(metaphlan --version 2>&1) | awk '{print \$3}')
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch marker_abundance.tsv.gz
    touch marker_presence.tsv.gz
    touch versions.yml
    """
}
