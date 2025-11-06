process METAPHLAN_UNKNOWN_VIRUSES_LISTS {
    tag "${meta.sample}"
    label 'process_medium'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'

    input:
    tuple val(meta), path(fastq)
    path metaphlan_db

    output:
    tuple val(meta), path('bowtie2.out.gz'), emit: metaphlan_bt2
    tuple val(meta), path('sam.bz2'), emit: metaphlan_sam
    path 'metaphlan_unknown_list.tsv', emit: metaphlan_unknown_list
    path 'metaphlan_unknown_list.tsv.gz', emit: metaphlan_unknown_list_gz
    path 'metaphlan_viruses_list.tsv', emit: metaphlan_viruses_list
    path 'metaphlan_viruses_list.tsv.gz', emit: metaphlan_viruses_list_gz
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    find .
    metaphlan --input_type fastq \\
        --index ${params.metaphlan_index} \\
        --bowtie2db metaphlan \\
        --samout sam.bz2 \\
        --bowtie2out bowtie2.out \\
        --nproc ${task.cpus} \\
        --unclassified_estimation \\
        --profile_vsc \\
        --vsc_breadth 0.75 \\
        --vsc_out metaphlan_viruses_list.tsv \\
        -o metaphlan_unknown_list.tsv \\
        ${fastq}

    gzip -c metaphlan_unknown_list.tsv > metaphlan_unknown_list.tsv.gz
    gzip -c metaphlan_viruses_list.tsv > metaphlan_viruses_list.tsv.gz
    gzip bowtie2.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(echo \$(metaphlan --version 2>&1) | awk '{print \$3}')
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | head -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch bowtie2.out.gz
    touch sam.bz2
    touch metaphlan_unknown_list.tsv
    touch metaphlan_unknown_list.tsv.gz
    touch metaphlan_viruses_list.tsv
    touch metaphlan_viruses_list.tsv.gz
    touch versions.yml
    """
}
