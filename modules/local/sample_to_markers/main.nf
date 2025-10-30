process SAMPLE_TO_MARKERS {
    tag "${meta.sample}"
    label 'process_low'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'

    input:
    tuple val(meta), path(metaphlan_sam)
    path metaphlan_db

    output:
    tuple val(meta), path("sample_to_markers"), emit: sample_to_markers
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir sample_to_markers

    sample2markers.py \\
        --input ${metaphlan_sam} \\
        --input_format bz2 \\
        --database ${params.store_dir}/metaphlan/\$(cat ${params.store_dir}/metaphlan/mpa_latest).pkl \\
        --nprocs ${task.cpus} \\
        --output_dir sample_to_markers

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sample2markers.py: \$(echo \$(sample2markers.py --version 2>&1) | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    mkdir sample_to_markers
    touch versions.yml
    """
}
