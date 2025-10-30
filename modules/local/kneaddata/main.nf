process KNEADDATA {
    tag "${meta.sample}"
    label 'process_medium'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'

    input:
    tuple val(meta), path(fastq)
    path kd_genome
    path kd_mouse

    output:
    tuple val(meta), path("kneaddata_output/out.fastq"), emit: fastq
    path "kneaddata_output/kneaddata_fastq_linecounts.txt", emit: line_counts
    path "kneaddata_output/out_kneaddata.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kneaddata --unpaired ${fastq} \\
        --reference-db ${params.organism_database} \\
        --output kneaddata_output  \\
        --trimmomatic /installed/Trimmomatic-0.39 \\
        --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:30' \\
        --bypass-trf \\
        --bowtie2-options='--very-fast' \\
        -t 16 -p 8

    cd kneaddata_output
    cat out_kneaddata.fastq | sed 's/^+.RR.*/+/g' > out.fastq
    rm out_kneaddata.fastq
    wc -l * | grep fastq > kneaddata_fastq_linecounts.txt
    """

    stub:
    """
    mkdir -p kneaddata_output
    touch kneaddata_output/out.fastq
    touch kneaddata_output/kneaddata_fastq_linecounts.txt
    touch kneaddata_output/out_kneaddata.log
    """
}
