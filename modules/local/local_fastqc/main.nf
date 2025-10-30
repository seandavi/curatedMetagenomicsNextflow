process LOCAL_FASTQC {
    tag "${meta.sample}"
    label 'process_medium'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'

    input:
    val meta

    output:
    tuple val(meta), path("out.fastq.gz"), emit: fastq
    path "*_fastqc/fastqc_data.txt", emit: fastqc_data
    path "fastq_line_count.txt", emit: line_count
    path "sampleinfo.txt", emit: sampleinfo
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo "file_paths: ${meta.fpaths}" > sampleinfo.txt
    echo "copying fastq files"
    for fpath in ${meta.fpaths.join(" ")}; do
        echo "copying \$fpath"
        ls -l \$fpath
        cp \$fpath .
    done
    ls -ld
    echo "copying done"
    gunzip *.gz
    wc -l *.fastq > fastq_line_count.txt

    echo "combining fastq files and gzipping"
    cat *.fastq | pv | pigz -p ${task.cpus} > out.fastq.gz

    echo "running fastqc"
    fastqc --extract out.fastq.gz

    echo "collecting version info"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(echo \$(fastqc --version 2>&1) | sed 's/^.*FastQC //')
    END_VERSIONS
    """

    stub:
    """
    touch out.fastq.gz
    touch fastq_line_count.txt
    touch sampleinfo.txt
    touch versions.yml
    mkdir -p ${meta.sample}_fastqc
    touch ${meta.sample}_fastqc/fastqc_data.txt
    """
}
