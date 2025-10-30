process FASTERQ_DUMP {
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
    echo "accessions: ${meta.accessions}" > sampleinfo.txt
    echo "starting fasterq-dump"
    for accession in ${meta.accessions.join(" ")}; do
        echo "downloading \$accession"
        curl -o \$accession.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/\$accession/\$accession  
        fasterq-dump --threads ${task.cpus} \\
            --skip-technical \\
            --force \\
            --split-files \$accession.sra 
    done
    ls -ld
    echo "fasterq-dump done"
    wc -l *.fastq > fastq_line_count.txt

    echo "combining fastq files and gzipping"
    cat *.fastq | pv | pigz -p ${task.cpus} > out.fastq.gz

    echo "running fastqc"
    fastqc --extract out.fastq.gz

    echo "collecting version info"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awscli: \$(aws --version 2>&1 | sed 's/aws-cli\\///' | awk '{print \$1}')
        fastqc: \$(echo \$(fastqc --version 2>&1) | sed 's/^.*FastQC //')
        fasterq-dump: \$(echo \$(fasterq-dump --version 2>&1) | head -2 | tail -1 | awk '{print \$3}')
    END_VERSIONS

    echo "finalizing fasterqc-dump"
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
