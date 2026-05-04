/*
 * Input preprocessing processes
 *
 * These processes normalize either remote SRA accessions or local FASTQ files
 * into the same output contract: a single gzipped FASTQ plus lightweight QC
 * artifacts used downstream.
 */

process fasterq_dump {
    label 'qc'
    label 'io_heavy'
    label 'download_retry'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/fasterq_dump", pattern: "{fastq_line_count.txt,*_fastqc/fastqc_data.txt,sampleinfo.txt,.command*}", mode: "${params.publish_mode}"

    cpus 8
    memory "16g"

    tag "${meta.sample}"

    input:
    val meta

    output:
    val(meta), emit: meta
    path "out.fastq.gz", emit: fastq
    path "*_fastqc/fastqc_data.txt", emit: fastqc_data
    path "fastq_line_count.txt"
    path ".command*"
    path "sampleinfo.txt"
    path "versions.yml"

    stub:
    """
    touch out.fastq.gz
    touch fastq_line_count.txt
    touch sampleinfo.txt
    touch versions.yml
    touch .command.run
    mkdir -p ${meta.sample}_fastqc
    touch ${meta.sample}_fastqc/fastqc_data.txt
    """

    script:
    """

    echo "accessions: ${meta.accessions}" > sampleinfo.txt
    echo "starting fasterq-dump"
    for accession in ${meta.accessions.join(" ")}; do
        echo "downloading \$accession"
        curl -o \$accession.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/\$accession/\$accession
        fasterq-dump --threads ${task.cpus} \
            --skip-technical \
            --force \
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
    versions:
        awscli: \$(aws --version)
        fastqc: \$( echo \$(fastqc --version 2>&1 ) | sed 's/^.*FastQC //' )
        fasterq-dump: \$( echo \$(fasterq-dump --version 2>&1 ) | head -2 | tail -1 | awk '{print \$3}')
    END_VERSIONS

    echo "finalizing fasterqc-dump"

    """
}

process local_fastqc {
    label 'qc'
    label 'io_heavy'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/local_fastqc", pattern: "{fastq_line_count.txt,*_fastqc/fastqc_data.txt,sampleinfo.txt,.command*}", mode: "${params.publish_mode}"

    cpus 8
    memory "16g"

    tag "${meta.sample}"

    input:
    val meta

    output:
    val(meta), emit: meta
    path "out.fastq.gz", emit: fastq
    path "*_fastqc/fastqc_data.txt", emit: fastqc_data
    path "fastq_line_count.txt"
    path ".command*"
    path "sampleinfo.txt"
    path "versions.yml"

    stub:
    """
    touch out.fastq.gz
    touch fastq_line_count.txt
    touch sampleinfo.txt
    touch versions.yml
    touch .command.run
    mkdir -p ${meta.sample}_fastqc
    touch ${meta.sample}_fastqc/fastqc_data.txt
    """

    script:
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
    versions:
        fastqc: \$( echo \$(fastqc --version 2>&1 ) | sed 's/^.*FastQC //' )
    END_VERSIONS
    """
}
