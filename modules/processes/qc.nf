/*
 * Per-sample QC reporting (ADR-0008)
 *
 * fastqc: FastQC on the host-decontaminated reads. Raw-read FastQC already runs
 *   in fasterq_dump/local_fastqc, so this adds a post-QC view of the reads that
 *   actually feed profiling. Runs in the base image (FastQC is already present),
 *   so no additional container is required.
 *
 * Per-sample only.
 */

process fastqc {
    label 'qc'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/fastqc", pattern: "{*_fastqc.zip,*_fastqc/fastqc_data.txt,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 2
    memory { 4.GB * task.attempt }

    input:
    val meta
    path reads

    output:
    val(meta), emit: meta
    path "*_fastqc/fastqc_data.txt", emit: fastqc_data
    path "*_fastqc.zip"
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    mkdir -p decontaminated_fastqc
    touch decontaminated_fastqc/fastqc_data.txt
    touch decontaminated_fastqc.zip
    touch .command.run
    touch versions.yml
    """

    script:
    """
    # Name the FastQC output after the sample stage so the report is clearly the
    # decontaminated reads rather than a generic "out".
    ln -s ${reads} decontaminated.fastq
    fastqc --extract decontaminated.fastq

    cat <<-END_VERSIONS > versions.yml
    versions:
        fastqc: \$( echo \$(fastqc --version 2>&1 ) | sed 's/^.*FastQC //' )
    END_VERSIONS
    """
}
