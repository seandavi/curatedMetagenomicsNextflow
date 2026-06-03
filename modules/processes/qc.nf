/*
 * Per-sample QC reporting (ADR-0008)
 *
 * fastqc: FastQC on the host-decontaminated reads (raw-read FastQC already runs
 *   in fasterq_dump/local_fastqc), giving a post-QC view of the reads that feed
 *   profiling.
 * multiqc: aggregate the raw and post-decontamination FastQC into one
 *   per-sample report. The two FastQC outputs are staged into raw/ and clean/
 *   subdirectories and MultiQC is run with --dirs so their (otherwise identical)
 *   sample names are disambiguated.
 *
 * Per-sample only — no cross-sample aggregation.
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

process multiqc {
    container 'docker://quay.io/biocontainers/multiqc:1.35--pyhdfd78af_0'

    label 'qc'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/multiqc", pattern: "{multiqc_report.html,multiqc_report_data/**,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 1
    memory { 4.GB * task.attempt }

    input:
    tuple val(meta),
          path(raw_fastqc, stageAs: 'raw/fastqc_data.txt'),
          path(clean_fastqc, stageAs: 'clean/fastqc_data.txt')

    output:
    val(meta), emit: meta
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", type: 'dir'
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    touch .command.run
    touch versions.yml
    """

    script:
    """
    # --dirs prepends the raw/ and clean/ directory names to the FastQC sample
    # names so the two reports (both named for the same read file) stay distinct.
    multiqc --dirs --dirs-depth 1 --force --no-ansi -n multiqc_report .

    cat <<-END_VERSIONS > versions.yml
    versions:
        multiqc: \$( multiqc --version 2>&1 | awk '{print \$NF}' )
    END_VERSIONS
    """
}
