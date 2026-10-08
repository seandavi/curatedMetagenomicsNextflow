/*
 * Rarefaction process — downsample FASTQ reads to a fixed depth using seqtk.
 *
 * Inputs and outputs mirror the kneaddata convention so that the rarefied
 * fastq channel can be dropped in wherever the full fastq channel is used.
 */

process rarefy_fastq {
    label 'qc'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/rarefied_data/rarefaction", pattern: "{rarefied.fastq,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    // Set in-body (like the sibling qc processes) rather than via a
    // conf/base.config withName, so the directive holds even if this process
    // is later imported under an alias. seqtk sample is single-threaded and
    // light (peak RSS 0.66 GiB in 2.2.x traces); memory and time escalate on
    // retry (1h x 4 attempts stays under Alpine's 24h cap). Sizes from
    // nextflow_telemetry docs/research/resource-tuning-2.3.0.md, like every
    // first-attempt request in this release.
    cpus 1
    memory { 1792.MB * task.attempt }
    time { 1.h * task.attempt }

    input:
    val meta
    path fastq

    output:
    val(meta), emit: meta
    path "rarefied.fastq", emit: fastq
    path ".command*"

    stub:
    """
    touch rarefied.fastq
    touch .command.run
    """

    script:
    """
    seqtk sample -s ${params.rarefy_seed} ${fastq} ${params.rarefy_reads} > rarefied.fastq
    """
}
