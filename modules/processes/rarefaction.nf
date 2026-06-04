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
    // is later imported under an alias. seqtk sample is light; the headroom is
    // for large inputs. Memory escalates on retry (Alpine caps time at 24h, so
    // no time escalation — memory only).
    cpus 2
    memory { 8.GB * task.attempt }

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
