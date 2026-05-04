/*
 * Finalization process
 *
 * This sentinel file is intentionally simple and stable because downstream
 * automation may use it to detect per-sample completion.
 */

process MARK_COMPLETE {
    label 'finalize'

    publishDir "${params.publish_dir}/${meta.sample}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 1
    memory "2g"

    input:
    val meta

    output:
    path "MARK_COMPLETE"

    stub:
    """
    echo "${meta.sample} \$(date -u +%Y-%m-%dT%H:%M:%SZ)" > MARK_COMPLETE
    """

    script:
    """
    echo "${meta.sample} \$(date -u +%Y-%m-%dT%H:%M:%SZ)" > MARK_COMPLETE
    """
}
