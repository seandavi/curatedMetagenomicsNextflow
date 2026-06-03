/*
 * Per-sample manifest
 *
 * Compiles a single machine-readable manifest.json at the sample root with
 * provenance (pipeline/Nextflow/container versions, command line, parameters,
 * inputs), read accounting (raw vs. host-decontaminated read/base/length/GC
 * statistics), the rarefaction parameters when that branch is active, and the
 * consolidated per-process tool versions.
 *
 * Read statistics are computed in pure Python (bin/build_manifest.py) so this
 * step needs no tool beyond the base image. See
 * docs/adr/0005-per-sample-manifest.md.
 */

process sample_manifest {
    label 'finalize'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}", pattern: "manifest.json", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 1
    memory "2g"

    input:
    tuple val(meta),
          path(raw_fastq),
          path(raw_versions, stageAs: 'versions_raw.yml'),
          path(kneaddata_fastq, stageAs: 'kneaddata.fastq'),
          path(kneaddata_versions, stageAs: 'versions_kneaddata.yml'),
          path(metaphlan_versions, stageAs: 'versions_metaphlan.yml')

    output:
    val(meta), emit: meta
    path "manifest.json", emit: manifest
    path ".command*"

    stub:
    """
    echo '{"sample_id": "${meta.sample}", "stub": true}' > manifest.json
    touch .command.run
    """

    script:
    def cmd = workflow.commandLine ? workflow.commandLine.replaceAll('"', "'") : ''
    def ids = meta.accessions ? meta.accessions.join(',') : (meta.fpaths ? meta.fpaths.join(',') : '')
    def mode = params.local_input ? 'local' : 'sra'
    """
    build_manifest.py \
        --sample '${meta.sample}' \
        --pipeline-name '${workflow.manifest.name}' \
        --pipeline-version '${workflow.manifest.version}' \
        --nextflow-version '${workflow.nextflow.version}' \
        --container '${task.container ?: ""}' \
        --command-line "${cmd}" \
        --run-name '${params.run_name ?: ""}' \
        --session-id '${workflow.sessionId ?: ""}' \
        --git-commit '${workflow.commitId ?: ""}' \
        --git-revision '${workflow.revision ?: ""}' \
        --start-time '${workflow.start ?: ""}' \
        --input-mode '${mode}' \
        --input-ids '${ids}' \
        --metaphlan-index '${params.metaphlan_index}' \
        --store-dir '${params.store_dir}' \
        --skip-rarefied ${params.skip_rarefied} \
        --rarefy-reads ${params.rarefy_reads} \
        --rarefy-seed ${params.rarefy_seed} \
        --skip-humann ${params.skip_humann} \
        --skip-gtdb ${params.skip_gtdb} \
        --raw-fastq ${raw_fastq} \
        --kneaddata-fastq ${kneaddata_fastq} \
        --output manifest.json
    """
}
