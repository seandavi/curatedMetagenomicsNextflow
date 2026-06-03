/*
 * Complementary read-based profiling: Kraken2 + Bracken (ADR-0006)
 *
 * Kraken2 classifies the host-decontaminated reads against a prebuilt PlusPF
 * index; Bracken then re-estimates abundances at species and genus level. This
 * complements MetaPhlAn's marker-based profiling with a whole-community,
 * read-count-based view.
 *
 * Resource, container, and concurrency directives live in the process bodies
 * (not conf/base.config) because these processes are imported under aliases for
 * the full and rarefied branches; in-body directives apply regardless of alias.
 *
 * The database is loaded into RAM (default Kraken2 mode — no --memory-mapping
 * and no staging to node-local scratch) for cluster portability. maxForks
 * throttles how many tasks read the DB off shared storage at once.
 */

process kraken2 {
    container 'docker://staphb/kraken2:2.1.3'

    label 'profiling'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/${meta.branch ? meta.branch + '/' : ''}kraken", pattern: "{*.report.txt.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 8
    memory { 32.GB * task.attempt }
    maxForks params.kraken_maxforks

    input:
    val meta
    path reads
    path kraken_db

    output:
    val(meta), emit: meta
    path "kraken2.report.txt", emit: report
    path "kraken2.report.txt.gz", emit: report_gz
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch kraken2.report.txt
    touch kraken2.report.txt.gz
    touch .command.run
    touch versions.yml
    """

    script:
    """
    # Default mode loads the DB into RAM; the per-read output is discarded
    # (Bracken consumes the report, not the per-read assignments).
    kraken2 \
        --db ${kraken_db} \
        --threads ${task.cpus} \
        --confidence ${params.kraken_confidence} \
        --report kraken2.report.txt \
        --output /dev/null \
        ${reads}

    gzip -c kraken2.report.txt > kraken2.report.txt.gz

    cat <<-END_VERSIONS > versions.yml
    versions:
        kraken2: \$( kraken2 --version 2>&1 | head -1 | sed 's/Kraken version //' )
    END_VERSIONS
    """
}

process bracken {
    container 'docker://staphb/bracken:2.9'

    label 'profiling'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/${meta.branch ? meta.branch + '/' : ''}kraken", pattern: "bracken*.txt.gz", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 2
    memory { 4.GB * task.attempt }

    input:
    val meta
    path report
    path kraken_db

    output:
    val(meta), emit: meta
    path "bracken.species.txt.gz", emit: bracken_species
    path "bracken.genus.txt.gz", emit: bracken_genus
    path "bracken.species.report.txt.gz"
    path "bracken.genus.report.txt.gz"
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch bracken.species.txt.gz
    touch bracken.genus.txt.gz
    touch bracken.species.report.txt.gz
    touch bracken.genus.report.txt.gz
    touch .command.run
    touch versions.yml
    """

    script:
    """
    bracken -d ${kraken_db} -i ${report} \
        -o bracken.species.txt -w bracken.species.report.txt \
        -r ${params.bracken_read_length} -l S
    bracken -d ${kraken_db} -i ${report} \
        -o bracken.genus.txt -w bracken.genus.report.txt \
        -r ${params.bracken_read_length} -l G

    gzip bracken.species.txt bracken.species.report.txt \
         bracken.genus.txt bracken.genus.report.txt

    cat <<-END_VERSIONS > versions.yml
    versions:
        bracken: \$( echo \$(bracken -v 2>&1) | awk '{print \$NF}' )
    END_VERSIONS
    """
}
