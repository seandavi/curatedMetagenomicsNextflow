/*
 * Resistome profiling: KMA against CARD (ADR-0012, supersedes ADR-0007)
 *
 * KMA (k-mer alignment) maps the host-decontaminated reads directly against a
 * KMA-indexed CARD reference (the homolog-model nucleotide FASTA; see
 * card_kma_db in databases.nf) and reports antimicrobial-resistance gene
 * template hits with depth/coverage statistics, assembly-free.
 *
 * Unlike the previous RGI/CARD resistome (full-branch only), KMA is cheap
 * enough to run on BOTH the full and rarefied branches, matching MetaPhlAn and
 * Kraken. It is imported under per-branch aliases, so container and resource
 * directives live in the process body (not conf/base.config) to apply
 * regardless of alias.
 *
 * NOTE: The KMA invocation is not exercised by stub tests (which do not run
 * containers); the first real run is the true validation. See
 * docs/adr/0012-resistome-kma-card.md.
 */

process resistome_kma {
    container 'docker://quay.io/biocontainers/kma:1.6.13--h118bc1c_0'

    label 'profiling'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/${meta.branch ? meta.branch + '/' : ''}resistome", pattern: "{*.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 8
    memory { 16.GB * task.attempt }

    input:
    val meta
    path reads
    path card_db

    output:
    val(meta), emit: meta
    path "card_kma.res.gz", emit: res
    path "card_kma.mapstat.gz", emit: mapstat
    path "card_kma.aln.gz"
    path "card_kma.fsa.gz"
    path "card_kma.frag.gz"
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch card_kma.res.gz
    touch card_kma.mapstat.gz
    touch card_kma.aln.gz
    touch card_kma.fsa.gz
    touch card_kma.frag.gz
    touch .command.run
    touch versions.yml
    """

    script:
    """
    # Reads are single-end throughout the pipeline. -ef adds the extended
    # mapstat table (fragment counts, depth) that is useful for ARG
    # quantification. KMA writes card_kma.frag.gz already gzipped.
    kma \
        -i ${reads} \
        -o card_kma \
        -t_db ${card_db}/card_kma_db \
        -t ${task.cpus} \
        -ef

    # Compress the remaining plain-text outputs to save space.
    gzip card_kma.res card_kma.fsa card_kma.aln card_kma.mapstat

    cat <<-END_VERSIONS > versions.yml
    versions:
        kma: \$( kma -v 2>&1 | head -n1 | sed 's/^KMA-//' )
    END_VERSIONS
    """
}
