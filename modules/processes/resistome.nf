/*
 * Resistome profiling: RGI against CARD (ADR-0007)
 *
 * RGI's read-based `bwt` workflow aligns the host-decontaminated reads against
 * the CARD reference and reports antimicrobial-resistance gene and allele
 * abundances, assembly-free. Runs on the full-depth branch only (ARGs are rare;
 * the rarefied depth yields a sparse, low-value resistome).
 *
 * Container and resources are set in the process body (not conf/base.config) to
 * stay consistent with the other tool-specific processes that override the base
 * image (e.g. kraken2/bracken).
 *
 * NOTE: The RGI load + bwt invocation follows the CARD documentation but is not
 * exercised by the stub tests (which do not run containers); the first real run
 * is the true validation. See docs/adr/0007-resistome-rgi-card.md.
 */

process resistome {
    container 'docker://quay.io/biocontainers/rgi:6.0.5--pyha8f3691_0'

    label 'profiling'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/${meta.branch ? meta.branch + '/' : ''}resistome", pattern: "{*.txt.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 8
    memory { 16.GB * task.attempt }

    input:
    val meta
    path reads
    path card_db

    output:
    val(meta), emit: meta
    path "rgi_bwt.gene_mapping_data.txt.gz", emit: gene_mapping
    path "rgi_bwt.allele_mapping_data.txt.gz", emit: allele_mapping
    path "rgi_bwt.*mapping_stats.txt.gz"
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch rgi_bwt.gene_mapping_data.txt.gz
    touch rgi_bwt.allele_mapping_data.txt.gz
    touch rgi_bwt.overall_mapping_stats.txt.gz
    touch .command.run
    touch versions.yml
    """

    script:
    """
    # Load CARD canonical data into a run-local database (--local writes a
    # localDB/ in the work dir, so no writable package install is required).
    rgi load --card_json ${card_db}/card.json --local

    # Build and load the bwt reference annotation; card_annotation writes a
    # version-stamped FASTA (card_database_v<CARD version>.fasta).
    rgi card_annotation -i ${card_db}/card.json
    card_fasta=\$(ls card_database_v*.fasta | grep -v '_all' | head -n 1)
    rgi load --card_json ${card_db}/card.json --card_annotation "\$card_fasta" --local

    rgi bwt \
        --read_one ${reads} \
        --aligner ${params.rgi_aligner} \
        --threads ${task.cpus} \
        --output_file rgi_bwt \
        --local \
        --clean

    gzip rgi_bwt.*.txt

    cat <<-END_VERSIONS > versions.yml
    versions:
        rgi: \$( rgi main --version 2>&1 | tail -n 1 )
    END_VERSIONS
    """
}
