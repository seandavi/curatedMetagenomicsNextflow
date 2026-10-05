/*
 * HUMAnN functional profiling (ADR-0016)
 *
 * Two per-sample processes, run on the full-depth host-decontaminated reads
 * only, driven by the selected bundle in conf/humann_bundles.config
 * (params.humann_bundle):
 *
 *   metaphlan_for_humann  the bundle's own MetaPhlAn pass (version + index
 *                         HUMAnN can read; independent of the main taxonomy
 *                         pass) producing the --taxonomic-profile
 *   humann                HUMAnN in the bundle's container, consuming it
 *
 * Containers, resources and maxForks are set in the process bodies (as in
 * kraken.nf): the container is bundle-dependent, so it is a closure resolved
 * when a task launches. Output declarations use globs and HUMAnN's native
 * filenames so bundles with a different file set (e.g. HUMAnN 4's
 * `_2_genefamilies` / `_3_reactions` / `_4_pathabundance`) need no change.
 *
 * Published under <sample>/humann/<bundle>/ with the profile that drove the
 * stratification under <sample>/humann/<bundle>/metaphlan/.
 */

process metaphlan_for_humann {
    container { params.humann_bundles[params.humann_bundle].metaphlan_container }

    label 'profiling'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/humann/${params.humann_bundle}/metaphlan", pattern: "{*.tsv,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 16
    memory { 32.GB * task.attempt }

    input:
    val meta
    path fastq
    path metaphlan_db

    output:
    val(meta), emit: meta
    tuple val(meta), path("metaphlan_rel_ab_w_read_stats.tsv"), emit: profile
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch metaphlan_rel_ab_w_read_stats.tsv
    touch .command.run
    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan_humann: stub
        bowtie2_humann: stub
    END_VERSIONS
    """

    script:
    def bundle = params.humann_bundles[params.humann_bundle]
    """
    metaphlan --input_type fastq \\
        --index ${bundle.metaphlan_index} \\
        ${bundle.metaphlan_db_option} ${metaphlan_db} \\
        --nproc ${task.cpus} \\
        -t rel_ab_w_read_stats \\
        -o metaphlan_rel_ab_w_read_stats.tsv \\
        ${fastq}

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan_humann: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2_humann: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}

process humann {
    container { params.humann_bundles[params.humann_bundle].humann_container }

    label 'functional_profile'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/humann/${params.humann_bundle}", pattern: "{out_*.tsv.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 16
    memory { 48.GB * task.attempt }
    maxForks params.humann_maxforks

    input:
    tuple val(meta), path(fastq), path(taxonomic_profile)
    path chocophlan_db
    path uniref_db
    path utility_mapping_db

    output:
    val(meta), emit: meta
    path "out_*.tsv.gz", emit: tables
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch out_genefamilies.tsv.gz
    touch out_pathabundance.tsv.gz
    touch out_pathcoverage.tsv.gz
    touch .command.run
    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: stub
    END_VERSIONS
    """

    script:
    def bundle = params.humann_bundles[params.humann_bundle]
    def utility_db = bundle.humann_utility_db_option ? "${bundle.humann_utility_db_option} ${utility_mapping_db}" : ''
    """
    humann -i ${fastq} \\
        -o '.' \\
        --output-basename out \\
        --verbose \\
        --nucleotide-database ${chocophlan_db} \\
        --taxonomic-profile ${taxonomic_profile} \\
        --protein-database ${uniref_db} \\
        ${utility_db} \\
        --threads ${task.cpus}

    # Renormalize the gene-family and pathway-abundance tables whatever the
    # HUMAnN release named them (out_genefamilies, out_2_genefamilies, ...).
    for table in out_*genefamilies.tsv out_*pathabundance.tsv; do
        base=\${table%.tsv}
        humann_renorm_table --input \$table --output \${base}_cpm.tsv --units cpm
        humann_renorm_table --input \$table --output \${base}_relab.tsv --units relab
    done

    # Split every table (including pathcoverage/reactions, where present) into
    # stratified and unstratified halves.
    for table in out_*.tsv; do
        humann_split_stratified_table -i \$table -o .
    done

    gzip out_*.tsv

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}
