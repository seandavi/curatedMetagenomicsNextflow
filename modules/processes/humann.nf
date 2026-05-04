/*
 * Functional profiling process
 */

process humann {
    label 'functional_profile'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/humann", mode: "${params.publish_mode}"
    cpus 16
    memory "48g"

    tag "${meta.sample}"

    input:
    val meta
    path fastq
    path marker_rel_ab_w_read_stats
    path chocophlan_db
    path uniref_db
    path utility_mapping_db

    output:
    val(meta), emit: meta
    path("out_genefamilies.tsv.gz")
    path("out_genefamilies_cpm.tsv.gz")
    path("out_genefamilies_relab.tsv.gz")
    path("out_genefamilies_stratified.tsv.gz")
    path("out_genefamilies_unstratified.tsv.gz")
    path("out_genefamilies_cpm_stratified.tsv.gz")
    path("out_genefamilies_relab_stratified.tsv.gz")
    path("out_genefamilies_cpm_unstratified.tsv.gz")
    path("out_genefamilies_relab_unstratified.tsv.gz")
    path("out_pathabundance.tsv.gz")
    path("out_pathabundance_cpm.tsv.gz")
    path("out_pathabundance_relab.tsv.gz")
    path("out_pathabundance_stratified.tsv.gz")
    path("out_pathabundance_unstratified.tsv.gz")
    path("out_pathabundance_cpm_stratified.tsv.gz")
    path("out_pathabundance_relab_stratified.tsv.gz")
    path("out_pathabundance_cpm_unstratified.tsv.gz")
    path("out_pathabundance_relab_unstratified.tsv.gz")
    path("out_pathcoverage_unstratified.tsv.gz")
    path("out_pathcoverage_stratified.tsv.gz")
    path("out_pathcoverage.tsv.gz")
    path ".command*"
    path "versions.yml"

    stub:
    """
    touch out_genefamilies.tsv.gz
    touch out_genefamilies_cpm.tsv.gz
    touch out_genefamilies_relab.tsv.gz
    touch out_genefamilies_stratified.tsv.gz
    touch out_genefamilies_unstratified.tsv.gz
    touch out_genefamilies_cpm_stratified.tsv.gz
    touch out_genefamilies_relab_stratified.tsv.gz
    touch out_genefamilies_cpm_unstratified.tsv.gz
    touch out_genefamilies_relab_unstratified.tsv.gz
    touch out_pathabundance.tsv.gz
    touch out_pathabundance_cpm.tsv.gz
    touch out_pathabundance_relab.tsv.gz
    touch out_pathabundance_stratified.tsv.gz
    touch out_pathabundance_unstratified.tsv.gz
    touch out_pathabundance_cpm_stratified.tsv.gz
    touch out_pathabundance_relab_stratified.tsv.gz
    touch out_pathabundance_cpm_unstratified.tsv.gz
    touch out_pathabundance_relab_unstratified.tsv.gz
    touch out_pathcoverage_unstratified.tsv.gz
    touch out_pathcoverage_stratified.tsv.gz
    touch out_pathcoverage.tsv.gz
    touch .command.run
    touch versions.yml
    """

    script:
    """
    humann -i ${fastq} \
        -o '.' \
        --verbose \
        --nucleotide-database ${chocophlan_db} \
        --taxonomic-profile ${marker_rel_ab_w_read_stats} \
        --protein-database ${uniref_db} \
        --utility-database ${utility_mapping_db} \
        --threads ${task.cpus}

    humann_renorm_table \
        --input out_pathabundance.tsv \
        --output out_pathabundance_cpm.tsv \
        --units cpm

    humann_renorm_table \
        --input out_genefamilies.tsv \
        --output out_genefamilies_cpm.tsv \
        --units cpm

    humann_renorm_table \
        --input out_genefamilies.tsv \
        --output out_genefamilies_relab.tsv \
        --units relab

    humann_renorm_table \
        --input out_pathabundance.tsv \
        --output out_pathabundance_relab.tsv \
        --units relab

    humann_split_stratified_table -i out_pathabundance.tsv -o .
    humann_split_stratified_table -i out_pathabundance_cpm.tsv -o .
    humann_split_stratified_table -i out_pathabundance_relab.tsv -o .
    humann_split_stratified_table -i out_pathcoverage.tsv -o .
    humann_split_stratified_table -i out_genefamilies.tsv -o .
    humann_split_stratified_table -i out_genefamilies_cpm.tsv -o .
    humann_split_stratified_table -i out_genefamilies_relab.tsv -o .
    gzip out_*tsv

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}
