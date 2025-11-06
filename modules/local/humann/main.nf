process HUMANN {
    tag "${meta.sample}"
    label 'process_high'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'

    input:
    tuple val(meta), path(fastq)
    path metaphlan_unknown_list
    path chocophlan_db
    path uniref_db

    output:
    path "out_genefamilies.tsv.gz", emit: genefamilies
    path "out_genefamilies_cpm.tsv.gz", emit: genefamilies_cpm
    path "out_genefamilies_relab.tsv.gz", emit: genefamilies_relab
    path "out_genefamilies_stratified.tsv.gz", emit: genefamilies_stratified
    path "out_genefamilies_unstratified.tsv.gz", emit: genefamilies_unstratified
    path "out_genefamilies_cpm_stratified.tsv.gz", emit: genefamilies_cpm_stratified
    path "out_genefamilies_relab_stratified.tsv.gz", emit: genefamilies_relab_stratified
    path "out_genefamilies_cpm_unstratified.tsv.gz", emit: genefamilies_cpm_unstratified
    path "out_genefamilies_relab_unstratified.tsv.gz", emit: genefamilies_relab_unstratified
    path "out_pathabundance.tsv.gz", emit: pathabundance
    path "out_pathabundance_cpm.tsv.gz", emit: pathabundance_cpm
    path "out_pathabundance_relab.tsv.gz", emit: pathabundance_relab
    path "out_pathabundance_stratified.tsv.gz", emit: pathabundance_stratified
    path "out_pathabundance_unstratified.tsv.gz", emit: pathabundance_unstratified
    path "out_pathabundance_cpm_stratified.tsv.gz", emit: pathabundance_cpm_stratified
    path "out_pathabundance_relab_stratified.tsv.gz", emit: pathabundance_relab_stratified
    path "out_pathabundance_cpm_unstratified.tsv.gz", emit: pathabundance_cpm_unstratified
    path "out_pathabundance_relab_unstratified.tsv.gz", emit: pathabundance_relab_unstratified
    path "out_pathcoverage_unstratified.tsv.gz", emit: pathcoverage_unstratified
    path "out_pathcoverage_stratified.tsv.gz", emit: pathcoverage_stratified
    path "out_pathcoverage.tsv.gz", emit: pathcoverage
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    humann -i ${fastq} \\
        -o '.' \\
        --verbose \\
        --metaphlan-options "-t rel_ab --index latest" \\
        --nucleotide-database ${chocophlan_db} \\
        --taxonomic-profile ${metaphlan_unknown_list} \\
        --protein-database ${uniref_db} \\
        --threads ${task.cpus}

    humann_renorm_table \\
        --input out_pathabundance.tsv \\
        --output out_pathabundance_cpm.tsv \\
        --units cpm

    humann_renorm_table \\
        --input out_genefamilies.tsv \\
        --output out_genefamilies_cpm.tsv \\
        --units cpm

    humann_renorm_table \\
        --input out_genefamilies.tsv \\
        --output out_genefamilies_relab.tsv \\
        --units relab

    humann_renorm_table \\
        --input out_pathabundance.tsv \\
        --output out_pathabundance_relab.tsv \\
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
    "${task.process}":
        humann: \$(echo \$(humann --version 2>&1) | awk '{print \$2}')
    END_VERSIONS
    """

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
    touch versions.yml
    """
}
