/*
 * Taxonomic and strain-level profiling processes
 */

process kneaddata {
    label 'host_filter'

    publishDir "${params.publish_dir}/${meta.sample}/kneaddata", pattern: "{kneaddata_output/kneaddata_fastq_linecounts.txt,kneaddata_output/out_kneaddata.log,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 8
    memory "30g"

    input:
    val meta
    path fastq
    path kd_genome
    path kd_mouse

    output:
    val(meta), emit: meta
    path "kneaddata_output/out.fastq", emit: fastq
    path "kneaddata_output/kneaddata_fastq_linecounts.txt"
    path "kneaddata_output/out_kneaddata.log"
    path ".command*"

    stub:
    """
    mkdir -p kneaddata_output
    touch kneaddata_output/out.fastq
    touch kneaddata_output/kneaddata_fastq_linecounts.txt
    touch kneaddata_output/out_kneaddata.log
    """

    script:
    """
    kneaddata --unpaired ${fastq} \
        --reference-db ${params.organism_database} \
        --output kneaddata_output  \
        --trimmomatic /installed/Trimmomatic-0.39 \
        --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:30' \
        --bypass-trf \
	--bowtie2-options='--very-fast' \
	-t ${task.cpus} -p ${task.cpus}

    cd kneaddata_output
    cat out_kneaddata.fastq | sed 's/^+.RR.*/+/g' > out.fastq
    rm out_kneaddata.fastq
    wc -l * | grep fastq > kneaddata_fastq_linecounts.txt
    """
}

process metaphlan_unknown_viruses_lists {
    label 'profiling'

    publishDir "${params.publish_dir}/${meta.sample}/metaphlan_lists", pattern: "{*tsv.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 16
    memory "47g"

    input:
    val meta
    path fastq
    path metaphlan_db


    output:
    val(meta), emit: meta
    path 'bowtie2.out.gz', emit: metaphlan_bt2
    path 'metaphlan_unknown_list.tsv', emit: metaphlan_unknown_list
    path 'metaphlan_unknown_list.tsv.gz', emit: metaphlan_unknown_list_gz
    path 'metaphlan_viruses_list.tsv', emit: metaphlan_viruses_list
    path 'metaphlan_viruses_list.tsv.gz', emit: metaphlan_viruses_list_gz
    path 'metaphlan.sam', emit: metaphlan_sam
    path ".command*"
    path "versions.yml"

    stub:
    """
    touch bowtie2.out.gz
    touch metaphlan.sam
    touch metaphlan_unknown_list.tsv
    touch metaphlan_unknown_list.tsv.gz
    touch metaphlan_viruses_list.tsv
    touch metaphlan_viruses_list.tsv.gz
    touch .command.run
    touch versions.yml
    """


    script:
    """
    find .
    metaphlan --input_type fastq \
        --index ${params.metaphlan_index} \
        --db_dir metaphlan \
        --mapout bowtie2.out \
        --nproc ${task.cpus} \
        --profile_vsc \
        -s metaphlan.sam \
        --vsc_breadth 0.75 \
        --vsc_out metaphlan_viruses_list.tsv \
        -o metaphlan_unknown_list.tsv \
        ${fastq}

    gzip -c metaphlan_unknown_list.tsv > metaphlan_unknown_list.tsv.gz
    gzip -c metaphlan_viruses_list.tsv > metaphlan_viruses_list.tsv.gz
    gzip bowtie2.out

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}

process metaphlan_unknown_list {
    label 'profiling'

    publishDir "${params.publish_dir}/${meta.sample}/metaphlan_lists", pattern: "{*tsv.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 16
    memory "30g"

    input:
    val meta
    path metaphlan_bt2
    path metaphlan_db


    output:
    val(meta), emit: meta
    path 'metaphlan_unknown_list.tsv', emit: metaphlan_unknown_list
    path 'metaphlan_unknown_list.tsv.gz', emit: metaphlan_unknown_list_gz
    path ".command*"
    path "versions.yml"

    stub:
    """
    touch metaphlan_unknown_list.tsv
    touch metaphlan_unknown_list.tsv.gz
    touch .command.run
    touch versions.yml
    """


    script:
    """
    metaphlan \
        --input_type mapout \
        --index ${params.metaphlan_index} \
        --db_dir metaphlan \
        --nproc ${task.cpus} \
        -o metaphlan_unknown_list.tsv \
        <( gunzip -c ${metaphlan_bt2} )

    gzip -c metaphlan_unknown_list.tsv > metaphlan_unknown_list.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}

process metaphlan_markers {
    label 'profiling'

    publishDir "${params.publish_dir}/${meta.sample}/metaphlan_markers/", pattern: "{*tsv.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 16
    memory "30g"

    input:
    val meta
    path metaphlan_bt2
    path metaphlan_db

    output:
    val meta, emit: meta
    path "marker_abundance.tsv.gz", emit: marker_abundance
    path "marker_presence.tsv.gz", emit: marker_presence
    path "marker_rel_ab_w_read_stats.tsv.gz", emit: marker_rel_ab_w_read_stats_gz
    path "marker_rel_ab_w_read_stats.tsv", emit: marker_rel_ab_w_read_stats
    path ".command*"
    path "versions.yml"

    stub:
    """
    touch marker_abundance.tsv.gz
    touch marker_presence.tsv.gz
    touch marker_rel_ab_w_read_stats.tsv.gz
    touch marker_rel_ab_w_read_stats.tsv
    touch .command.run
    touch versions.yml
    """

    script:
    """
    metaphlan --input_type mapout \
        --index ${params.metaphlan_index} \
        --db_dir metaphlan \
        -t marker_pres_table \
        --nproc ${task.cpus} \
        -o marker_presence.tsv \
        <( gunzip -c ${metaphlan_bt2} )
    metaphlan --input_type mapout \
        --index ${params.metaphlan_index} \
        --db_dir metaphlan \
        --nproc ${task.cpus} \
        -t marker_ab_table \
        -o marker_abundance.tsv \
        <( gunzip -c ${metaphlan_bt2} )
    metaphlan --input_type mapout \
        --index ${params.metaphlan_index} \
        --db_dir metaphlan \
        -t rel_ab_w_read_stats \
        --nproc ${task.cpus} \
        -o marker_rel_ab_w_read_stats.tsv \
        <( gunzip -c ${metaphlan_bt2} )
    gzip *.tsv
    gunzip -c marker_rel_ab_w_read_stats.tsv.gz > marker_rel_ab_w_read_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}

process sample_to_markers {
    label 'profiling'

    publishDir "${params.publish_dir}/${meta.sample}/strainphlan_markers/", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 16
    memory "32g"

    input:
    val meta
    path metaphlan_sam
    path metaphlan_db

    output:
    val meta, emit: meta
    path "metaphlan.json.bz2", emit: sample_to_markers_json_bz2
    path ".command*"
    path "versions.yml"

    stub:
    """
    touch metaphlan.json.bz2
    touch .command.run
    touch versions.yml
    """

    script:
    """
    sample2markers.py \
        --input ${metaphlan_sam} \
        --input_format sam \
        --database ${metaphlan_db}/${params.metaphlan_index}.pkl \
        --nprocs ${task.cpus} \
        --output_dir .

    cat <<-END_VERSIONS > versions.yml
    versions:
        sample2markers.py: \$( echo \$(sample2markers.py --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}
