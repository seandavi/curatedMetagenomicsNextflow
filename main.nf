#!/usr/bin/env nextflow
// params.runs
// params.metaphlan_index="latest"
// params.chocophlan_db
// params.uniref_db
// params.

nextflow.preview.dsl=2
params.publish_dir="gs://temp-testing/nf_testing/${params.runs}/"
params.uniref="uniref50_diamond"
params.chocophlan="full"

process fasterq_dump {
    publishDir "${params.publish_dir}"

    tag "$srr"
    
    input:
      val srr
    output:
    path "*.fastq", emit: fastq_files
    path "*_fastqc/fastqc_data.txt", emit: fastqc_data
    script:
      """
      fasterq-dump \
          --skip-technical \
          --force \
          --threads ${task.cpus} \
          --split-files ${srr} 
      fastqc --extract *.fastq
      """
}


process concat_fastq {
    input:
      path x
    output:
    path 'out.fastq', emit: fastq
      path 'wordcount.fastq'
    script:
      """
      cat $x >> out.fastq
      wc ${x} > wordcount.fastq
      """
}


process install_metaphlan_db {
    storeDir "$HOME/bigdb/metaphlan/"

    output:
    path 'mpa_v30_CHOCOPhlAn_201901.*', emit: metaphlan_db
    
    script:
      """
      metaphlan --install --index latest --bowtie2db .
      """
}

process metaphlan_bugs_list {
    publishDir "${params.publish_dir}/metaphlan"

    cpus 8
    
//    input:
//    path metaphlan_db
    
    input:
    path fastq

    
    output:
    path 'sam.bz2', emit: metaphlan_sam
    path 'bowtie2.out', emit: metaphlan_bt2
    path 'metaphlan_bugs_list.tsv', emit: metaphlan_bugs_list
    
    script:
    """
    metaphlan --input_type fastq \
        --index ${params.metaphlan_index} \
        --bowtie2db /Users/sdavis2/bigdb/metaphlan/ \
        --samout sam.bz2 \
        --bowtie2out bowtie2.out \
        --nproc ${task.cpus} \
        -o metaphlan_bugs_list.tsv \
        ${fastq}
    """
}

process metaphlan_markers {
    publishDir "${params.publish_dir}/metaphlan"

    cpus 8
    
//    input:
//    path metaphlan_db
    
    input:
    path metaphlan_bt2

    
    output:
    path "marker_abundance.tsv", emit: marker_abundance
    path "marker_presence.tsv", emit: marker_presence

    script:
    """
    metaphlan --input_type bowtie2out \
        --index ${params.metaphlan_index} \
        --bowtie2db /Users/sdavis2/bigdb/metaphlan/ \
        -t marker_pres_table \
        -o marker_presence.tsv \
        ${metaphlan_bt2}
    metaphlan --input_type bowtie2out \
        --index ${params.metaphlan_index} \
        --bowtie2db /Users/sdavis2/bigdb/metaphlan/ \
        -t marker_ab_table \
        -o marker_abundance.tsv \
        ${metaphlan_bt2}
    """
}


process chocophlan_db {
    storeDir "$HOME/bigdb/"

    output:
    path "chocophlan/*", emit: chocophlan_db

    script:
    """
    humann_databases --download chocophlan ${params.chocophlan} .
    """
}


process uniref_db {
    storeDir "$HOME/bigdb/"

    output:
    path "${params.uniref}/*", emit: uniref_db

    script:
    """
    humann_databases --download uniref ${params.uniref} ${params.uniref}
    """
}


process humann {
    publishDir 'humann'
    cpus 8

    input:
    path fastq
    path metaphlan_bugs_list // metaphlan_bugs_list.tsv

    output:
    "humann/**"

    script:
    """
    humann -i ${fastq} \
        -o 'humann' \
        --nucleotide-database /Users/sdavis2/bigdb/chocophlan/chocophlan/ \
        --taxonomic-profile ${metaphlan_bugs_list} \
        --protein-database /Users/sdavis2/bigdb/${params.uniref} \
        --metaphlan-options "--bowtie2db /Users/sdavis2/bigdb/metaphlan" \
        --threads ${task.cpus}
    """
}


workflow {
    fasterq_dump(Channel.from(params.runs))
    concat_fastq(fasterq_dump.out.fastq_files.collect())
    install_metaphlan_db()
    uniref_db()
    chocophlan_db()
    metaphlan_bugs_list(concat_fastq.out.fastq) 
    metaphlan_markers(metaphlan_bugs_list.out.metaphlan_bt2)
    humann(concat_fastq.out.fastq, metaphlan_bugs_list.out.metaphlan_bugs_list)
    // concat_fastq(fasterq_dump.out.collect()) | wc_fastq | view
}
