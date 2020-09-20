#!/usr/bin/env nextflow

nextflow.preview.dsl=2

env_command = "env"
res = env_command.execute().text

println res


process fasterq_dump {
    cpus 4
    memory "16g"
    errorStrategy 'retry'
    maxRetries 4    

    tag "${srr}"

    input:
      tuple val(samp), val(srr)


    output:
      tuple val(samp), path("*.fastq")
    // path "*_fastqc/fastqc_data.txt", emit: fastqc_data

    script:
      """
      fasterq-dump \
          --skip-technical \
          --force \
          --threads ${task.cpus} \
          --split-files ${srr}
      #cat *.fastq > out.fastq
      #fastqc --extract *.fastq
      """
}


process concat_fastq {
    cpus 1
    memory "1g"

    tag "${params.samp}"

    input:
      tuple samp, path(x)
    output:
      val samp, emit: samp 
      path 'out.fastq', emit: fastq
      path 'wordcount.fastq'
    script:
      """
      ls -lah .
      echo ${x}
      cat ${x} >> out.fastq
      wc out.fastq > wordcount.fastq
      """
}


process install_metaphlan_db {
    cpus 1
    memory '32g'

    storeDir "${params.store_dir}"

    output:
    path 'metaphlan', emit: metaphlan_db, type: 'dir'

    script:
      """
      metaphlan --install --index latest --bowtie2db metaphlan
      """
}

process metaphlan_bugs_list {
    publishDir "${params.publish_dir}/${workflow.sessionId}/metaphlan"

    tag "${params.samp}"

    time "1d"
    cpus 16
    memory "32g"

//    input:
//    path metaphlan_db

    input:
    val samp
    path fastq
    path metaphlan_db


    output:
    path 'bowtie2.out.gz', emit: metaphlan_bt2
    path 'metaphlan_bugs_list.tsv', emit: metaphlan_bugs_list

    script:
    """
    find .
    metaphlan --input_type fastq \
        --index ${params.metaphlan_index} \
        --bowtie2db metaphlan \
        --samout sam.bz2 \
        --bowtie2out bowtie2.out \
        --nproc ${task.cpus} \
        -o metaphlan_bugs_list.tsv \
        ${fastq}
    gzip bowtie2.out
    """
}

process metaphlan_markers {
    publishDir "${params.publish_dir}/${workflow.sessionId}/metaphlan"
    
    tag "${params.samp}"
    
    cpus 1
    memory "16g"

    input:
    val samp
    path metaphlan_bt2
    path metaphlan_db

    output:
    path "marker_abundance.tsv.gz", emit: marker_abundance
    path "marker_presence.tsv.gz", emit: marker_presence

    script:
    """
    metaphlan --input_type bowtie2out \
        --index ${params.metaphlan_index} \
        --bowtie2db metaphlan \
        -t marker_pres_table \
        -o marker_presence.tsv \
        <( gunzip -c ${metaphlan_bt2} )    
    metaphlan --input_type bowtie2out \
        --index ${params.metaphlan_index} \
        --bowtie2db metaphlan \
        -t marker_ab_table \
        -o marker_abundance.tsv \
        <( gunzip -c ${metaphlan_bt2} )
    gzip *.tsv
    """
}


process chocophlan_db {
    cpus 1
    memory "1g"
    time "1d"

    storeDir "${params.store_dir}"

    output:
    path "chocophlan", emit: chocophlan_db, type: 'dir'

    script:
    """
    humann_databases --update-config no --download chocophlan ${params.chocophlan} .
    """
}


process uniref_db {
    cpus 1
    memory "1g"
    time "1d"

    storeDir "${params.store_dir}"

    output:
    path "uniref", emit: uniref_db, type: 'dir'

    script:
    """
    humann_databases --update-config no --download uniref ${params.uniref} .
    """
}


process humann {
    publishDir "${params.publish_dir}/${workflow.sessionId}/humann"
    cpus 16

    tag "${params.samp}"

    time "2d"
    memory "64g"

    input:
    val samp
    path fastq
    path metaphlan_bugs_list // metaphlan_bugs_list.tsv
    path chocophlan_db
    path uniref_db

    output:
    path "out_*.tsv.gz"

    script:
    """
    humann -i ${fastq} \
        -o '.' \
        --verbose \
        --metaphlan-options "-t rel_ab --index latest" \
        --nucleotide-database ${chocophlan_db} \
        --taxonomic-profile ${metaphlan_bugs_list} \
        --protein-database ${uniref_db} \
        --threads ${task.cpus}

    humann_renorm_table \
        --input out_pathabundance.tsv \
        --output out_cpm_pathabundance.tsv \
        --units cpm

    humann_renorm_table \
        --input out_genefamilies.tsv \
        --output out_cpm_genefamilies.tsv \
        --units cpm

    humann_renorm_table \
        --input out_genefamilies.tsv \
        --output out_relab_genefamilies.tsv \
        --units relab

    humann_renorm_table \
        --input out_pathabundance.tsv \
        --output out_relab_pathabundance.tsv \
        --units relab

    gzip out_*tsv
    """
}

process check {
    input:
      tuple val(samp), val(srr)
    output:
      tuple( val(samp), path('*.txt'))

    script:
    """
    echo ${samp} >> ${srr}_samp.txt
    echo ${srr} >> ${srr}_samp.txt
    """
}



workflow {
    runs = params.runs.split(";")
    samp = params.samp
    study = params.study
    // Channel.from(samp).combine(Channel.from(runs)) | fasterq_dump | groupTuple | map { it -> [ it[0], it[1].flatten ] } | view
    Channel.from(samp).combine(Channel.from(runs)) | fasterq_dump | groupTuple | map { it -> [ it[0], it[1].flatten() ] } | concat_fastq 
    install_metaphlan_db()
    uniref_db()
    chocophlan_db()
    metaphlan_bugs_list(
        concat_fastq.out.samp,
        concat_fastq.out.fastq,
        install_metaphlan_db.out.metaphlan_db.collect())
    metaphlan_markers(
        concat_fastq.out.samp,
        metaphlan_bugs_list.out.metaphlan_bt2,
        install_metaphlan_db.out.metaphlan_db.collect())
    humann(
        concat_fastq.out.samp,
        concat_fastq.out.fastq,
        metaphlan_bugs_list.out.metaphlan_bugs_list,
        chocophlan_db.out.chocophlan_db,
        uniref_db.out.uniref_db)
}
