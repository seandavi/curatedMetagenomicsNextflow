#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process fasterq_dump {
    cpus 8
    memory "16g"
    errorStrategy 'retry'
    maxRetries 3    

    errorStrategy 'retry'
    maxRetries 3

    tag "$srr"

    input:
      tuple val(samp), val(srr)


    output:
    path "out.fastq", emit: fastq
    val(samp), emit: samp
    // path "*_fastqc/fastqc_data.txt", emit: fastqc_data

    script:
      """
      fasterq-dump \
          --skip-technical \
          --force \
          --threads ${task.cpus} \
          --split-files ${srr.join(' ')}
      cat *.fastq > out.fastq
      #fastqc --extract *.fastq
      """
}


process concat_fastq {
    cpus 1
    memory "1g"

    input:
      tuple samp, path(x)
    output:
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
    publishDir "${params.publish_dir}/${samp}/metaphlan"

    time "1d"
    cpus 32
    memory "32g"

//    input:
//    path metaphlan_db

    input:
    path fastq
    path metaphlan_db


    output:
    path 'bowtie2.out', emit: metaphlan_bt2
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

    """
}

process metaphlan_markers {
    publishDir "${params.publish_dir}/${samp}/metaphlan"

    cpus 4
    memory "16g"

    input:
    path metaphlan_bt2
    path metaphlan_db
    val samp

    output:
    path "marker_abundance.tsv", emit: marker_abundance
    path "marker_presence.tsv", emit: marker_presence

    script:
    """
    metaphlan --input_type bowtie2out \
        --index ${params.metaphlan_index} \
        --bowtie2db metaphlan \
        -t marker_pres_table \
        -o marker_presence.tsv \
        ${metaphlan_bt2}
    metaphlan --input_type bowtie2out \
        --index ${params.metaphlan_index} \
        --bowtie2db metaphlan \
        -t marker_ab_table \
        -o marker_abundance.tsv \
        ${metaphlan_bt2}
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
    publishDir "${params.publish_dir}/${samp}/humann"
    cpus 32
    time "7d"
    memory "32g"

    input:
    val samp
    path fastq
    path metaphlan_bugs_list // metaphlan_bugs_list.tsv
    path chocophlan_db
    path uniref_db

    output:
    path "out*.tsv"
    path "files.txt"

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
    find . > files.txt
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
    data = Channel.fromPath(params.csv_file)
       .splitCsv(header: true, sep: "\t")
       .map { tuple( it.sampleID, it.NCBI_accession.tokenize(';')) }
    // data | check | groupTuple | concat_fastq
    fasterq_dump(data) | transpose | groupTuple | concat_fastq
    install_metaphlan_db()
    uniref_db()
    chocophlan_db()
    metaphlan_bugs_list(
        fasterq_dump.out.samp,
        fasterq_dump.out.fastq,
        install_metaphlan_db.out.metaphlan_db.collect())
    metaphlan_markers(
        fasterq_dump.out.samp,
        metaphlan_bugs_list.out.metaphlan_bt2,
        install_metaphlan_db.out.metaphlan_db.collect())
    humann(
        fasterq_dump.out.samp,
        fasterq_dump.out.fastq,
        metaphlan_bugs_list.out.metaphlan_bugs_list,
        chocophlan_db.out.chocophlan_db,
        uniref_db.out.uniref_db)
}
