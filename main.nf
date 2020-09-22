#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//env_command = "env"
//res = env_command.execute().text

//println res


process fasterq_dump {
    publishDir "${params.publish_dir}/${workflow.sessionId}/${rowhash}/fasterq_dump", pattern: "{fastq_line_count.txt,*_fastqc/fastqc_data.txt,sampleinfo.txt}"
    
    cpus 4
    memory "16g"
    errorStrategy 'retry'
    maxRetries 4

    tag "${rowhash}"

    input:
    tuple val(samp), val(srr), val(rowhash)

    output:
    val(rowhash)
    path "out.fastq.gz", emit: fastq
    path "*_fastqc/fastqc_data.txt", emit: fastqc_data
    path "fastq_line_count.txt"


    script:
      """
      echo "sample: ${samp}\naccessions: ${srr}\nrowhash: ${rowhash}" > sampleinfo.txt
      fasterq-dump \
          --skip-technical \
          --force \
          --threads ${task.cpus} \
          --split-files ${srr.join(" ")}
      cat *.fastq | gzip > out.fastq.gz
      gunzip -c out.fastq.gz | wc -l > fastq_line_count.txt
      rm *.fastq
      seqtk sample -s100 out.fastq.gz 50000 > out_sample.fastq
      fastqc --extract out_sample.fastq
      rm out_sample.fastq
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
    publishDir "${params.publish_dir}/${workflow.sessionId}/${rowhash}/metaphlan_bugs_list", pattern: "*tsv.gz"

    tag "${rowhash}"

    time "1d"
    cpus 16
    memory "32g"


    input:
    val rowhash
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

    # gzip metaphlan_bugs_list.tsv
    gzip bowtie2.out
    """
}

process metaphlan_markers {
    publishDir "${params.publish_dir}/${workflow.sessionId}/${rowhash}/metaphlan_markers"
    
    tag "${rowhash}"
    
    cpus 1
    memory "16g"

    input:
    val rowhash
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
    publishDir "${params.publish_dir}/${workflow.sessionId}/${rowhash}/humann"
    cpus 16

    tag "${rowhash}"

    time "2d"
    memory "64g"

    input:
    val rowhash
    path fastq
    path metaphlan_bugs_list // metaphlan_bugs_list.tsv
    path chocophlan_db
    path uniref_db

    output:
    // lots of files....
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

def generate_row_tuple(row) {
    sampleID=row.sampleID;
    accessions=row.NCBI_accession.split(';');
    // Create a hash of sampleID and joined accessions for
    // use as a unique id.
    rowhash = "${sampleID} ${accessions.sort().join(' ')}".md5().toString()
    return tuple(sampleID, accessions, rowhash)
}

workflow {
    // Channel.from(samp).combine(Channel.from(runs)) | fasterq_dump | groupTuple | map { it -> [ it[0], it[1].flatten ] } | view
    samples = Channel.fromPath(params.metadata_tsv)
        .splitCsv(header: true, sep: "\t")
        .map { row -> generate_row_tuple(row) }    

    fasterq_dump(samples)

    install_metaphlan_db()
    uniref_db()
    chocophlan_db()
    
    metaphlan_bugs_list(
        fasterq_dump.out[0],
        fasterq_dump.out.fastq,
        install_metaphlan_db.out.metaphlan_db.collect())
    metaphlan_markers(
        fasterq_dump.out[0],
        metaphlan_bugs_list.out.metaphlan_bt2,
        install_metaphlan_db.out.metaphlan_db.collect())
    humann(
        fasterq_dump.out[0],
        fasterq_dump.out.fastq,
        metaphlan_bugs_list.out.metaphlan_bugs_list,
        chocophlan_db.out.chocophlan_db,
        uniref_db.out.uniref_db)
}
