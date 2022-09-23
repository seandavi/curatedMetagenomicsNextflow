#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process fasterq_dump {
    publishDir "${params.publish_dir}/${meta.sample}/fasterq_dump", pattern: "{fastq_line_count.txt,*_fastqc/fastqc_data.txt,sampleinfo.txt,.command*}"
    
    // maxForks 80
    cpus 8
    memory "16g"

    tag "${meta.sample}"

    input:
    val meta

    output:
    val(meta), emit: meta
    path "out.fastq.gz", emit: fastq
    path "*_fastqc/fastqc_data.txt", emit: fastqc_data
    path "fastq_line_count.txt"
    path ".command*"
    path "sampleinfo.txt"
    path "fastq-scan.json"
    path "versions.yml"

    script:
    """

    echo "accessions: ${meta.accessions}" > sampleinfo.txt
    echo "starting fasterq-dump"
    for accession in ${meta.accessions.join(" ")}; do
        echo "downloading \$accession"
        aws s3 cp s3://sra-pub-run-odp/sra/\$accession/\$accession . \
            --no-sign-request
        mv \$accession \$accession.sra
        fasterq-dump --threads ${task.cpus} \
            --skip-technical \
            --force \
            --split-files \$accession.sra 
    done
    ls -ld
    echo "fasterq-dump done"
    wc -l *.fastq > fastq_line_count.txt

    echo "combining fastq files and gzipping"
    cat *.fastq | pv | pigz -p ${task.cpus} > out.fastq.gz

    echo "sampling fastq file for fastqc"
    seqtk sample -s100 out.fastq.gz 50000 > out_sample.fastq

    echo "running fastqc"
    fastqc --extract out_sample.fastq

    echo "running fastq-scan"
    pigz -p ${task.cpus} -dc out.fastq.gz | fastq-scan > fastq-scan.json


    echo "collecting version info"
    cat <<-END_VERSIONS > versions.yml
    versions:
        awscli: \$(aws --version)
        fastq-scan: \$( echo \$(fastq-scan -v 2>&1) | sed 's/^.*fastq-scan //' )
        pigz: \$( echo \$(pigz --version 2>&1 ) | sed 's/^.*pigz //' )
        fastqc: \$( echo \$(fastqc --version 2>&1 ) | sed 's/^.*FastQC //' )
        fasterq-dump: \$( echo \$(fasterq-dump --version 2>&1 ) | head -2 | tail -1 | awk '{print \$3}')
    END_VERSIONS

    echo "finalizing fasterqc-dump" 

    """
}



process install_metaphlan_db {
    cpus 4
    memory '8g'

    storeDir "${params.store_dir}"

    output:
    path 'metaphlan', emit: metaphlan_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    script:
    """
    metaphlan --install --index latest --bowtie2db metaphlan

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS

    """
}

process metaphlan_bugs_list {
    publishDir "${params.publish_dir}/${meta.sample}/metaphlan_bugs_list", pattern: "{*tsv.gz,.command*}"

    tag "${meta.sample}"
    
    cpus 16
    memory { 16.GB * task.attempt }
    
    input:
    val meta
    path fastq
    path metaphlan_db


    output:
    val(meta), emit: meta
    path 'bowtie2.out.gz', emit: metaphlan_bt2
    path 'metaphlan_bugs_list.tsv', emit: metaphlan_bugs_list
    path 'metaphlan_bugs_list.tsv.gz', emit: metaphlan_bugs_list_gz
    path ".command*"
    path "versions.yml"

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

    gzip -c metaphlan_bugs_list.tsv > metaphlan_bugs_list.tsv.gz
    gzip bowtie2.out

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}

process metaphlan_markers {
    publishDir "${params.publish_dir}/${meta.sample}/metaphlan_markers/"
    
    tag "${meta.sample}"

    cpus 2
    memory "8g"

    input:
    val meta
    path metaphlan_bt2
    path metaphlan_db

    output:
    val meta, emit: meta
    path "marker_abundance.tsv.gz", emit: marker_abundance
    path "marker_presence.tsv.gz", emit: marker_presence
    path ".command*"
    path "versions.yml"

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

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        bowtie2: \$( echo \$(bowtie2 --version 2>&1 ) | awk '{print \$3}')
    END_VERSIONS
    """
}


process chocophlan_db {
    cpus 1
    memory "1g"

    storeDir "${params.store_dir}"

    output:
    path "chocophlan", emit: chocophlan_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    script:
    """
    humann_databases --update-config no --download chocophlan ${params.chocophlan} .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}


process uniref_db {
    cpus 1
    memory "1g"

    storeDir "${params.store_dir}"

    output:
    path "uniref", emit: uniref_db, type: 'dir'
    path ".command*"
    path "versions.yml"

    script:
    """
    humann_databases --update-config no --download uniref ${params.uniref} .

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}


process humann {
    publishDir "${params.publish_dir}/${meta.sample}/humann"
    cpus 32
    memory { 32.GB + 16.GB * task.attempt }

    tag "${meta.sample}"

    input:
    val meta
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
    path ".command*"
    path "versions.yml"

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

    cat <<-END_VERSIONS > versions.yml
    versions:
        humann: \$( echo \$(humann --version 2>&1 ) | awk '{print \$2}')
    END_VERSIONS
    """
}

def generate_row_tuple(row) {
    accessions=row.NCBI_accession.split(';');
    study_id = row.study_name;
    sample_id = row.sample_id;
    sample_encoded = "${study_id}::${sample_id}".bytes.encodeBase64().toString()
    // Create a hash of sampleID and joined accessions for
    // use as a unique id.
    // rowhash = "${accessions.sort().join(' ')}".md5().toString()
    return [sample: sample_encoded, accessions: accessions, meta: row]
}


workflow {
    samples = Channel
        .fromPath(params.metadata_tsv)
        .splitCsv(header: true, quote: '"', sep:'\t')
        .map { row -> generate_row_tuple(row) }
    // for debugging: 
    // samples.view()

    fasterq_dump(samples)

    install_metaphlan_db()
    uniref_db()
    chocophlan_db()
    
    metaphlan_bugs_list(
        fasterq_dump.out.meta,
        fasterq_dump.out.fastq,
        install_metaphlan_db.out.metaphlan_db.collect())
    metaphlan_markers(
        metaphlan_bugs_list.out.meta,
        metaphlan_bugs_list.out.metaphlan_bt2,
        install_metaphlan_db.out.metaphlan_db.collect())
    humann(
        fasterq_dump.out.meta,
        fasterq_dump.out.fastq,
        metaphlan_bugs_list.out.metaphlan_bugs_list,
        chocophlan_db.out.chocophlan_db,
        uniref_db.out.uniref_db)
}
