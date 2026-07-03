/*
 * Input preprocessing processes
 *
 * These processes normalize either remote SRA accessions or local FASTQ files
 * into the same output contract: a single gzipped FASTQ plus lightweight QC
 * artifacts used downstream.
 */

process fasterq_dump {
    label 'qc'
    label 'io_heavy'
    label 'download_retry'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/fasterq_dump", pattern: "{fastq_line_count.txt,*_fastqc/fastqc_data.txt,sampleinfo.txt,.command*}", mode: "${params.publish_mode}"

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
    path "versions.yml", emit: versions

    stub:
    """
    touch out.fastq.gz
    touch fastq_line_count.txt
    touch sampleinfo.txt
    touch versions.yml
    touch .command.run
    mkdir -p ${meta.sample}_fastqc
    touch ${meta.sample}_fastqc/fastqc_data.txt
    """

    script:
    """

    echo "accessions: ${meta.accessions}" > sampleinfo.txt
    # Read acquisition: ENA-first (published fastq), SRA fallback (curl .sra + fasterq-dump).
    # ENA serves usable fastq even for runs whose SRA object lacks a QUALITY column
    # (those make fasterq-dump exit 3); trying ENA first removes that whole failure class.
    # Runs ENA does not serve (ingestion lag, submitted-only, controlled access) fall
    # back to the original SRA path, so coverage is unchanged. See ADR-0012.
    fetch_ena() {
        acc=\$1
        report=\$(curl -sf --max-time 60 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=\${acc}&result=read_run&fields=fastq_ftp,fastq_md5") || return 1
        row=\$(printf '%s\\n' "\$report" | tail -n +2 | head -n1)
        urls=\$(printf '%s' "\$row" | cut -f1)
        md5s=\$(printf '%s' "\$row" | cut -f2)
        [ -z "\$urls" ] && return 1
        i=1
        for u in \$(printf '%s' "\$urls" | tr ';' ' '); do
            fn=\$(basename "\$u")
            echo "  ENA \$fn"
            curl -sf --max-time 3600 -o "\$fn" "https://\$u" || return 1
            m=\$(printf '%s' "\$md5s" | cut -d';' -f\$i)
            if [ -n "\$m" ]; then
                echo "\$m  \$fn" | md5sum -c - || return 1
            fi
            i=\$((i + 1))
        done
        return 0
    }

    fetch_sra() {
        acc=\$1
        echo "  SRA \$acc"
        curl -sf --max-time 3600 -o \$acc.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/\$acc/\$acc || return 1
        fasterq-dump --threads ${task.cpus} --skip-technical --force --split-files \$acc.sra
    }

    echo "starting read acquisition (ENA-first, SRA fallback)"
    for accession in ${meta.accessions.join(" ")}; do
        echo "acquiring \$accession"
        if ! fetch_ena \$accession; then
            echo "  ENA unavailable for \$accession; falling back to SRA"
            fetch_sra \$accession || { echo "ERROR: ENA and SRA both failed for \$accession"; exit 3; }
        fi
    done

    # ENA delivers *.fastq.gz; fasterq-dump delivers *.fastq. Unify to *.fastq.
    for f in *.fastq.gz; do [ -e "\$f" ] && gunzip -f "\$f"; done

    ls -ld
    echo "read acquisition done"
    wc -l *.fastq > fastq_line_count.txt

    echo "combining fastq files and gzipping"
    cat *.fastq | pv | pigz -p ${task.cpus} > out.fastq.gz

    echo "running fastqc"
    fastqc --extract out.fastq.gz


    echo "collecting version info"
    cat <<-END_VERSIONS > versions.yml
    versions:
        awscli: \$(aws --version)
        fastqc: \$( echo \$(fastqc --version 2>&1 ) | sed 's/^.*FastQC //' )
        fasterq-dump: \$( echo \$(fasterq-dump --version 2>&1 ) | head -2 | tail -1 | awk '{print \$3}')
    END_VERSIONS

    echo "finalizing fasterqc-dump"

    """
}

process local_fastqc {
    label 'qc'
    label 'io_heavy'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/local_fastqc", pattern: "{fastq_line_count.txt,*_fastqc/fastqc_data.txt,sampleinfo.txt,.command*}", mode: "${params.publish_mode}"

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
    path "versions.yml", emit: versions

    stub:
    """
    touch out.fastq.gz
    touch fastq_line_count.txt
    touch sampleinfo.txt
    touch versions.yml
    touch .command.run
    mkdir -p ${meta.sample}_fastqc
    touch ${meta.sample}_fastqc/fastqc_data.txt
    """

    script:
    """

    echo "file_paths: ${meta.fpaths}" > sampleinfo.txt
    echo "copying fastq files"
    for fpath in ${meta.fpaths.join(" ")}; do
        echo "copying \$fpath"
        ls -l \$fpath
        cp \$fpath .
    done
    ls -ld
    echo "copying done"
    gunzip *.gz
    wc -l *.fastq > fastq_line_count.txt

    echo "combining fastq files and gzipping"
    cat *.fastq | pv | pigz -p ${task.cpus} > out.fastq.gz

    echo "running fastqc"
    fastqc --extract out.fastq.gz


    echo "collecting version info"
    cat <<-END_VERSIONS > versions.yml
    versions:
        fastqc: \$( echo \$(fastqc --version 2>&1 ) | sed 's/^.*FastQC //' )
    END_VERSIONS
    """
}
