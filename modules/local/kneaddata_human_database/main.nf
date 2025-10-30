process KNEADDATA_HUMAN_DATABASE {
    label 'process_single'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'
    storeDir "${params.store_dir}"

    output:
    path "human_genome", emit: kd_genome, type: "dir"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo ${PWD}
    mkdir -p human_genome
    kneaddata_database --download human_genome bowtie2 human_genome
    """

    stub:
    """
    mkdir -p human_genome
    touch human_genome/hg37dec_v0.1.1.bt2
    """
}
