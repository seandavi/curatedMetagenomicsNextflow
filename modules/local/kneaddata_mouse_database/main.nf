process KNEADDATA_MOUSE_DATABASE {
    label 'process_single'

    container 'docker://seandavi/curatedmetagenomics:metaphlan4.1.0'
    storeDir "${params.store_dir}"

    output:
    path "mouse_C57BL", emit: kd_mouse, type: "dir"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo ${PWD}
    mkdir -p mouse_C57BL
    kneaddata_database --download mouse_C57BL bowtie2 mouse_C57BL
    """

    stub:
    """
    mkdir -p mouse_C57BL
    touch mouse_C57BL/mouse_C57BL_6NJ_Bowtie2_v0.1.bt2
    """
}
