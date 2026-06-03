/*
 * GTDB taxonomy conversion
 *
 * MetaPhlAn 4 reports taxonomy using its own SGB-based clade names. This step
 * translates a MetaPhlAn relative-abundance profile into a GTDB-taxonomy
 * profile using the official SGB->GTDB assignment table (downloaded and cached
 * by the sgb_to_gtdb_db process). The conversion itself is performed by the
 * vendored bin/sgb_to_gtdb_profile.py helper, which:
 *   - substitutes 1:1 SGB->GTDB mappings directly,
 *   - bins (sums) n:1 mappings into the shared GTDB taxon,
 *   - aggregates abundances up the GTDB lineage.
 *
 * Outputs are published under the same per-sample/branch directory layout used
 * by the rest of the pipeline, in a `gtdb` subdirectory.
 */

process metaphlan_to_gtdb {
    label 'profiling'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/${meta.branch ? meta.branch + '/' : ''}gtdb", pattern: "{*tsv.gz,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    cpus 1
    memory "2g"

    input:
    val meta
    path metaphlan_profile
    path sgb2gtdb_db

    output:
    val(meta), emit: meta
    path 'gtdb_profile.tsv', emit: gtdb_profile
    path 'gtdb_profile.tsv.gz', emit: gtdb_profile_gz
    path ".command*"
    path "versions.yml"

    stub:
    """
    touch gtdb_profile.tsv
    touch gtdb_profile.tsv.gz
    touch .command.run
    touch versions.yml
    """

    script:
    """
    mapping_file=\$(find ${sgb2gtdb_db} -name '*SGB2GTDB.tsv' | head -n 1)
    if [ -z "\$mapping_file" ]; then
        echo "ERROR: no SGB2GTDB mapping table found in ${sgb2gtdb_db}" >&2
        exit 1
    fi

    sgb_to_gtdb_profile.py \
        -i ${metaphlan_profile} \
        -o gtdb_profile.tsv \
        -m "\$mapping_file"

    gzip -c gtdb_profile.tsv > gtdb_profile.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    versions:
        metaphlan: \$( echo \$(metaphlan --version 2>&1 ) | awk '{print \$3}')
        sgb2gtdb_mapping: \$(basename \$mapping_file)
    END_VERSIONS
    """
}
