/*
 * GTDB taxonomy conversion
 *
 * MetaPhlAn 4 reports taxonomy using its own SGB-based clade names. This step
 * translates a MetaPhlAn relative-abundance profile into a GTDB-taxonomy
 * profile using the official SGB->GTDB assignment table (downloaded and cached
 * by the sgb_to_gtdb_db process). The conversion itself is performed by the
 * vendored bin/cmgd_sgb_to_gtdb.py helper, which:
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
    # sgb_to_gtdb_db downloads the table to a known name (the basename of
    # params.sgb2gtdb_url), so reference it directly — no need to search the
    # staged dir (which is a symlink `find` won't descend anyway). -s also
    # catches an empty/failed download.
    mapping_file="${sgb2gtdb_db}/${file(params.sgb2gtdb_url).name}"
    if [ ! -s "\$mapping_file" ]; then
        echo "ERROR: SGB2GTDB mapping table missing or empty: \$mapping_file" >&2
        exit 1
    fi

    # Uniquely named (not 'sgb_to_gtdb_profile.py') so the container's baked
    # /usr/local/bin copy — an older version without -m — cannot shadow the
    # pipeline's vendored bin/ script on PATH.
    cmgd_sgb_to_gtdb.py \
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
