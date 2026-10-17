/*
 * BAQLaVa viral profiling (ADR-0020)
 *
 * Full-depth branch only, gated by !skip_baqlava (which requires HUMAnN and a
 * bundle with a baqlava_container). The input is the reads the bundle's
 * `humann` task left unaligned against the sample's own ChocoPhlAn pangenomes,
 * given to BAQLaVa with --humann-depleted-fasta. That is the same depleted read
 * set BAQLaVa's own depletion step would make (HUMAnN with
 * --bypass-translated-search on the same profile), without a second HUMAnN
 * nucleotide pass, and against the real ChocoPhlAn: BAQLaVa's own step calls
 * HUMAnN without --nucleotide-database, so in the image it would deplete
 * against HUMAnN's bundled DEMO database.
 *
 * BAQLaVa's own HUMAnN calls take HUMAnN's default utility database, which
 * in the image is utility_DEMO, and HUMAnN refuses a non-demo input against
 * it. So the bundle's real utility mapping database is appended to both
 * passthrough parameter sets, after BAQLaVa's own defaults (configs/baqlava.cfg
 * at the pinned commit), which passing the option replaces.
 *
 * The input is linked as <sample>_bowtie2_unaligned.fa because BAQLaVa names
 * its outputs after the input file (minus that suffix), so the published file
 * is <sample>/baqlava/<sample>_BAQLaVa_profile.txt. The intermediates
 * (_tempfile_*, logs) are not published.
 *
 * Container and resources are in the process body, as for humann.
 */

process baqlava {
    container { params.humann_bundles[params.humann_bundle].baqlava_container }

    label 'functional_profile'

    publishDir "${params.publish_dir ?: "${params.publish_base_dir}/${workflow.manifest.name}/${workflow.manifest.version}"}/${meta.sample}/baqlava", pattern: "{*_BAQLaVa_profile.txt,.command*}", mode: "${params.publish_mode}"

    tag "${meta.sample}"

    // ponytail: first guess (bowtie2 on a 1.6 GB index + DIAMOND on 0.3 GB,
    // over already-depleted reads); re-size from the Bridges-2 test traces.
    cpus 16
    memory { 24.GB * task.attempt }
    time { 6.h * task.attempt }

    input:
    tuple val(meta), path(unaligned)
    path baqlava_db
    path utility_mapping_db

    output:
    val(meta), emit: meta
    path "${meta.sample}_BAQLaVa_profile.txt", emit: profile
    path ".command*"
    path "versions.yml", emit: versions

    stub:
    """
    touch ${meta.sample}_BAQLaVa_profile.txt
    touch .command.run
    cat <<-END_VERSIONS > versions.yml
    versions:
        baqlava: stub
        baqlava_image: stub
    END_VERSIONS
    """

    script:
    def image = params.humann_bundles[params.humann_bundle].baqlava_container
    """
    nucdb=\$(dirname "\$(find -L ${baqlava_db} -name '*.1.bt2' | head -1)")
    protdb=\$(dirname "\$(find -L ${baqlava_db} -name '*.dmnd' | head -1)")

    ln -s ${unaligned} ${meta.sample}_bowtie2_unaligned.fa
    baqlava -i ${meta.sample}_bowtie2_unaligned.fa \\
        --humann-depleted-fasta \\
        -o . \\
        --nucdb "\$nucdb" \\
        --protdb "\$protdb" \\
        --threads ${task.cpus} \\
        --humann-passthrough-parameters-nucleotide='--bowtie-options="--very-sensitive" --nucleotide-identity-threshold 0.0 --nucleotide-query-coverage-threshold 90.0 --utility-database ${utility_mapping_db}' \\
        --humann-passthrough-parameters-translated='--diamond-options="--top 1 --outfmt 6" --evalue 1.0 --translated-alignment diamond --translated-identity-threshold 80.0 --translated-query-coverage-threshold 90.0 --utility-database ${utility_mapping_db}'

    # BAQLaVa reports its config version (1.2.0-devel); the image tag pins
    # the build that actually ran (#117).
    cat <<-END_VERSIONS > versions.yml
    versions:
        baqlava: \$( baqlava --version 2>&1 | awk '{print \$2; exit}' )
        baqlava_image: ${image.replaceFirst('^docker://', '')}
    END_VERSIONS
    """
}
