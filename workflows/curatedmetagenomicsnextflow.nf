/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FASTERQ_DUMP                  } from '../modules/local/fasterq_dump/main'
include { LOCAL_FASTQC                  } from '../modules/local/local_fastqc/main'
include { KNEADDATA                     } from '../modules/local/kneaddata/main'
include { INSTALL_METAPHLAN_DB          } from '../modules/local/install_metaphlan_db/main'
include { METAPHLAN_UNKNOWN_VIRUSES_LISTS } from '../modules/local/metaphlan_unknown_viruses_lists/main'
include { METAPHLAN_UNKNOWN_LIST        } from '../modules/local/metaphlan_unknown_list/main'
include { METAPHLAN_MARKERS             } from '../modules/local/metaphlan_markers/main'
include { SAMPLE_TO_MARKERS             } from '../modules/local/sample_to_markers/main'
include { CHOCOPHLAN_DB                 } from '../modules/local/chocophlan_db/main'
include { UNIREF_DB                     } from '../modules/local/uniref_db/main'
include { KNEADDATA_HUMAN_DATABASE      } from '../modules/local/kneaddata_human_database/main'
include { KNEADDATA_MOUSE_DATABASE      } from '../modules/local/kneaddata_mouse_database/main'
include { HUMANN                        } from '../modules/local/humann/main'

/*
========================================================================================
    HELPER FUNCTIONS
========================================================================================
*/

def generate_row_tuple(row) {
    accessions = row.NCBI_accession.split(';')
    sample_id = row.sample_id
    return [sample:sample_id, accessions: accessions, meta: row]
}

def generate_row_tuple_local(row) {
    fpaths = row.file_paths.split(';')
    sample_id = row.sample_id
    return [sample:sample_id, fpaths: fpaths, meta: row]
}

def generate_sample_metadata_single_sample(sample_id, run_ids) {
    accessions = run_ids.split(';')
    return [sample: sample_id, accessions: accessions, meta: null]
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CURATEDMETAGENOMICSNEXTFLOW {

    // Determine input source and create samples channel
    samples = null
    
    // Allow EITHER input/metadata_tsv or run_ids/sample_id
    def input_file = params.input ?: params.metadata_tsv
    
    if (input_file == null) {
        if (params.run_ids == null || params.sample_id == null) {
            error "Either input/metadata_tsv or run_ids/sample_id must be provided"
        } else {
            samples = Channel.value(
                generate_sample_metadata_single_sample(
                    params.sample_id, 
                    params.run_ids
                )
            )
            FASTERQ_DUMP(samples)
            fastq_meta = FASTERQ_DUMP.out.fastq
        }
    } else {
        if (params.local_input) {
            samples = Channel
                .fromPath(input_file)
                .splitCsv(header: true, quote: '"', sep:'\t')
                .map { row -> generate_row_tuple_local(row) }
            
            LOCAL_FASTQC(samples)
            fastq_meta = LOCAL_FASTQC.out.fastq
        } else {
            samples = Channel
                .fromPath(input_file)
                .splitCsv(header: true, quote: '"', sep:'\t')
                .map { row -> generate_row_tuple(row) }
            
            FASTERQ_DUMP(samples)
            fastq_meta = FASTERQ_DUMP.out.fastq
        }
    }

    // Download/install databases
    INSTALL_METAPHLAN_DB()
    UNIREF_DB()
    CHOCOPHLAN_DB()
    KNEADDATA_HUMAN_DATABASE()
    KNEADDATA_MOUSE_DATABASE()
    
    // Quality control with KneadData
    KNEADDATA(
        fastq_meta,
        KNEADDATA_HUMAN_DATABASE.out.kd_genome.collect(),
        KNEADDATA_MOUSE_DATABASE.out.kd_mouse.collect()
    )

    // MetaPhlAn taxonomic profiling
    METAPHLAN_UNKNOWN_VIRUSES_LISTS(
        KNEADDATA.out.fastq,
        INSTALL_METAPHLAN_DB.out.metaphlan_db.collect()
    )

    METAPHLAN_MARKERS(
        METAPHLAN_UNKNOWN_VIRUSES_LISTS.out.metaphlan_bt2,
        INSTALL_METAPHLAN_DB.out.metaphlan_db.collect()
    )

    SAMPLE_TO_MARKERS(
        METAPHLAN_UNKNOWN_VIRUSES_LISTS.out.metaphlan_sam,
        INSTALL_METAPHLAN_DB.out.metaphlan_db.collect()
    )

    // HUMAnN functional profiling (optional)
    if (!params.skip_humann) {
        HUMANN(
            KNEADDATA.out.fastq,
            METAPHLAN_UNKNOWN_VIRUSES_LISTS.out.metaphlan_unknown_list,
            CHOCOPHLAN_DB.out.chocophlan_db,
            UNIREF_DB.out.uniref_db
        )
    }
}
