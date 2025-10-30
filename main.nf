#!/usr/bin/env nextflow
/*
========================================================================================
    seandavi/curatedmetagenomicsnextflow
========================================================================================
    Github : https://github.com/seandavi/curatedmetagenomicsnextflow
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

// Print help message if requested
if (params.help) {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
        nextflow run seandavi/curatedmetagenomicsnextflow --input samplesheet.tsv --outdir results
    
    Mandatory arguments:
        --input [file]                  Path to input samplesheet (TSV format)
        --outdir [path]                 Path to the output directory
    
    Optional arguments:
        --skip_humann                   Skip HUMAnN functional profiling [default: false]
        --local_input                   Use local FASTQ files instead of SRA download [default: false]
        --store_dir [path]              Directory to store reference databases [default: 'databases']
        --metaphlan_index [str]         MetaPhlAn index version [default: 'latest']
        --chocophlan [str]              ChocoPhlAn database version [default: 'full']
        --uniref [str]                  UniRef database version [default: 'uniref90_diamond']
        --organism_database [str]       Organism database for KneadData [default: 'human_genome']
    
    Profile options:
        -profile [str]                  Configuration profile to use. Available: test, docker, singularity, local, google, anvil, alpine, unitn
    """.stripIndent()
    exit 0
}

// Validate mandatory parameters
def input_file = params.input ?: params.metadata_tsv
if (!input_file && (!params.run_ids || !params.sample_id)) {
    error "ERROR: Please provide either --input/--metadata_tsv or both --run_ids and --sample_id"
}

// Set output directory
def outdir = params.outdir ?: params.publish_dir ?: 'results'
params.outdir = outdir

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { CURATEDMETAGENOMICSNEXTFLOW } from './workflows/curatedmetagenomicsnextflow'

//
// WORKFLOW: Run main workflow
//
workflow {
    CURATEDMETAGENOMICSNEXTFLOW ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
