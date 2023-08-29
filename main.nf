#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/extractmag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/extractmag
    Website: https://nf-co.re/extractmag
    Slack  : https://nfcore.slack.com/channels/extractmag
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



// Validate input parameters
WorkflowExtractmag.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from './subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


WorkflowMain.initialise(workflow, params, log)


// Check if --input file is empty
ch_input = file(params.input, checkIfExists: true)

if (ch_input.isEmpty()) {exit 1, "File provided with --input is empty: ${ch_input.getName()}!"}

// Read in ids from --input file
Channel
    .from(file(params.input, checkIfExists: true))
    .splitCsv(header:true, sep:'\t', strip:true)
    .filter{ it.R1 == null }
    .map { it.accession }
    .set { ch_sra }

ch_sra.view()

// Read in ids from --input file
Channel
    .from(file(params.input, checkIfExists: true))
    .splitCsv(header:true, sep:'\t', strip:true)
    .filter{ it.R1 == null }
    .filter{ it.filtering == "" }
    .map { create_sample_channel(it) }
    .set { ch_samples_no_filtering }

Channel
    .from(file(params.input, checkIfExists: true))
    .splitCsv(header:true, sep:'\t', strip:true)
    .filter{ it.R1 == null }
    .filter{ it.filtering != "" }
    .map { create_sample_channel(it) }
    .set { ch_samples_filtering }

ch_samples_no_filtering.view()
ch_samples_filtering.view()

/*
Channel
    .from(file(params.input, checkIfExists: true))
    .splitCsv(header:true, sep:'\t', strip:true)
    .filter{ it.R1 != null }
    .set { ch_local_fastq }

ch_local_fastq.view()
*/

Channel
    .from(file(params.input, checkIfExists: true))
    .splitCsv(header:true, sep:'\t', strip:true)
    .filter{ it.filtering != '' }
    .map { it.filtering }
    .unique()
    .set { ch_references }

// ch_references.collect().view()

def create_sample_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.accession =  row.accession
    meta.filtering =  row.filtering
    meta.assembly  =  row.assembly

    return meta
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACTMAG   } from './workflows/extractmag'
include { DOWNLOAD_SRA } from './subworkflows/local/download_sra'
include { DOWNLOAD_REF } from './subworkflows/local/download_ref'
include { FILTER_READS } from './subworkflows/local/filter_reads'
include { ASSEMBLY as ASSEMBLY_FILT     } from './subworkflows/local/assembly'
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS as DUMP_NO_FILT } from './subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools/main'
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS as DUMP_FILT } from './subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools/main'


def map_join(channel_a, channel_b, key){
    channel_a
        .map{ it -> [it[key], it] }
        .cross(channel_b.map{it -> [it[key], it]})
        .map { it[0][1] + it[1][1] }
}

//
// WORKFLOW: Run main nf-core/extractmag analysis pipeline
//
workflow NFCORE_EXTRACTMAG {


    // NO FILTERING
    DUMP_NO_FILT(
        ch_samples_no_filtering.map{ [ it, it.accession ] }, // [meta, accession]
        []
        )


    // FILTERING
    // download and index ref genomes
    DOWNLOAD_REF(ch_references)

    // download and filter reads
    DUMP_FILT(
        ch_samples_filtering.map{ [ it, it.accession ] }, // [meta, accession]
        []
    )

    // prepare channel for filtering
    DOWNLOAD_REF.out.indexed_ref
        .map{ [ "filtering":it[0], "ref_index":it[1]] }
        .set{ref}
    
    // DUMP_FILT.out.reads.view()

    map_join(ch_samples_filtering, ref, "filtering").map{ [ ["id": it.id, "accession":it.accession, "filtering":it.filtering, "assembly":it.assembly], it.ref_index] }
    .join(DUMP_FILT.out.reads)
    .set{filtering_input} //
    
    FILTER_READS(filtering_input)

    // [id:ERR2810282, accession:ERR2810282, filtering:GCA_000349085.1]
    // [GCA_000349085.1, /scratch/hdd3/tpillone/2022_07_extractmag/work/cc/bf439b9f1185c9c6271b6dc4ad676b/bwamem2]
    //ch_samples.join(DOWNLOAD_REF.out.indexed_ref, remainder: true, by: [2, 0]).view()


    //ch_samples.join(DUMP_NO_FILT.out.reads).set{extractmap_input}

    //extractmap_input.view()

    // combine sra and local channels
    // DOWNLOAD_SRA.out.reads.view()
    FILTER_READS.out.reads.view()
    // run main workflow with all samples
    EXTRACTMAG (FILTER_READS.out.reads.mix(DUMP_NO_FILT.out.reads))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_EXTRACTMAG ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
