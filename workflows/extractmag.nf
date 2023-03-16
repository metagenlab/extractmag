/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

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
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                       } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                      } from '../modules/nf-core/modules/multiqc/main'
include { FASTQC as TRIM_FASTQC        } from '../modules/nf-core/modules/fastqc/main'  addParams( options: [:] )
include { FASTP                        } from '../modules/nf-core/modules/fastp/main'  addParams( options: [:] )
include { QUAST                        } from '../modules/nf-core/modules/quast/main'  addParams( options: [:] )
include { CHECKM_LINEAGEWF             } from '../modules/nf-core/modules/checkm/lineagewf/main' addParams( options: [:] )
include { CUSTOM_SRATOOLSNCBISETTINGS  } from '../modules/nf-core/modules/custom/sratoolsncbisettings/main'
include { SRATOOLS_PREFETCH            } from '../modules/nf-core/modules/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP         } from '../modules/nf-core/modules/sratools/fasterqdump/main'
include { SRA_IDS_TO_RUNINFO           } from '../modules/local/sra_ids_to_runinfo'
include { SRA_RUNINFO_TO_FTP           } from '../modules/local/sra_runinfo_to_ftp'
include { BLOOCOO                      } from '../modules/local/bloocoo'
include { PIGZ                         } from '../modules/local/compress'
include { GATBMINIA                    } from '../modules/local/gatb_minia'
include { MEGAHIT                      } from '../modules/nf-core/modules/megahit/main'
include { SPADES                       } from '../modules/nf-core/modules/spades/main'
include { FATOGFA                      } from '../modules/local/miniatogfa'
include { MEGAHITFATOFASTG             } from '../modules/local/megahittofastg'
include { FASTGTOGFA                   } from '../modules/local/bandage_fastg2gfa'
include { TAXONOMY as MINIA_TAXONOMY } from '../subworkflows/local/taxonomy'
include { TAXONOMY as MEGAHIT_TAXONOMY } from '../subworkflows/local/taxonomy'
include { TAXONOMY as SPADES_TAXONOMY  } from '../subworkflows/local/taxonomy'
include { ROCK                         } from '../modules/local/rock'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow EXTRACTMAG {
    take:
    sra_ids  // channel: [ val(meta), val(id) ]

    main:

    ch_versions = Channel.empty()

    sra_ids
    .collect().view()

    //
    // MODULE: Get SRA run information for public database ids
    //
    SRA_IDS_TO_RUNINFO (
        sra_ids,
        params.ena_metadata_fields ?: ''
    )
    ch_versions = ch_versions.mix(SRA_IDS_TO_RUNINFO.out.versions.first())

    //
    // MODULE: Parse SRA run information, create file containing FTP links and read into workflow as [ meta, [reads] ]
    //
    SRA_RUNINFO_TO_FTP (
        SRA_IDS_TO_RUNINFO.out.tsv
    )
    ch_versions = ch_versions.mix(SRA_RUNINFO_TO_FTP.out.versions.first())

    SRA_RUNINFO_TO_FTP
        .out
        .tsv
        .splitCsv(header:true, sep:'\t')
        .map {
            meta ->
                meta.single_end = meta.single_end.toBoolean()
                [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
        }
        .unique()
        .map { meta, reads -> [ meta, meta.run_accession ] }
        .set { ch_sra_reads }

    

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    //
    // configure SRA
    //
    CUSTOM_SRATOOLSNCBISETTINGS()
    def settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings
    ch_versions = ch_versions.mix( CUSTOM_SRATOOLSNCBISETTINGS.out.versions )

    //
    // Download reads
    //

    SRATOOLS_PREFETCH ( ch_sra_reads, settings )
    ch_versions = ch_versions.mix( SRATOOLS_PREFETCH.out.versions.first())

    SRATOOLS_FASTERQDUMP ( SRATOOLS_PREFETCH.out.sra, settings )
    ch_versions = ch_versions.mix( SRATOOLS_FASTERQDUMP.out.versions.first() )

    //
    // FastQC raw reads
    //
    FASTQC (
        SRATOOLS_FASTERQDUMP.out.reads // channel: [ val(meta), [ reads ] ]
    )

    /*
     * Trim Reads
     */
    FASTP(SRATOOLS_FASTERQDUMP.out.reads, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions.first().ifEmpty(null))

    /*
     * Run FastQC on trimmed reads
     */
    TRIM_FASTQC(FASTP.out.reads)
    ch_versions = ch_versions.mix(TRIM_FASTQC.out.versions.first().ifEmpty(null))

    /*
     * READ CORRECTION
     */
    BLOOCOO(FASTP.out.reads)
    ch_versions = ch_versions.mix(BLOOCOO.out.versions.first().ifEmpty(null))

    /*
     * ASSEMBLY
     * - MINIA
     * - MEGAHIT
     * - SPADES
     */
    GATBMINIA(BLOOCOO.out.reads)
    ch_versions = ch_versions.mix(GATBMINIA.out.versions.first().ifEmpty(null))

    MEGAHIT(BLOOCOO.out.reads)
    ch_versions = ch_versions.mix(MEGAHIT.out.versions.first().ifEmpty(null))

    ROCK(BLOOCOO.out.reads)
    ROCK.out.reads.map { meta, reads -> [meta, reads, [], []] }
        .set { ch_spades }

    
    // SPADES(ch_spades, [])
    // ch_versions = ch_versions.mix(SPADES.out.versions.first().ifEmpty(null))

    /* convert fasta to gfa */
    FATOGFA(GATBMINIA.out.contigs, GATBMINIA.out.final_k)
    MEGAHITFATOFASTG(MEGAHIT.out.contigs, MEGAHIT.out.k_contigs)

    FASTGTOGFA(MEGAHITFATOFASTG.out.fastg)

    /* taxo annotation gfa */
    // gtdb taxnonomy
    // root;d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Streptomycetales;f__Streptomycetaceae;g__Streptomyces;s__Streptomyces aureoverticillatus
    MEGAHIT_TAXONOMY(FASTGTOGFA.out.gfa, "MEGAHIT_")
    MINIA_TAXONOMY(FATOGFA.out.gfa, "MINIA_")
    //SPADES_TAXONOMY(SPADES.out.gfa, "SPADES_")


    //
    // MODULE: Extract target genomes based on taxonomy and assembly graph data
    //


    //
    // MODULE: MultiQC
    //
    
    workflow_summary    = WorkflowExtractmag.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
	ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
	//ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1][0]}.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
