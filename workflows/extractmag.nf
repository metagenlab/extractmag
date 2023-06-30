

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                       } from '../modules/nf-core/fastqc/main'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'
include { FASTQC as TRIM_FASTQC        } from '../modules/nf-core/fastqc/main'  addParams( options: [:] )
include { FASTP                        } from '../modules/nf-core/fastp/main'  addParams( options: [:] )
include { QUAST                        } from '../modules/nf-core/quast/main'  addParams( options: [:] )
include { CHECKM_LINEAGEWF             } from '../modules/nf-core/checkm/lineagewf/main' addParams( options: [:] )
include { BLOOCOO                      } from '../modules/local/bloocoo'
include { PIGZ                         } from '../modules/local/compress'
include { GATBMINIA                    } from '../modules/local/gatb_minia'
include { MEGAHIT                      } from '../modules/nf-core/megahit/main'
include { SPADES                       } from '../modules/nf-core/spades/main'
include { FATOGFA                      } from '../modules/local/miniatogfa'
include { MEGAHITFATOFASTG             } from '../modules/local/megahittofastg'
include { FASTGTOGFA                   } from '../modules/local/bandage_fastg2gfa'
include { TAXONOMY as MINIA_TAXONOMY } from '../subworkflows/local/taxonomy'
include { TAXONOMY as MEGAHIT_TAXONOMY } from '../subworkflows/local/taxonomy'
include { TAXONOMY as SPADES_TAXONOMY  } from '../subworkflows/local/taxonomy'
include { ROCK                         } from '../modules/local/rock'
include { MOTUS_PROFILE                } from '../modules/nf-core/motus/profile/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Info required for completion email and summary
def multiqc_report = []
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)


workflow EXTRACTMAG {
    take:
    reads  // channel: [ val(meta), val(id) ]

    main:

    ch_versions = Channel.empty()


    //
    // FastQC raw reads
    //
    FASTQC (
        reads // channel: [ val(meta), [ reads ] ]
    )

    /*
     * Trim Reads
     */
    FASTP(reads, false, false)
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

    // motus taxnonomy profile

    MOTUS_PROFILE(FASTP.out.reads, params.motus_db)

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

    // ROCK(BLOOCOO.out.reads)
    // ROCK.out.reads.map { meta, reads -> [meta, reads, [], []] }
    //    .set { ch_spades }

    
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
