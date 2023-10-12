

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                       } from '../modules/nf-core/fastqc/main'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'
include { FASTQC as TRIM_FASTQC        } from '../modules/nf-core/fastqc/main'  addParams( options: [:] )
include { FASTP                        } from '../modules/nf-core/fastp/main'  addParams( options: [:] )
include { QUAST                        } from '../modules/local/quast'  addParams( options: [:] )
include { CHECKM_LINEAGEWF             } from '../modules/nf-core/checkm/lineagewf/main' addParams( options: [:] )
include { BLOOCOO                      } from '../modules/local/bloocoo'
include { NGMERGE                      } from '../modules/local/ngmerge'
include { FLASH                        } from '../modules/nf-core/flash/main'
include { PIGZ                         } from '../modules/local/compress'
include { MOTUS_PROFILE                } from '../modules/nf-core/motus/profile/main'
include { COVERM                       } from '../modules/local/coverm'
include { CONTIGS_STATS                } from '../modules/local/contigs_stats'
include { ASSEMBLY                     } from '../subworkflows/local/assembly'
include { BANDAGE_IMAGE                } from '../modules/local/bandage_image'
include { CHECKM2                      } from '../modules/local/checkm2'
include { BANDAGE_IMAGE as BANDAGE_IMAGE_SUBG           } from '../modules/local/bandage_image'
include { GTDBTK_CLASSIFYWF } from '../modules/nf-core/gtdbtk/classifywf/main'

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
    reads  // channel: [ val(meta), [R1, R2]  ]

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
    // BLOOCOO(FASTP.out.reads)
    //ch_versions = ch_versions.mix(BLOOCOO.out.versions.first().ifEmpty(null))

    //
    // MODULE: Extract target genomes based on taxonomy and assembly graph data
    //

    // FLASH(BLOOCOO.out.reads)
    ASSEMBLY(FASTP.out.reads)
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions.first().ifEmpty(null))


    FASTP.out.reads.join(ASSEMBLY.out.contigs).set{reads_and_contigs}

    reads_and_contigs.view()

    // get depth stats
    COVERM(reads_and_contigs)
    QUAST(ASSEMBLY.out.contigs)
    
    // contigs stats 

    ASSEMBLY.out.contigs.join(ASSEMBLY.out.gfa).join(ASSEMBLY.out.taxonomy_contigs).join(COVERM.out.tsv).join(ASSEMBLY.out.contigs_paths).join(ASSEMBLY.out.contigs_spades)
    .multiMap{meta, contigs, gfa, taxo_contigs, coverm, paths, contigs_spades ->
    ch1: [meta, contigs, gfa, taxo_contigs, coverm]
    ch2: paths
    ch3: contigs_spades}.set{contigs_data}

    CONTIGS_STATS(contigs_data)
    CONTIGS_STATS.out.gfa.join(CONTIGS_STATS.out.target_cols).set{bandage_in}
    
    // plot assembly (sub)graph of target taxon only
    BANDAGE_IMAGE_SUBG(bandage_in)
    // ASSEMBLY.out.gfa.join(CONTIGS_STATS.out.subgraphs_cols).set{bandage_in2}
    // BANDAGE_IMAGE(bandage_in2)

    GTDBTK_CLASSIFYWF(CONTIGS_STATS.out.target_fasta, tuple("r214", params.gtdb), params.gtdb_mash)

    CHECKM2(CONTIGS_STATS.out.target_fasta)

    // plot GC vs Cov 
    // full graph

    // subset graph
    
    
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
