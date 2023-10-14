

include { SRA_IDS_TO_RUNINFO           } from '../../modules/local/sra_ids_to_runinfo'
include { SRA_RUNINFO_TO_FTP           } from '../../modules/local/sra_runinfo_to_ftp'
include { CUSTOM_SRATOOLSNCBISETTINGS  } from '../../modules/nf-core/custom/sratoolsncbisettings/main'
include { SRATOOLS_PREFETCH            } from '../../modules/nf-core/sratools/prefetch/main'
include { SRATOOLS_FASTERQDUMP         } from '../../modules/nf-core/sratools/fasterqdump/main'
include { clean_work_dirs as clean_sra } from '../../modules/local/clean_work.nf'

workflow DOWNLOAD_SRA {
    
    take:
    sra_ids
    
    main:

    ch_versions = Channel.empty()


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
    // configure SRA
    //
    CUSTOM_SRATOOLSNCBISETTINGS()
    def settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings
    ch_versions = ch_versions.mix( CUSTOM_SRATOOLSNCBISETTINGS.out.versions )

    //
    // Download reads
    //

    // SRATOOLS_PREFETCH ( ch_sra_reads, settings )
    ch_versions = ch_versions.mix( SRATOOLS_PREFETCH.out.versions.first())

    // SRATOOLS_FASTERQDUMP ( SRATOOLS_PREFETCH.out.sra, settings )
    ch_versions = ch_versions.mix( SRATOOLS_FASTERQDUMP.out.versions.first() )

    // SRATOOLS_PREFETCH.out.sra.join(SRATOOLS_FASTERQDUMP.out.reads).set{sra_to_clean}
    // SRATOOLS_PREFETCH.out.sra.view()
    // SRATOOLS_FASTERQDUMP.out.reads.view()
    //sra_to_clean.view()

    if( params.delete_intermediates ) {                                                               
        clean_sra(sra_to_clean) 
        clean =  clean_sra.out                                                                        
    }
    else {
        clean =  null  
    }

    
    emit:
    reads = SRATOOLS_FASTERQDUMP.out.reads
    clean = clean
} 
