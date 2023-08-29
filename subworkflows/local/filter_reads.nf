

include { GET_UNMAPPED                } from '../../modules/local/get_unmapped.nf'

workflow FILTER_READS {
    
    take:
    reads
    
    
    main:

    ch_versions = Channel.empty()

    // map to reference and extract unmapped 
    GET_UNMAPPED(reads)

    emit:
    reads = GET_UNMAPPED.out.reads

} 
