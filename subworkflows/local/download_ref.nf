
include { BWAMEM2_INDEX                } from '../../modules/nf-core/bwamem2/index/main.nf'
include { NCBI_DATASETS                } from '../../modules/local/ncbidatasets.nf'


workflow DOWNLOAD_REF {
    
    take:
    references
    
    main:

    ch_versions = Channel.empty()

    // download ref 
    NCBI_DATASETS(references)

    // index
    BWAMEM2_INDEX(NCBI_DATASETS.out.fna)

    emit:
    indexed_ref = BWAMEM2_INDEX.out.index

} 
