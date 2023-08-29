
include { GFA2FA           } from '../../modules/local/gfa2fa'
include { BARRNAP          } from '../../modules/local/barrnap'
include { CAT              } from '../../modules/local/cat_bat'
include { SSEARCH              } from '../../modules/local/ssearch'
include { PARSE_SSEARCH    } from '../../modules/local/parse_ssearch'

workflow TAXONOMY {
    take:
    gfa
    prefix

    main:

    GFA2FA(gfa) 

    CAT(GFA2FA.out.fasta, params.cat_db, params.cat_tax)

    BARRNAP(GFA2FA.out.fasta)
    SSEARCH(BARRNAP.out.rrna_16s, params.db_16S)
    PARSE_SSEARCH(SSEARCH.out.ssearch)

    emit:
    rrna_taxonomy = PARSE_SSEARCH.out.ssearch_tsv 
    contig2classification = CAT.out.contig2classification
    contigs = GFA2FA.out.fasta

}
