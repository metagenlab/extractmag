
include { GATBMINIA                    } from '../../modules/local/gatb_minia'
include { MEGAHIT                      } from '../../modules/nf-core/megahit/main'
include { SPADES                       } from '../../modules/nf-core/spades/main'
include { FATOGFA                      } from '../../modules/local/miniatogfa'
include { MEGAHITFATOFASTG             } from '../../modules/local/megahittofastg'
include { FASTGTOGFA                   } from '../../modules/local/bandage_fastg2gfa'
include { TAXONOMY as MINIA_TAXONOMY   } from '../../subworkflows/local/taxonomy'
include { TAXONOMY as MEGAHIT_TAXONOMY } from '../../subworkflows/local/taxonomy'
include { TAXONOMY as SPADES_TAXONOMY  } from '../../subworkflows/local/taxonomy'
include { ROCK                         } from '../../modules/local/rock'


/*
* ASSEMBLY
* - MINIA
* - MEGAHIT
* - SPADES
*/



workflow ASSEMBLY {
    
    take:
    reads
    
    main:

    ch_versions = Channel.empty()

    if ( params.assembly == 'spades' ) {

        // motus taxnonomy profile
        // MOTUS_PROFILE(FASTP.out.reads, params.motus_db)
        
        //ROCK(BLOOCOO.out.reads)
        // ROCK.out.reads.map { meta, reads -> [meta, reads, [], []] }
        //    .set { ch_spades }

        reads.map { meta, reads -> [meta, reads, [], []] }
        .set { ch_spades }

        SPADES(ch_spades, [])
        ch_versions = ch_versions.mix(SPADES.out.versions.first().ifEmpty(null))
        
        SPADES_TAXONOMY(SPADES.out.gfa, "SPADES_")
    }

    else if ( params.assembly == 'minia' ) {

        GATBMINIA(reads)
        ch_versions = ch_versions.mix(GATBMINIA.out.versions.first().ifEmpty(null))
        /* convert fasta to gfa */
        FATOGFA(GATBMINIA.out.contigs, GATBMINIA.out.final_k)
        MINIA_TAXONOMY(FATOGFA.out.gfa,"MINIA_")

    }
    else if ( params.assembly == 'megahit' ) {

        MEGAHIT(reads)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first().ifEmpty(null))
        MEGAHITFATOFASTG(MEGAHIT.out.contigs, MEGAHIT.out.k_contigs)

        FASTGTOGFA(MEGAHITFATOFASTG.out.fastg)
        /* taxo annotation gfa */
        // gtdb taxnonomy
        // root;d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Streptomycetales;f__Streptomycetaceae;g__Streptomyces;s__Streptomyces aureoverticillatus
        MEGAHIT_TAXONOMY(FASTGTOGFA.out.gfa, "MEGAHIT_")

    }
    else 
    {

        println "Unknown assembler"
    }


    if ( params.assembly == 'spades' ) {
        SPADES_TAXONOMY.out.contigs.set{contigs}
        taxonomy_contigs = SPADES_TAXONOMY.out.contig2classification
        taxonomy_rrna = SPADES_TAXONOMY.out.rrna_taxonomy
        gfa = SPADES.out.gfa
        contigs_paths = SPADES.out.contigs_paths
        contigs_spades = SPADES.out.contigs
     }
    if ( params.assembly == 'megahit' ) {
        MEGAHIT_TAXONOMY.out.contigs.set{contigs}
        taxonomy_contigs = MEGAHIT_TAXONOMY.out.contig2classification
        taxonomy_rrna = MEGAHIT_TAXONOMY.out.rrna_taxonomy
        gfa = FASTGTOGFA.out.gfa
        contigs_paths = null
        contigs_spades = null 
      }
    
    if ( params.assembly == 'minia' ) {
        MINIA_TAXONOMY.out.contigs.set{contigs}
        taxonomy_contigs = MINIA_TAXONOMY.out.contig2classification
        taxonomy_rrna = MINIA_TAXONOMY.out.rrna_taxonomy
        gfa = FATOGFA.out.gfa
        contigs_paths = null
        contigs_spades = null
    }



    emit:

    versions = ch_versions
    contigs = contigs
    taxonomy_contigs 
    taxonomy_rrna
    gfa
    contigs_paths
    contigs_spades
    


} 
