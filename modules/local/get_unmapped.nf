
process GET_UNMAPPED {
    publishDir "reference_genomes"
    cpus 20
    container 'docker://metagenlab/seqtk-samtools-bwa:latest'

    input:
        tuple val(meta), path(reference_genome_index), path(reads)

    output:
        path "unmapped_reads.lst", emit: unmapped_lst
        tuple val(meta), path("*_filt.fastq.gz"), emit: reads

    script:
    R1_unmapped = reads[0].name.replace(".fastq", "_filt.fastq")
    R2_unmapped = reads[1].name.replace(".fastq", "_filt.fastq")

    """
        # 4 unmapped 
        # 12 both reads unmapped
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
        bwa-mem2 mem -t ${task.cpus} \$INDEX ${reads} | samtools view -f 12 \
            | cut -f 1 | sort | uniq > unmapped_reads.lst
        seqtk subseq ${reads[0]} unmapped_reads.lst | pigz > ${R1_unmapped}
        seqtk subseq ${reads[1]} unmapped_reads.lst | pigz > ${R2_unmapped}
    """
}
