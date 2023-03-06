process GATBMINIA {
    tag "$meta.id"
    label 'process_high'

    container "docker://cimendes/gatb-minia-pipeline:31.07.2020-1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*final.contigs.fa')                            , emit: contigs
    path "*final_k.txt"                                                   , emit: final_k
    path "versions.yml"                                                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz

        gatb ${args} -1 ${prefix}_1.fastq.gz --nb-cores ${task.cpus} -o ${prefix} > gatb.log 2>&1

        readlink ${prefix}_final.contigs.fa | sed 's/.*_//' | tr -dc '0-9' > ${prefix}_final_k.txt 

        rm -r *.unitigs* *.h5 || true
        rm *list_reads* || true

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatbminia: ${task.container}
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        
        gatb ${args} -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz --nb-cores ${task.cpus} -o ${prefix} > gatb.log 2>&1

        link=\$(readlink ${prefix}_final.contigs.fa)
        echo \$link | sed 's/.*_//' | tr -dc '0-9' > final_k.txt
        rm ${prefix}_final.contigs.fa
        cp \$link ${prefix}_final.contigs.fa

        rm -rf *unitigs* *.h5 || true
        rm *list_reads* || true

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatbminia: ${task.container}
        END_VERSIONS
        """
    }
}
