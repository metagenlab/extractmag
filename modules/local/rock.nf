
process ROCK {
    tag "$meta.id"
    label 'process_medium'

    container "metagenlab/rock:2.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.norm.fastq.gz') , optional:true, emit: reads
    path "versions.yml"                                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        /opt/conda/bin/pigz -d --keep --stdout ${reads[0]} > ${prefix}_1.fastq
        echo ${prefix}_1.fastq > input.txt
        echo ${prefix}_1.norm.fastq > output.txt

        /opt/conda/bin/rock  -i input.txt -o output.txt ${args}
        /opt/conda/bin/pigz ${prefix}_1.norm.fastq

        rm *fastq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rock: \$( rock -h | grep ROCK | sed 's/Cop.*//' )
        END_VERSIONS
        """
    } else {
        """
        /opt/conda/bin/pigz -d --keep --stdout ${reads[0]} > ${prefix}_1.fastq
        /opt/conda/bin/pigz -d --keep --stdout ${reads[1]} > ${prefix}_2.fastq
        
        echo "${prefix}_1.fastq,${prefix}_2.fastq" > input.txt
        echo "${prefix}_1.norm.fastq,${prefix}_2.norm.fastq" > output.txt

        /opt/conda/bin/rock  -i input.txt -o output.txt ${args}
        /opt/conda/bin/pigz ${prefix}_1.norm.fastq ${prefix}_2.norm.fastq 

        rm *fastq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rock: \$( rock -h | grep ROCK | sed 's/Cop.*//' )
        END_VERSIONS
        """
    }
}
