
process BLOOCOO {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bloocoo=1.0.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bloocoo:1.0.7--0' :
        'quay.io/biocontainers/bloocoo:1.0.7--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.corrected.fastq') , optional:true, emit: reads
    path "versions.yml"                                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        Bloocoo $args -nb-cores $task.cpus -file ${prefix}_1.fastq.gz
        mv ${prefix}_1.fastq_corrected.fastq ${prefix}_1.corrected.fastq 
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bloocoo: \$( Bloocoo -version | grep Bloocoo | sed 's/\\*//g' )
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        Bloocoo $args -nb-cores $task.cpus -file ${prefix}_1.fastq.gz,${prefix}_2.fastq.gz
        mv ${prefix}_1.fastq_corrected.fastq ${prefix}_1.corrected.fastq 
        mv ${prefix}_2.fastq_corrected.fastq ${prefix}_2.corrected.fastq 
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bloocoo: \$( Bloocoo -version | grep Bloocoo | sed 's/\\*//g' )
        END_VERSIONS
        """
    }
}
