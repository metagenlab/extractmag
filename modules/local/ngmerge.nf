process NGMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ngmerge:15.2.0::0.3--ha92aebf_1" : null)
    // quay.io/biocontainers/bandage:0.8.1--hc9558a2_2
    container 'quay.io/biocontainers/ngmerge:0.3--ha92aebf_1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.merge.fastq.gz') , optional:true, emit: reads_merged
    path "versions.yml"                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz

    NGmerge -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -o ${prefix}.merge.fastq.gz -n ${task.cpus} -f ${prefix}.nomerge

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NGmerge: \$( NGmerge --version )
    END_VERSIONS
        """

}
