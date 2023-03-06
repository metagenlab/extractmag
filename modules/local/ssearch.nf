
process SSEARCH {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/fasta3:36.3.8h--hec16e2b_0"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path('*_ssearch.txt')      , emit: ssearch
    path "versions.yml"                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    ssearch36  -b 1 -m 10 -n -z 3 ${fasta} ${db} > ${prefix}_ssearch.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ssearch36: \$( ssearch36 | grep version | sed 's/version: //' )
    END_VERSIONS
    """

}

