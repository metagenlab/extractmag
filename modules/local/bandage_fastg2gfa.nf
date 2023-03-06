process FASTGTOGFA {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/biocontainers/bandage:0.8.1--hc9558a2_2"

    input:
    tuple val(meta), path(fastg)

    output:
    tuple val(meta), path("*gfa")      , emit: gfa

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Bandage reduce ${fastg} ${prefix}.gfa
    """

}
