process FATOGFA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(contigs)
    path final_k

    output:
    tuple val(meta), path("*gfa")      , emit: gfa

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    k=\$(cat ${final_k})
    file=${contigs}
    convertToGFA.py ${contigs} ${prefix}.gfa \$k
    """

}
