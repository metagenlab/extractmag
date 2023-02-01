process GFA2FA {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(gfa)
    val(folder)

    publishDir "results/gfa/${folder}"

    output:
    tuple val(meta), path("*fasta")      , emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '/^S/{print ">"\$2"\\n"\$3}' ${gfa} | fold > ${prefix}.fasta
    """

}
