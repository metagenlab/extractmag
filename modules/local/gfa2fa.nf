process GFA2FA {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*fasta")      , emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if file --mime-type -L "${gfa}" | grep -q gzip; then
        echo "${gfa} is gzipped"
        gunzip -c ${gfa} >  ${prefix}.tmp.fa
        gfa=${prefix}.tmp.fa
    else
        echo "${gfa} is not gzipped"
        gfa=${gfa}
    fi

    awk '/^S/{print ">"\$2"\\n"\$3}' \$gfa | fold > ${prefix}.fasta
    rm -f *tmp*
    """

}
