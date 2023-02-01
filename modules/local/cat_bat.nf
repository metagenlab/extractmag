process CAT {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/biocontainers/cat:5.2.3--hdfd78af_1"

    input:
    tuple val(meta), path(fasta)
    path(dbpath)
    path(taxodb)
    val(exprefix)

    output:
    tuple val(meta), path('*.CAT.contig2classification.txt')      , emit: contig2classification
    path "versions.yml"                                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    if [[ ${fasta} == *gz ]]
    then
        gunzip -k ${fasta};
        fasta_file=\$(echo ${fasta} | sed 's/.gz//');
    else
        fasta_file=${fasta};
    fi

    CAT contigs --no_stars -c \$fasta_file -n ${task.cpus} -d ${dbpath} -t ${taxodb}
    mv out.CAT.contig2classification.txt ${exprefix}${prefix}.CAT.contig2classification.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$( CAT -v | sed 's/ (.*//' | sed 's/CAT //' )
    END_VERSIONS
    """

}