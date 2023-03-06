

process BARRNAP {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/barrnap:0.9--3"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*barrnap_16S.fasta')      , emit: rrna_16s
    tuple val(meta), path('*barrnap_23S.fasta')      , emit: rrna_23s
    path "versions.yml"                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    barrnap --threads $task.cpus ${fasta} --outseq ${prefix}_barrnap.fasta
    awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' ${prefix}_barrnap.fasta > ${prefix}_barrnap_single.fasta
    grep ">16S" -A 1 ${prefix}_barrnap_single.fasta > ${prefix}_barrnap_16S.fasta
    grep ">23S" -A 1 ${prefix}_barrnap_single.fasta > ${prefix}_barrnap_23S.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barrnap: \$( barrnap --version | sed 's/barrnap //' )
    END_VERSIONS
    """

}

