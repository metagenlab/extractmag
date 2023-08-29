
process CONTIGS_STATS {
    tag "$meta.id"
    label 'process_low'

    container "metagenlab/biopython:1.79c"

    input:
    tuple val(meta), path(fasta), path(gfa), path(contigs2taxonomy), path(coverm)

    output:
    tuple val(meta), path("*contigs.tsv"), emit: contigs_stats_tsv
    tuple val(meta), path("*subgraphs_cols.csv"), emit: subgraphs_cols
    tuple val(meta), path("*subgraphs_summary.tsv"), emit: subgraphs_summary
    tuple val(meta), path("*target_cols.csv"), emit: target_cols
    tuple val(meta), path("*.gfa"), emit: gfa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    contigs_summary.py -f ${fasta} -c ${coverm} -t ${contigs2taxonomy} -g ${gfa} -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "python": \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
