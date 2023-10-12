
process CONTIGS_STATS {
    tag "$meta.id"
    label 'process_low'

    container "metagenlab/biopython:1.79c"

    input:
    tuple val(meta), path(fasta), path(gfa), path(contigs2taxonomy), path(coverm) 
    path(contigs_path)
    path(contigs_fa)

    output:
    tuple val(meta), path("*contigs.tsv"), emit: contigs_stats_tsv
    tuple val(meta), path("*subgraphs_cols.csv"), emit: subgraphs_cols
    tuple val(meta), path("*subgraphs_summary.tsv"), emit: subgraphs_summary
    tuple val(meta), path("*target_cols.csv"), emit: target_cols
    tuple val(meta), path("*chromosome.fasta"), emit: target_fasta
    tuple val(meta), path("*unconnected_contigs.fasta"), emit: unconnected_contigs
    tuple val(meta), path("*.gfa"), emit: gfa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def spades_paths = contigs_path.name != 'NO_FILE' ? "-p $contigs_path" : ''
    def spades_contigs = contigs_fa.name != 'NO_FILE' ? "-s $contigs_fa" : ''
    """
    contigs_summary.py -f ${fasta} -c ${coverm} -t ${contigs2taxonomy} -g ${gfa} -o ${prefix} ${spades_paths} ${spades_contigs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "python": \$(/opt/conda/bin/python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
