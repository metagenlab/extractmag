process PARSE_SSEARCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::biopython=1.79" : null)
    container "metagenlab/biopython:1.79"

    input:
    tuple val(meta), path(ssearch_out)

    output:
    tuple val(meta), path("*_ssearch.tsv")      , emit: ssearch_tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ssearch_summary.py -i ${ssearch_out} -o ${prefix}_ssearch.tsv
    """

}
