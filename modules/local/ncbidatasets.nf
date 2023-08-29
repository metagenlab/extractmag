process NCBI_DATASETS {
    tag "$genome_accession"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ncbi-datasets-pylib::14.27.0--pyhdfd78af_0" : null)
    container 'staphb/ncbi-datasets:15.2.0'

    input:
    val(genome_accession)

    output:
    tuple val(genome_accession), path("ncbi_dataset/*/*/*fna")  , emit: fna

    script:
    """
    datasets download --filename ncbi_dataset.zip genome accession --include genome ${genome_accession}
    unzip *zip
    """

}
