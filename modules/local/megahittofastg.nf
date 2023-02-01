process MEGAHITFATOFASTG {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::megahit=1.2.9 conda-forge::pigz=2.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' :
        'quay.io/biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' }"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(k_contigs)

    output:
    tuple val(meta), path("*fastg")      , emit: fastg

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # SRX018017_SRR038425.contigs.fa.gz
    k=\$(zcat ${prefix}.contigs.fa.gz | head -n1 | sed 's/_.*//' | sed 's/>k//') || true
    echo \$k > k.txt
    zcat k\$k.contigs.fa.gz > k\$k.contigs.fa
    megahit_core contig2fastg \$k k\$k.contigs.fa > MEGAHIT_${prefix}_k\$k.fastg
    """

}
