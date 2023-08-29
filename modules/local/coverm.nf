


process COVERM {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/coverm:0.6.1--h1535e20_5"

    input:
    tuple val(meta), path(reads), path(contigs)
    
    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz 

        coverm contig --output-format sparse --single ${prefix}.fastq.gz --reference ${contigs} --threads $task.cpus  -m mean covered_fraction covered_bases variance length count rpkm > ${prefix}_coverm.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz

        coverm contig --output-format sparse --coupled ${prefix}_1.fastq.gz ${prefix}_1.fastq.gz --reference ${contigs} --threads $task.cpus  -m mean covered_fraction covered_bases variance length count rpkm > ${prefix}_coverm.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            coverm: \$(coverm -V)
        END_VERSIONS
        """

    }
}
