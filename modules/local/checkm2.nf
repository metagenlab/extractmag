process CHECKM2 {
    tag "$meta.id"
    label 'process_high'
    maxForks 2

    container "quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0"

    input:
    tuple val(meta), path(fasta)
    
    output:
    path "*quality_report.tsv",       emit: tsv
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    checkm2 predict --threads ${task.cpus} --input *fasta --output-directory checkm --force --general --database_path ${params.checkm2_database} --remove_intermediates
    mv checkm/quality_report.tsv ${meta.id}.quality_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS

    """
}
