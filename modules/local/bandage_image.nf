process BANDAGE_IMAGE {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::bandage_ng=2022.09' : null)
    container 'https://depot.galaxyproject.org/singularity/bandage_ng:2022.09'

    input:
    tuple val(meta), path(gfa), path(color_file)

    output:
    tuple val(meta), path('*.png'), emit: png
    tuple val(meta), path('*.svg'), emit: svg
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for gfa in `ls *gfa`;do 
        Bandage image \$gfa \${gfa/gfa/png} $args --colors ${color_file} --nodewidth 5 --iter 4 --depwidth 1
        Bandage image \$gfa \${gfa/gfa/svg} $args --colors ${color_file} --nodewidth 5 --iter 4 --depwidth 1
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bandage: \$(echo \$(Bandage --version 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
