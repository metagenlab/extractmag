
process QUAST {

    container "registry.hub.docker.com/staphb/quast:latest"

    input:
        tuple val(meta), file(contigs)
        
    output:
        file "${meta.id}_assembly_metrics/*"
    
    script:
    """
        quast.py -o ${meta.id}_assembly_metrics ${contigs}
    """
}
