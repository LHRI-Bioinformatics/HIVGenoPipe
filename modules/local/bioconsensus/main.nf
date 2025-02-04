process BIOCONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::biopython=1.71" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75' :
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(fasta)

    output:
    path("*_consensus.fasta"), emit: consensus_fasta

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bioalign_consensus.py \\
        ${fasta} > ${prefix}_consensus.fasta
    """
}
