// custom script by Brad to collapse contigs from mafft output
process MAFFT_CONSENSUS {
    errorStrategy 'ignore'
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::biopython=1.71" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75' :
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mafft.hybrid_consensus.fasta"), emit: hybrid_fasta
    file("logs/mafft/*.{log,err}")


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p mafft logs/mafft
    log_file=logs/mafft/${prefix}.mafft_consensus.log
    err_file=logs/mafft/${prefix}.mafft_consensus.err

    consensusFromMAFFT.py \\
        ${fasta} ${prefix}_mafft_hybrid_consensus ${prefix}.mafft.hybrid_consensus.fasta 2>> \$err_file >> \$log_file \\

    """
}
