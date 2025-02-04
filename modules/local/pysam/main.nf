process PYSAM{
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pysam=0.19.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' :
        'quay.io/biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: pysam_stats

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pysam_test.py ${bam} $args

    """

}
