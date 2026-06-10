process PYSAMSTATS_PARSER{
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"
    input:
    tuple val(meta), path(stats_tsv)

    output:
    tuple val(meta), path("*.tsv"), optional: true, emit: pysam_stats_results
    tuple val(meta), path("*parsed.tsv"),   optional: true, emit: parsed_tsv
    tuple val(meta), path("*Quality*.tsv"), optional: true, emit: quality_report
    tuple val(meta), path("*.fasta"), optional: true, emit: fastas
    tuple val(meta), path("*Amb*.fasta"), optional: true, emit: amb_fastas
    tuple val(meta), path("*consensus.fasta"), optional: true, emit: consensus_fasta
    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pysamstats_parse.py ${stats_tsv} ${prefix} $args

    """

}
