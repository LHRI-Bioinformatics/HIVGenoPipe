process REPORT_STANFORD_SUMMARY{
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(json_files)

    output:
    tuple val(meta), path("*.tsv"), emit: dr_summaries

    script: // This script is bundled with the pipeline, in nf-core/hivgenopipe/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stanford_summary.py -j ${json_files}
    """
}
