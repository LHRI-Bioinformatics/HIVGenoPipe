process REPORT_LOG_REPORT{

    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"
    input:
    val ready
    path(execution_trace_log)

    output:
    path("*error.log"), emit: trinity_error_log

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${execution_trace_log} | grep FAILED > all.error.log

    """

}
