process REPORT_SIMILARITY_MATRIX{
    label 'process_low'

    conda "conda-forge::python=3.8.3 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    path(similarity_stats)

    output:
    path("*.png"), emit: similarity_matrix_png

    script: // This script is bundled with the pipeline, in nf-core/hivgenopipe/bin/
    def args = task.ext.args ?: ''
    """
    similarity_matrix.py -i ${similarity_stats}

    """

}
