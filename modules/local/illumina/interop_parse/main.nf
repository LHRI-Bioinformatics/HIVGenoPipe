process ILLUMINA_INTEROP_PARSE{
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    path(interop_csv)

    output:
    path "interop_stats.csv"        , emit: interop_stats_csv

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
   
    """
    interop_parse.py --i ${interop_csv}

    """

}
