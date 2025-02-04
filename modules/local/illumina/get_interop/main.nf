process ILLUMINA_GET_INTEROP{
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::illumina-interop=:1.2.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/illumina-interop:1.2.4--hdbdd923_2' :
        'quay.io/biocontainers/illumina-interop:1.2.4--hdbdd923_2' }"

    input:
    path(rundir)

    output:
    path "interop.csv"        , emit: interop_csv

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    """
    interop_summary --csv=1 ${rundir} > interop.csv

    """

}
