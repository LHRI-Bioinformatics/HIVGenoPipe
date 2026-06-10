process REPORT_GAGPOL_QUALITY{
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11.0   conda-forge::biopython=1.80 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'quay.io/biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    tuple val(meta), path(amb_jsons), path(parsed_pysamstats), path(final_consensus)

    output:
    tuple val(meta), path("*.csv"), emit: gagpol_quality
    tuple val(meta), path("*truncated.csv"), emit: truncated_pysamstats

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    json_gagpol_stats.py -j ${amb_jsons} -p ${parsed_pysamstats} -s${final_consensus} -o ${prefix}.gagpol_quality.csv  $args

    """

}
