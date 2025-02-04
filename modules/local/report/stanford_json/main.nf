process REPORT_STANFORD_JSON{
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::sierrapy==0.4.3--pyh7cba7a3_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sierrapy:0.4.2--pyh7cba7a3_0' :
        'quay.io/biocontainers/sierrapy:0.4.3--pyh7cba7a3_0'}"

    input:
    tuple val(meta), path(fasta_files)

    output:
    tuple val(meta), path("*.json"), emit: dr_json

    script: // This script is bundled with the pipeline, in nf-core/hivgenopipe/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stanford_json.py -f ${fasta_files}
    """
}
