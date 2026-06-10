process REPORT_SAMTOOLS_DEPTH_AND_STATS{
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(samtools_stats), path(pysamstats_file)

    output:
    tuple val(meta), path("*.csv"), emit: samtools_depth_and_stats


    script: 
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_type = task.ext.prefix ?: "${meta.sample_type}"
    """
    cat ${samtools_stats} | grep ^SN | cut -f 2- > summary.stats
    samtools_stats_and_pysam_depth.py summary.stats ${pysamstats_file} ${prefix} ${sample_type}
    """

}
