process REPORT_NGS_STATS{
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::pysam=0.16.0.1 conda-forge::biopython=1.78' :  null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0"
    }

    input:
    tuple val(meta), path(trim_stats), path(sealstats), path(samtools_depth_and_stats), path(contig_coverage)
    // , , path(self_align_bam), path(pysamstats_fastas), path(samtools_depth_and_stats)

    output:
    tuple val(meta), path("*.txt"), emit: ngs_stats
    //add to confirm completion of process so trace logs can show errors
    val true                               , emit: stats_complete

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ngs_stats.py -t ${trim_stats} -b ${sealstats} \\
    -d ${samtools_depth_and_stats} -c ${contig_coverage} -o control_subwf


    """
   // -s ${self_align_bam} -l 100 -p ${pysamstats_fastas} -d ${samtools_depth_and_stats}

}
