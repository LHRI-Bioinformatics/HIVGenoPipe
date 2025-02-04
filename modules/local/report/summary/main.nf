process REPORT_SUMMARY{
    // tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pandas=1.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    path grouped_stats_files

    output:
    path("*.tsv"), emit: final_report

    script: // This script is bundled with the pipeline, in nf-core/hivgenopipe/bin/
    def args = task.ext.args ?: ''
    // def sample_type = "${meta.sample_type}"
    def files_list       = grouped_stats_files.join(' ')
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """

    final_report_format.py ${grouped_stats_files}


    """

    // # cat ${interop_stats} >> interop_stats.txt
    // # for file in ${files_list}; do
    // if [ "${sample_type}" == 'control_positive' ];
    // then
    //     echo "control positive"
    //     cat ${files_list} >> positive_final_report.txt
    // elif [ "${sample_type}" == 'positive' ] || [ "${sample_type}" == 'negative' ] || ["${sample_type}" == 'test' ];
    // then
    //     echo "not control positive"
    //     cat ${files_list} >> final_report.txt
    // else
    //     echo "no valid type"
    // fi
    // # done
}
