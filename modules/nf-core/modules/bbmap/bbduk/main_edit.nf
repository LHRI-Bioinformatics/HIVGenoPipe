process BBMAP_BBDUK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0':
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path('*.reference.fa'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    tuple val(meta), path('*.other.fa'), optional:true, emit: removed_reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        in=${reads} \\
        outm=${prefix}.assembly.reference.fa \\
        outu=${prefix}.assembly.other.fa \\
        threads=$task.cpus \\
        $args \\
        ref=${reference} \\
        &> ${prefix}.bbduk.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}