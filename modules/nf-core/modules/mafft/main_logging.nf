process MAFFT {
    errorStrategy 'ignore'
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mafft=7.508"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.508--hec16e2b_0':
        'quay.io/biocontainers/mafft:7.508--hec16e2b_0' }"

    input:
    tuple val(meta), path(fasta)
    path(reference)

    output:
    tuple val(meta), path("*.fas"), emit: fas, optional: true
    file("mafft/logs/*.err")
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref = reference ? "$reference" : ''
    """
    mkdir -p mafft mafft/logs
    err_file=mafft/logs/${prefix}.mafft.err

    # time stamp + capturing tool versions
    date | tee -a \$err_file > /dev/null

    mafft \\
        --thread ${task.cpus} \\
        ${args} \\
        --add ${fasta} \\
        ${ref} \\
        2>> \$err_file > ${prefix}.fas \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}