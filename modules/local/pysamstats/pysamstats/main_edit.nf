process PYSAMSTATS{
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pysamstats=1.1.2-10 bioconda::samtools=1.15.1" : null)
    container "quay.io/biocontainers/pysamstats:1.1.2--py38h7cf9df2_10"

    input:
    tuple val(meta), path(bam), path(index)
    path(hybrid_fasta)
    // tuple val(meta2), path(hybrid_fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: pysamstats_tsv

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pysamstats --type variation_strand -d $args -D 1000000 --output ${prefix}_counts.tsv -f ${hybrid_fasta} ${bam}

    """

}
