process BBTOOLS_REFORMAT {
    publishDir "${params.outdir}", mode: 'copy'
    tag "$meta.id"
    label 'process_high'
    cpus 4
    container 'staphb/bbtools:latest'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("reformat_fasta/*.fasta"), emit: reformat_fasta
    file("logs/reformat_fasta/*.reformat_fasta.{log,err}")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir -p logs/reformat_fasta reformat_fasta
        log_file=logs/reformat_fasta/${prefix}.reformat_fasta.log
        err_file=logs/reformat_fasta/${prefix}.reformat_fasta.err

        date | tee -a \$log_file \$err_file >> /dev/null

        reformat.sh in=${fastq} out=reformat_fasta/${prefix}.fasta

        """
}
