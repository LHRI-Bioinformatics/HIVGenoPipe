process BBTOOLS_SEAL {
    tag "$meta.id"
    label 'process_high'

    container 'staphb/bbtools:latest'

    publishDir "${params.outdir}", mode: 'copy'
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/filter_reads/*.{log,err}"
    echo false



    input:
    tuple val(meta), path(reads)
    path(reference_fasta)

    output:
    tuple val(meta), path("filter_reads/*.reference.{R1,R2}.fastq"), emit: filtered_reads
    file("filter_reads/*")
    tuple val(meta), file("filter_reads/*.sealstats.txt"), emit: sealstats
    file("logs/filter_reads/*.{log,err}")
    tuple val(meta), env(total_read_count), env(reference_read_count), env(reference_read_pct), emit: read_stats


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir -p filter_reads logs/filter_reads
        log_file=logs/filter_reads/${prefix}.log
        err_file=logs/filter_reads/${prefix}.err

        # time stamp + capturing tool versions
        date | tee -a \$log_file \$err_file > /dev/null
        # dedupe
        dedupe.sh in=${reads[0]} in2=${reads[1]} out=${prefix}.deduped.fastq.gz ac=f
        # reformat
        reformat.sh in=${prefix}.deduped.fastq.gz out1=${prefix}_R1_001.deduped.fastq.gz out2=${prefix}_R2_001.deduped.fastq.gz
        # seal command
        seal.sh in=${prefix}_R1_001.deduped.fastq.gz in2=${prefix}_R2_001.deduped.fastq.gz ref=${reference_fasta} out=filter_reads/${prefix}.reference.R1.fastq $args \\
        out2=filter_reads/${prefix}.reference.R2.fastq \\
        outu=filter_reads/${prefix}.other.fastq stats=filter_reads/${prefix}.sealstats.txt 2>> \$err_file >> \$log_file

        total_read_count=\$(grep "#Total" filter_reads/${prefix}.sealstats.txt | awk 'BEGIN{ FS=OFS="\t" }{ print \$2 }' )
        if [ -z "\$total_read_count" ] ; then total_read_count="NA" ; fi

        reference_read_count=\$(grep "#Matched" filter_reads/${prefix}.sealstats.txt | awk 'BEGIN{ FS=OFS="\t" }{ print \$2 }')
        if [ -z "\$reference_read_count" ] ; then reference_read_count="NA" ; fi

        reference_read_pct=\$(grep "#Matched" filter_reads/${prefix}.sealstats.txt | awk 'BEGIN{ FS=OFS="\t" }{ print \$3 }')
        if [ -z "\$reference_read_pct" ] ; then reference_read_pct="NA" ; fi

        """
}
