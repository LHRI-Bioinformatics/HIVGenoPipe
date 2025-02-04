/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { BWA_INDEX                       } from '../../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM                         } from '../../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_INDEX                  } from '../../modules/nf-core/modules/samtools/index/main'
include { PYSAMSTATS                      } from '../../modules/local/pysamstats/pysamstats/main_edit'
include { PYSAMSTATS_PARSER               } from '../../modules/local/pysamstats/parser/main'
include { REPORT_CHECK_COVERAGE           } from '../../modules/local/report/check_coverage/main'
include { SAMTOOLS_DEPTH                  } from '../../modules/nf-core/modules/samtools/depth/main'
include { SAMTOOLS_STATS                  } from '../../modules/nf-core/modules/samtools/stats/main_edit'
include { REPORT_SAMTOOLS_DEPTH_AND_STATS } from '../../modules/local/report/samtools_depth_and_stats/main'
include { REPORT_NGS_STATS                } from '../../modules/local/report/ngs_stats/main_new_copy_control'

workflow CONTROL_TEST {
    take:
    ch_ref_fasta           // channel: [ fasta ] NL4.3
    ch_cleaned_reads       // channel: [ meta, reads ] from bbtools
    ch_read_stats          // channel: [ meta, trim_stats, seal_stats ]
    // add read stats from trimmomatic and bbtools as input for another ngs_stats module here
    // adjust output for this subwf to be ngs_stats output

    main:

    ch_versions = Channel.empty()

    BWA_INDEX(
        ch_ref_fasta
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    BWA_MEM(
        ch_cleaned_reads,
        BWA_INDEX.out.index,
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    SAMTOOLS_INDEX(
        BWA_MEM.out.bam
        )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())


  //Map the bai to the bam for use in pysamstats
    ch_pysamstats_input = Channel.empty()
    BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai -> [ meta, bam, bai ]
        }
        .set { ch_bam_bai }

    // ch_pysamstats_input = ch_bam_bai.concat(ch_ref_fasta)


    PYSAMSTATS(
        ch_bam_bai,
        ch_ref_fasta
    )


    PYSAMSTATS_PARSER(
        PYSAMSTATS.out.pysamstats_tsv
    )

    SAMTOOLS_STATS (
        ch_bam_bai
    )

    SAMTOOLS_DEPTH (
        BWA_MEM.out.bam
    )

    REPORT_CHECK_COVERAGE (
        PYSAMSTATS_PARSER.out.amb_fastas,
        ch_ref_fasta

    )

    ch_depth_and_stats = Channel.empty()
    ch_depth_and_stats = SAMTOOLS_STATS.out.stats.join(SAMTOOLS_DEPTH.out.tsv)

    REPORT_SAMTOOLS_DEPTH_AND_STATS (
        ch_depth_and_stats
    )



    ch_control_ngs_stats = ch_read_stats.join(REPORT_SAMTOOLS_DEPTH_AND_STATS.out.samtools_depth_and_stats)
    ch_control_ngs_stats = ch_control_ngs_stats.join(REPORT_CHECK_COVERAGE.out.contig_coverage_stats)

    REPORT_NGS_STATS (
        ch_control_ngs_stats
    )

    ch_control_output = REPORT_NGS_STATS.out.ngs_stats
        .map {
            meta, file ->
            new_meta = meta.clone()
            new_meta.sample_type = "control_positive"
            [new_meta, file]
        }


    emit:
    report            = ch_control_output         // channel: [ val(meta),  ngs_stats ]


    versions          = ch_versions                     // channel: [ versions.yml ]
}
