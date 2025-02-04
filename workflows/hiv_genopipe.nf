/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTestpipeline.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check run dir

ch_rundir = params.rundir ? Channel.fromPath(params.rundir) : Channel.empty()

// Check if control ref fasta

ch_positive_control_ref_fasta = params.positive_control_ref ? Channel.fromPath(params.positive_control_ref) : Channel.empty()

// Check if metadata provided

ch_metadata = params.metadata ? Channel.fromPath(params.metadata) : Channel.empty()
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { INTEROP_STATS } from '../subworkflows/local/interop_stats'
include { CONTROL_TEST  } from '../subworkflows/local/control_test'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                      } from '../modules/nf-core/modules/fastqc/main'
include { TRIMMOMATIC                                 } from '../modules/nf-core/modules/trimmomatic/main'
include { BBTOOLS_SEAL                                } from '../modules/local/bbtools/seal/main'
include { TRINITY                                     } from '../modules/nf-core/modules/trinity/main'
include { BBMAP_BBDUK                                 } from '../modules/nf-core/modules/bbmap/bbduk/main_edit'
include { MAFFT                                       } from '../modules/nf-core/modules/mafft/main_logging'
include { MAFFT_CONSENSUS                             } from '../modules/local/mafft/consensus/main'
include { BWA_INDEX as BWA_INDEX_1                    } from '../modules/nf-core/modules/bwa/index/main_edit'
include { BWA_INDEX as BWA_INDEX_2                    } from '../modules/nf-core/modules/bwa/index/main_edit'
include { BWA_INDEX as BWA_INDEX_REF                  } from '../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM as BWA_MEM_1                        } from '../modules/nf-core/modules/bwa/mem/main_edit'
include { BWA_MEM as BWA_MEM_2                        } from '../modules/nf-core/modules/bwa/mem/main_edit'
include { BWA_MEM as BWA_MEM_REF                      } from '../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_INDEX                              } from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_STATS                              } from '../modules/nf-core/modules/samtools/stats/main_edit'
include { PYSAMSTATS                                  } from '../modules/local/pysamstats/pysamstats/main'
include { PYSAMSTATS_PARSER                           } from '../modules/local/pysamstats/parser/main'
include { REPORT_CHECK_COVERAGE                       } from '../modules/local/report/check_coverage/main'
include { SAMTOOLS_CONSENSUS                          } from '../modules/nf-core/modules/samtools/consensus/main'
include { SAMTOOLS_DEPTH as SAMTOOLS_DEPTH_1          } from '../modules/nf-core/modules/samtools/depth/main'
include { SAMTOOLS_DEPTH as SAMTOOLS_DEPTH_REF        } from '../modules/nf-core/modules/samtools/depth/main'
include { REPORT_SAMTOOLS_DEPTH_AND_STATS             } from '../modules/local/report/samtools_depth_and_stats/main'
include { REPORT_NGS_STATS as REPORT_NGS_STATS_ALL    } from '../modules/local/report/ngs_stats/main_new_copy'
include { REPORT_NGS_STATS as REPORT_NGS_STATS_PASS   } from '../modules/local/report/ngs_stats/main_new'
include { REPORT_STANFORD_JSON                        } from '../modules/local/report/stanford_json/main'
include { REPORT_STANFORD_SUMMARY                     } from '../modules/local/report/stanford_summary/main'
include { REPORT_SUMMARY                              } from '../modules/local/report/summary/main'
include { REPORT_PLOT_DEPTH as REPORT_PLOT_DEPTH_1    } from '../modules/local/report/plot_depth/main'
include { REPORT_PLOT_DEPTH as REPORT_PLOT_DEPTH_REF  } from '../modules/local/report/plot_depth/main'
include { REPORT_SIMILARITY                           } from '../modules/local/report/similarity_score/main'
include { REPORT_SIMILARITY_MATRIX                    } from '../modules/local/report/similarity_score/matrix/main'
include { DATABASE_IMPORT                             } from '../modules/local/database_import.nf'
include { MULTIQC                                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Info required for completion email and summary
def multiqc_report = []

workflow HIVGENOPIPE {
// input stuff
    channel
    .value( params.fasta )
    .map { filename -> file(filename, checkIfExists: true) }
    // .map { row -> [ [id:'ch_fasta'], row] }
    .set { ch_fasta }
    ch_versions = Channel.empty()

    // ch_fasta = Channel.fromPath("/mnt/LHRI/Members/ldotrang/HIVGenoPipe/HIVGenoPipe_Dev/bin/resources/*.fasta") //map { row -> [ [id:'ch_fasta'], row] }
    // ch_index = Channel.fromPath( "/mnt/LHRI/Members/ldotrang/HIVGenoPipe/HIVGenoPipe_Dev/bin/resources/*" )
    // ch_index.view()

    //
    //
    //IF THE RUNDIR CHANNEL IS NOT NULL
    //SUBWORKFLOW: Read illumina dir and generate run level stats
    if (params.rundir != null)
        {
            INTEROP_STATS(ch_rundir)
        }
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files

    INPUT_CHECK (
        ch_input
    )
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    TRIMMOMATIC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())


    BBTOOLS_SEAL (
        TRIMMOMATIC.out.trimmed_reads,
        ch_fasta
    )
    //Branch positive controls here after bbtools
    ch_control_samples = TRIMMOMATIC.out.trimmed_reads
        .branch { meta, reads ->
            positives : meta.sample_type == "positive"
                return [meta, reads]
            negatives : meta.sample_type == "negative"
                return [meta, reads]
        }
    ch_control_samples_stats = ch_ngs_stats = TRIMMOMATIC.out.summary.join(BBTOOLS_SEAL.out.sealstats)
    ch_control_samples_stats = ch_control_samples_stats.branch { meta, trim_stats, seal_stats ->
            positives : meta.sample_type == "positive"
                return [meta, trim_stats, seal_stats]
            negatives : meta.sample_type == "negative"
                return [meta, trim_stats, seal_stats]
        }

    if (params.positive_control_ref)
    {
        channel
        .value( params.positive_control_ref )
        .map { filename -> file(filename, checkIfExists: true) }
        .set { ch_positive_control_ref }

        CONTROL_TEST(
            ch_positive_control_ref,
            ch_control_samples.positives,
            ch_control_samples_stats.positives
            )
        }

    //Collect stats for all reads before we get failures in Trinity
    ch_ngs_stats_all = Channel.empty()
    ch_ngs_stats_all = TRIMMOMATIC.out.summary.join(BBTOOLS_SEAL.out.sealstats)
    REPORT_NGS_STATS_ALL (
        ch_ngs_stats_all
    )

    TRINITY (
        TRIMMOMATIC.out.trimmed_reads
    )
    ch_versions = ch_versions.mix(TRINITY.out.versions.first())

    BBMAP_BBDUK (
        TRINITY.out.transcript_fasta,
        ch_fasta
    )

    MAFFT (
        BBMAP_BBDUK.out.reads,
        ch_fasta

    )
    ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    MAFFT_CONSENSUS (
        MAFFT.out.fas
    )

    BWA_INDEX_1 (
        MAFFT_CONSENSUS.out.hybrid_fasta
    )
    ch_bwa_mem_input = TRIMMOMATIC.out.trimmed_reads.join(BWA_INDEX_1.out.index)

    BWA_MEM_1 (
        ch_bwa_mem_input,
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM_1.out.versions.first())

    SAMTOOLS_CONSENSUS (
        BWA_MEM_1.out.bam
    )

    BWA_INDEX_2 (
        SAMTOOLS_CONSENSUS.out.consensus
    )
    ch_bwa_mem_input_2 = TRIMMOMATIC.out.trimmed_reads.join(BWA_INDEX_2.out.index)

    BWA_MEM_2 (
        ch_bwa_mem_input_2,
        true
    )

    SAMTOOLS_INDEX (
        BWA_MEM_2.out.bam
    )
  //Map the bai to the bam for use in pysamstats
    BWA_MEM_2.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai -> [ meta, bam, bai ]
        }
        .set { ch_bam_bai }
    ch_pysamstats_input = ch_bam_bai.join(SAMTOOLS_CONSENSUS.out.consensus)

    SAMTOOLS_STATS (
        ch_bam_bai
    )

    PYSAMSTATS (
        ch_pysamstats_input
    )

    PYSAMSTATS_PARSER (
        PYSAMSTATS.out.pysamstats_tsv
    )

    REPORT_CHECK_COVERAGE (
        PYSAMSTATS_PARSER.out.amb_fastas,
        ch_fasta
    )

    REPORT_STANFORD_JSON (
        PYSAMSTATS_PARSER.out.fastas
    )

    REPORT_STANFORD_SUMMARY (
        REPORT_STANFORD_JSON.out.dr_json
    )

    SAMTOOLS_DEPTH_1 (
        BWA_MEM_2.out.bam
    )

    REPORT_PLOT_DEPTH_1 (
        SAMTOOLS_DEPTH_1.out.tsv
    )

    ch_depth_and_stats = Channel.empty()
    ch_depth_and_stats = SAMTOOLS_STATS.out.stats.join(SAMTOOLS_DEPTH_1.out.tsv)

    REPORT_SAMTOOLS_DEPTH_AND_STATS (
        ch_depth_and_stats
    )

    //create channel for all files needed for ngs_stats
    ch_ngs_stats = Channel.empty()
    ch_ngs_stats = TRIMMOMATIC.out.summary.join(BBTOOLS_SEAL.out.sealstats)
    ch_ngs_stats = ch_ngs_stats.join(BWA_MEM_2.out.bam)
    ch_ngs_stats = ch_ngs_stats.join(PYSAMSTATS_PARSER.out.fastas)
    ch_ngs_stats = ch_ngs_stats.join(REPORT_SAMTOOLS_DEPTH_AND_STATS.out.samtools_depth_and_stats)
    ch_ngs_stats = ch_ngs_stats.join(REPORT_CHECK_COVERAGE.out.contig_coverage_stats)

    REPORT_NGS_STATS_PASS (
        ch_ngs_stats
    )

    //Realign cleaned reads to reference to check coverage
    BWA_INDEX_REF (
        ch_fasta
    )

    BWA_MEM_REF (
        TRIMMOMATIC.out.trimmed_reads,
        BWA_INDEX_REF.out.index,
        true
    )

    SAMTOOLS_DEPTH_REF (
        BWA_MEM_REF.out.bam
    )

    REPORT_PLOT_DEPTH_REF (
        SAMTOOLS_DEPTH_REF.out.tsv
    )
    // end realign part
        ch_similarity = PYSAMSTATS_PARSER.out.fastas.collect{it[1]}.ifEmpty([])

    REPORT_SIMILARITY (
        ch_similarity
    )
    REPORT_SIMILARITY_MATRIX (
        REPORT_SIMILARITY.out.similarity_stats
    )

    ch_try_combine_ngs_stats = Channel.empty()
    ch_try_combine_ngs_stats = REPORT_NGS_STATS_ALL.out.ngs_stats.concat(REPORT_NGS_STATS_PASS.out.ngs_stats).groupTuple(by: 0)
    ch_try_combine_ngs_stats = ch_try_combine_ngs_stats.collect{it[1]}

    ch_final_report = Channel.empty()
    ch_final_report = ch_final_report.mix(ch_try_combine_ngs_stats.ifEmpty([]))
    if (params.positive_control_ref)
        {
            ch_final_report = ch_final_report.mix(CONTROL_TEST.out.report.collect{it[1]}.ifEmpty([]))
        }
    if (params.rundir)
        {
            ch_final_report = ch_final_report.mix(INTEROP_STATS.out.interop_stats.ifEmpty([]))
        }
    if(params.metadata)
        {
            ch_final_report = ch_final_report.mix(ch_metadata.ifEmpty([]))
        }

    REPORT_SUMMARY(
        ch_final_report.collect()
    )

    db_samtools_stats = SAMTOOLS_STATS.out.stats.collect{it[1]}
    db_jsons = REPORT_STANFORD_JSON.out.dr_json.collect{it[1]}

// optionally import data to MySQL db via sqlalchemy (/bin/db_import.py)
// if (params.rundir)
//     {
//     DATABASE_IMPORT (
//         ch_rundir,
//         ch_final_report.collect(),
//         db_samtools_stats,
//         db_jsons
//     )
//     }
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowTestpipeline.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
