
include { ILLUMINA_GET_INTEROP } from '../../modules/local/illumina/get_interop/main'
include { ILLUMINA_INTEROP_PARSE } from '../../modules/local/illumina/interop_parse/main'

// include { module to parse csv for stats } from ''

workflow INTEROP_STATS {
    take:
    rundir // dir: /path/to/run

    main:
    
    ILLUMINA_GET_INTEROP(
        rundir
    )

    ILLUMINA_INTEROP_PARSE(
        ILLUMINA_GET_INTEROP.out.interop_csv
    )


    // MODULE_TO_PARSE_FOR_STATS (interop.csv)
    // end with csv of interop stats we need

    emit:
    interop_stats       = ILLUMINA_INTEROP_PARSE.out.interop_stats_csv                        // channel: [ [csv] ]

}
