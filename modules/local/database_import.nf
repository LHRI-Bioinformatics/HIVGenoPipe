process DATABASE_IMPORT {

container = null
conda = "/your/local/env/here"

input:
path(rundir)
path(stats_files)
path(samtools_stats)
path(json_files)

output:
path("database_import/logs/db_import.err"), optional:true, emit: db_error
path("database_import/logs/db_import.log"), optional:true, emit: db_log

script:
"""
mkdir -p database_import/logs
    log_file=database_import/logs/db_import.log
    err_file=database_import/logs/db_import.err

for file in ${samtools_stats}
do
    cat \$file | grep ^SN | cut -f 2- > \$file.summary_only
done

db_import.py ${rundir}/ ${stats_files} -s *.summary_only -j ${json_files} 2>> \$err_file >> \$log_file

"""













}
