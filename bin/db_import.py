#!/usr/bin/env python3

import os
import json
import pandas as pd
import xml.etree.ElementTree as ET
from datetime import datetime
from collections import OrderedDict
import sys
import re
from os.path import exists
import argparse
import sqlalchemy
from sqlalchemy import create_engine

# start with final report script, remove all QC passing checks, just get numbers needed for db

############### dfSeqRun ###############

#function to format interop, return this section
    # interop db import is done directly from file so no setup needed here
    # has average insert size and average read length as mean for all samples
    # comes from test section?

def get_interop_stats(interop_file, average_insert_size, average_read_length):

    interop_value_list = []
    for file in interop_file:
        df=pd.read_csv(file)
        # print(df)
        interop_value_list = df['Value'].tolist()

    cluster_density = interop_value_list[0]
    Q30_score = interop_value_list[1]
    cluster_pf = interop_value_list[2]
    interop_yield = interop_value_list[3]

    # use to make new formatted df
    interop_QC_dict = {
    "Yield": [interop_yield],
    "ClusterDensity": [cluster_density],
    "ClustersPF": [cluster_pf],
    "Q30": [Q30_score],
    "average_insert_size": [average_insert_size],
    "average_read_length": [average_read_length],
    }

    df_interop = pd.DataFrame(interop_QC_dict)



    return(df_interop)

#function to format positive control test sample(s), return this section
    # add QC ranges for this section
def get_pos_control_stats(pos_control_file):

    # make QC dict or list here for adding to df with formatting

    dict_of_pos_control_dfs = {}
    # loop through in case of multiple positive controls
    for file in pos_control_file:
        sample_name = str(file).rsplit("_", 1)[0]
        # print("Sample name " + sample_name)
        dict_of_pos_control_dfs[sample_name] = pd.read_csv(file, delimiter= "\t")
        # remove extra stats
        dict_of_pos_control_dfs[sample_name] = dict_of_pos_control_dfs[sample_name].drop(columns=(['Sample','Total Raw Read1', 'Total HIV Read1(% of raw)','Sample Type','Average Insert Size (bp)', 'Average Read Length (bp)',
                                                                                                'Bases Mapped Count', 'Average Quality Score' ]))
        # move error rate metric
        column_to_move = dict_of_pos_control_dfs[sample_name].pop("Error Rate (Read Variation %)")
        dict_of_pos_control_dfs[sample_name].insert(7, "Error Rate (Read Variation %)", column_to_move )

        #rename for db
        dict_of_pos_control_dfs[sample_name] = dict_of_pos_control_dfs[sample_name].rename(columns={'sampleId':'pos_ctrl', 'Total HIV Read1':'pos_ctrl_filtered_hiv_reads', 'Median Depth (reads)':'pos_ctrl_median_depth',
                                                                                                    'Depth MAD (reads)':'pos_ctrl_depth_mad', 'Contig Coverage (%)':'pos_ctrl_coverage', 'Error Rate (Read Variation %)':'pos_ctrl_error_rate',
                                                                                                    'Contig Error Rate at 5% (%)':'pos_ctrl_error_rate_Amb5', 'Contig Error Rate at 15% (%)':'pos_ctrl_error_rate_Amb15'})

        # check column names, drop extra amb errors
        list_of_columns = list(dict_of_pos_control_dfs[sample_name].columns)
        list_to_drop = []
        for metric in list_of_columns:
            # only need 5 and 15%
            if metric.startswith('Contig Error'):
                list_to_drop.append(metric)
        dict_of_pos_control_dfs[sample_name]=dict_of_pos_control_dfs[sample_name].drop(columns=(list_to_drop))

    return(dict_of_pos_control_dfs)


# function to format negative control info
    # only need filtered HIV reads and Contig coverage
def get_neg_control_stats(negative_control_file):
    for file in negative_control_file:
        df=pd.read_csv(file, delimiter= "\t")
        df = df.drop(columns=(['Total Raw Read1', 'Total HIV Read1(% of raw)', 'Sample']))

        df = df.rename(columns={'sampleId':'neg_ctrl', 'Total HIV Read1':'neg_ctrl_filtered_hiv_reads','Contig Coverage (%)':'neg_ctrl_coverage'})

    return(df)

def get_runparams(rta_file, runparam_file,seq_run_id):
    dfCompletion=pd.read_csv(rta_file,header=None)

    #### the etree reading in is messy right now, need to find better way to read in or manually format xml file with etree package
    runParametersXML = ET.parse(runparam_file)
    runParameters = runParametersXML.getroot()
    # https://github.com/bartongroup/AlmostSignificant/blob/master/almostSignificant/dataLoading/dataLoading.py
    runParametersDict = {}
    #run id
    for child in runParameters.iterfind("RunID"):
        runParametersDict["RunID"] = child.text
    #run start date
    for child in runParameters.iterfind("RunStartDate"):
        runParametersDict["RunStartDate"] = child.text
    #run param ver
    for child in runParameters.iterfind("RunParametersVersion"):
        runParametersDict["RunParametersVersion"] = child.text
    #scanner id
    for child in runParameters.iterfind("ScannerID"):
        runParametersDict["ScannerID"] = child.text
    #run num
    for child in runParameters.iterfind("RunNumber"):
        runParametersDict["RunNumber"] = child.text
    #fpga
    for child in runParameters.iterfind("FPGAVersion"):
        runParametersDict["FPGAVersion"] = child.text
    #MCSVersion
    for child in runParameters.iterfind("MCSVersion"):
        runParametersDict["MCSVersion"] = child.text
    #RTAVersion
    for child in runParameters.iterfind("RTAVersion"):
        runParametersDict["RTAVersion"] = child.text
    #Barcode
    for child in runParameters.iterfind("Barcode"):
        runParametersDict["Barcode"] = child.text
    #PR2BottleBarcode
    for child in runParameters.iterfind("PR2BottleBarcode"):
        runParametersDict["PR2BottleBarcode"] = child.text
    #ReagentKitPartNumberEntered
    for child in runParameters.iterfind("ReagentKitPartNumberEntered"):
        runParametersDict["ReagentKitPartNumberEntered"] = child.text
    #ReagentKitVersion
    for child in runParameters.iterfind("ReagentKitVersion"):
        runParametersDict["ReagentKitVersion"] = child.text
    #ReagentKitBarcode
    for child in runParameters.iterfind("ReagentKitBarcode"):
        runParametersDict["ReagentKitBarcode"] = child.text
    #PreviousPR2BottleBarcode
    for child in runParameters.iterfind("PreviousPR2BottleBarcode"):
        runParametersDict["PreviousPR2BottleBarcode"] = child.text
    #ExperimentName
    for child in runParameters.iterfind("ExperimentName"):
        runParametersDict["ExperimentName"] = child.text
    #RunDescription
    for child in runParameters.iterfind("RunDescription"):
        runParametersDict["RunDescription"] = child.text
    #LibraryPrepKit
    for child in runParameters.iterfind("LibraryPrepKit"):
        runParametersDict["LibraryPrepKit"] = child.text
    #IndexKit
    for child in runParameters.iterfind("IndexKit"):
        runParametersDict["IndexKit"] = child.text
    #Chemistry
    for child in runParameters.iterfind("Chemistry"):
        runParametersDict["Chemistry"] = child.text
    #Username
    for child in runParameters.iterfind("Username"):
        runParametersDict["Username"] = child.text
    #SampleSheetName
    for child in runParameters.iterfind("SampleSheetName"):
        runParametersDict["SampleSheetName"] = child.text
    #FocusMethod
    for child in runParameters.iterfind("FocusMethod"):
        runParametersDict["FocusMethod"] = child.text
    #SelectedRun
    for child in runParameters.iterfind("SelectedRun"):
        runParametersDict["SelectedRun"] = child.text
    #ModuleName
    for child in runParameters.iterfind("ModuleName"):
        runParametersDict["ModuleName"] = child.text
    #ModuleVersion
    for child in runParameters.iterfind("ModuleVersion"):
        runParametersDict["ModuleVersion"] = child.text
    #MAPRunId
    for child in runParameters.iterfind("MAPRunId"):
        runParametersDict["MAPRunId"] = child.text


    ## parts needed: ['RunID','RunStartDate','RunParametersVersion','ScannerID','RunNumber','FPGAVersion','MCSVersion','RTAVersion','Barcode','PR2BottleBarcode','ReagentKitPartNumberEntered','ReagentKitVersion','ReagentKitBarcode',
        #'PreviousPR2BottleBarcode','ExperimentName','RunDescription','LibraryPrepKit','IndexKit','Chemistry','Username','SampleSheetName','FocusMethod','SelectedRun','ModuleName','ModuleVersion','MAPRunId']

    dfRunParams = pd.DataFrame([runParametersDict])

    dfRunParams = dfRunParams.rename(columns={'RunID':'seq_run_id', 'RunStartDate':'assay_date'})
    ## make check for columns because old format

    # dfRunParams.columns=['seq_run_id','assay_date','RunParametersVersion','ScannerID','RunNumber','FPGAVersion','MCSVersion','RTAVersion','Barcode','PR2BottleBarcode','ReagentKitPartNumberEntered','ReagentKitVersion','ReagentKitBarcode','PreviousPR2BottleBarcode','ExperimentName','RunDescription','LibraryPrepKit','IndexKit','Chemistry','Username','SampleSheetName','FocusMethod','SelectedRun','ModuleName','ModuleVersion','MAPRunId']
    dfRunParams['seq_run_id']=seq_run_id
    dfRunParams['assay_date']=pd.to_datetime(dfRunParams['assay_date'],format='%y%m%d')
    dfRunParams['completion_date']=pd.to_datetime(dfCompletion[0])

    print("THIS IS RUNPARAMS")
    print(dfRunParams)

    return dfRunParams

# seqrun df includes neg control, pos control, interop, and run params
## get controls and interop from final report script
## add run params controls (use exact file names)
def build_dfseqrun(df_interop, df_neg, df_pos, runparams_df):
    # need to merge these together
    dfSeqRun=pd.concat([runparams_df,df_interop,df_neg,df_pos],axis=1)

    return(dfSeqRun)

################ dfSeqRunSample ###############

def getSamtoolsStats(samtools_stats_files):
    ## can we add grep ^SN | cut -f 2- in our new module for this db step then we can just copy old script
    dfSamToolsStats=None
    for File in samtools_stats_files:
        sample_name = os.path.basename(File)
        seq_run_sample_id = str(sample_name).rsplit(".")[0]
        #remove _T suffix
        seq_run_sample_id = str(seq_run_sample_id).rsplit("_",1)[0]

        SamtoolsStatsDF = pd.read_csv(File, sep="\t", usecols=[0, 1], names=['Metric', seq_run_sample_id], header=None)
        SamtoolsStatsDF['Metric']=SamtoolsStatsDF['Metric'].str.rstrip(':')
        SamtoolsStatsDFTransposed=SamtoolsStatsDF.transpose()
        columns=SamtoolsStatsDFTransposed.loc['Metric']
        SamtoolsStatsDFTransposed.columns = columns
        SamtoolsStatsDFTransposed.drop('Metric',inplace=True)
        if dfSamToolsStats is None:
            dfSamToolsStats=SamtoolsStatsDFTransposed
        else:
            dfSamToolsStats=pd.concat([dfSamToolsStats,SamtoolsStatsDFTransposed])

    print("THIS IS SAMTOOLS STATS")
    print(dfSamToolsStats)

    return dfSamToolsStats

    # get column names from ngs_stats script



#function to format each sample, need to get negative control stats from sample type here because negative samples dont have separate stats
def get_test_sample_stats(test_sample_files):
    dict_of_test_dfs = {}
    for file in test_sample_files:
        # cut off the _T1
        sample_name = str(file).rsplit("_", 2)[0]
        print(sample_name)
        dict_of_test_dfs[sample_name] = pd.read_csv(file, delimiter= "\t")
        # add this because the negatives dont have sample type added from pipeline
        if file.startswith('negative'):
            dict_of_test_dfs[sample_name]['Sample Type'] = 'negative'

    build_df = pd.concat(dict_of_test_dfs.values(), ignore_index=True)
    build_df = build_df.groupby(build_df["Sample"])
    final_df = build_df.first()



    final_df = final_df.reset_index()

    # cut off the _T1
    final_df['Sample'] = final_df['Sample'].str.rsplit("_", n=1).str[0]
    final_df['sampleId'] = final_df['sampleId'].str.rsplit("_", n=1).str[0]


    average_insert_size = final_df['Average Insert Size (bp)'].mean()
    average_read_length = final_df['Average Read Length (bp)'].mean()

    print("THIS IS TEST SAMPLE STATS")
    print(final_df)

    return(final_df, average_insert_size, average_read_length)

def build_dfseqrunsample(test_samples_df, samtools_df,seq_run_id_long):

    dfTestSamples = test_samples_df.copy()
    dfTestSamples['seq_run_id'] = seq_run_id_long
    # time script is run, OK to keep?
    analysis_date=pd.to_datetime(datetime.now().strftime('%Y-%m-%d'))

    # seq_run_sample_id = rundir + ~ + full sample ID (up to _S#)
    # sample_id = P number only if available
    # seq_run_id = rundir only
    print("THIS IS TEST SAMPLES BEFORE RENAME")
    print(dfTestSamples)

    dfSeqRunSample = dfTestSamples[['seq_run_id','Sample', 'sampleId','Total Raw Read1', 'Total HIV Read1','Filtered HIV Read1','Short HIV Read1',
                                    'Non-properly paired HIV Read1','consensus Length', "Numconsensus N's",'Contig Coverage (%)','Median Depth (reads)', 'Depth MAD (reads)',]].copy()

    dfSeqRunSample = dfSeqRunSample.rename(columns={'Sample':'seq_run_sample_id', 'sampleId':'sample_id', 'Total Raw Read1':'raw_pairs', 'Total HIV Read1':'hiv_pairs',
                                                    'Filtered HIV Read1':'filtered_hiv_pairs','Short HIV Read1':'short_hiv_pairs','Non-properly paired HIV Read1':'non_prop_paired_hiv_pairs',
                                                    'consensus Length':'consensus_length', "Numconsensus N's":'num_consensus_n','Contig Coverage (%)':'coverage',
                                                    'Median Depth (reads)':'median_depth', 'Depth MAD (reads)':'depth_mad'})
    print("THIS IS SEQRUNSAMPLE before samtools merge")
    print(dfSeqRunSample)
    # merge samtools_df
    samtools_df.reset_index(inplace=True)
    samtools_df = samtools_df.rename(columns = {'index':'seq_run_sample_id'})
    print("THIS IS SAMTOOLS DF BEFORE MERGE")
    print(samtools_df)
    # print(dfSamToolsStats)
    # print(dfSamToolsStats.columns)
    dfSeqRunSample = dfSeqRunSample.merge(samtools_df, on=['seq_run_sample_id'], how="left")

    # add analyst name based on name of sample
    dfAttributes = dfSeqRunSample['seq_run_sample_id'].str.split('-', expand = True)
    if len(dfAttributes.columns)==5:
        dfAttributes=dfAttributes[[2,4]]
        dfAttributes.columns=['amplicon','analyst']
        dfAttributes['analyst']=dfAttributes['analyst'].str.split('_').str[0].str.replace("0","")
        # print(dfAttributes)
        dfSeqRunSample=pd.concat([dfSeqRunSample,dfAttributes],axis=1)

    # add analysis date
    dfSeqRunSample['analysis_date']=analysis_date
    #adjust seq run id to short
    seqRunIDParts=seq_run_id_long.split("_")
    seq_run_id = seqRunIDParts[0] + "_" + seqRunIDParts[1]
    dfSeqRunSample['seq_run_sample_id']=seq_run_id+"~"+dfSeqRunSample['seq_run_sample_id']

    newColumns=[]
    for column in dfSeqRunSample.columns:
        newColumns.append(column.replace(" (%)","").replace(" ","_").lower())
    dfSeqRunSample.columns=newColumns

    return dfSeqRunSample

############### dfHIVSequences ###############

def build_dfhivseq(test_samples_df, seq_run_id):
    # seq_id example: 230208_M70177~3405-211347604-E-157-Ayub_S9.consensus
    # seq_run_id example: 230208_M70177~3405-211347604-E-157-Ayub_S9

    dfTestSamples = test_samples_df.copy()

    dfTestSamples['Sample']=seq_run_id+"~"+dfTestSamples['Sample']

    # currently coming from _stats file time, for now just put sys time
    analysis_date=pd.to_datetime(datetime.now().strftime('%Y-%m-%d'))
    dfTestSamples['analysis_date']=analysis_date

    dfAmb5HIVSequences=dfTestSamples[['Sample', 'Amb5 Length', 'NumAmb5 N\'s', 'Amb5 Seq','analysis_date']].copy()
    dfAmb5HIVSequences.columns=['seq_run_sample_id', 'nt_length', 'num_nocalls', 'NT','analysis_date']
    dfAmb5HIVSequences['seq_id']=dfAmb5HIVSequences['seq_run_sample_id']+'.Amb5'
    dfAmb5HIVSequences['AmbFreq']='Amb5'

    dfAmb15HIVSequences=dfTestSamples[['Sample', 'Amb15 Length', 'NumAmb15 N\'s', 'Amb15 Seq','analysis_date']].copy()
    dfAmb15HIVSequences.columns=['seq_run_sample_id', 'nt_length', 'num_nocalls', 'NT','analysis_date']
    dfAmb15HIVSequences['seq_id']=dfAmb5HIVSequences['seq_run_sample_id']+'.Amb15'
    dfAmb15HIVSequences['AmbFreq']='Amb15'

    dfConsensusHIVSequences=dfTestSamples[['Sample', 'consensus Length', 'Numconsensus N\'s', 'consensus Seq','analysis_date']].copy()
    dfConsensusHIVSequences.columns=['seq_run_sample_id', 'nt_length', 'num_nocalls', 'NT','analysis_date']
    dfConsensusHIVSequences['seq_id']=dfAmb5HIVSequences['seq_run_sample_id']+'.consensus'
    dfConsensusHIVSequences['AmbFreq']='consensus'

    dfHIVSequences=pd.concat([dfAmb15HIVSequences,dfAmb5HIVSequences,dfConsensusHIVSequences])
    dfHIVSequences=dfHIVSequences.dropna(subset=['NT'])


    return dfHIVSequences

############### dfProteinSequences ###############

# pipe output from stanford_summary step
# seq_id = 240326_M07814~3727-P247579-1266-A-A12-Ayub_S11.consensus (from dfHIVSequences)
def build_stanford_dfs(jsonFile, seq_run_id):
    analysis_date=os.path.getmtime(jsonFile)
    # print("analysis_date: {}".format(analysis_date))
    # print("Reading characterization data from Stanford")
    with open(jsonFile, "r") as jsonIn:
        # Writing data to a file
        data = json.load(jsonIn)
        jsonIn.close()  # to change file access modes

    hivdb_version=data['currentVersion']['text']
    sierra_version=data['currentProgramVersion']['text']

    dfMutations=pd.DataFrame(columns = ['seq_id','protein','reference','mutation','mutation_site','codon','insertedNAs','isInsertion','isDeletion','isIndel','isAmbiguous','isApobecMutation','isApobecDRM','hasStop','isUnusual','mutation_type','remarks','sierra_version'])

    dfDRMutations=pd.DataFrame(columns = ['seq_id', 'drug_class', 'drug_name', 'drug_abbr', 'score', 'level', 'interpretation', 'mutations','hivdb_version'])

    dfProteinSequences=pd.DataFrame(columns = ['seq_id', 'protein', 'alignedNT','alignedAA', 'clade', 'matchPcnt','sierra_version'])

    for sequence in data['sequenceAnalysis']:
        # need to reformat seq_id
        seqName=sequence['inputSequence']['header'] #3627-P245364-A-11A-Ayub_S15_T1.Amb15
        seqname_keep = seqName.rsplit("_",1)
        seqname_keep_amb_part = seqname_keep[-1]
        amb_part = seqname_keep_amb_part.split(".")
        seqName = seq_run_id + "~" + str(seqname_keep[0]) + "." + str(amb_part[-1])

        clade=sequence['bestMatchingSubtype']['displayWithoutDistance']

        for gene in sequence['alignedGeneSequences']:
            # print(gene)
            dfProteinSequences.loc[len(dfProteinSequences)] = [
                        seqName,
                        gene['gene']['name'],
                        gene['alignedNAs'],
                        gene['alignedAAs'],
                        clade,
                        gene['matchPcnt'],
                        sierra_version]

            for mutation in gene['mutations']:
                # print("{}\t{}({})\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                comments=[]
                for comment in mutation["comments"]:
                    comments.append(comment["text"])
                dfMutations.loc[len(dfMutations)] = [
                        seqName,
                        # "{}({})".format(gene['gene']['name'],gene['gene']['length']),
                        gene['gene']['name'],
                        # gene['firstAA'],
                        # gene['lastAA'],
                        # "{}{}{}".format(mutation['reference'],mutation['position'],mutation['AAs']),
                        mutation['reference'],
                        mutation['AAs'],
                        mutation['position'],
                        mutation['triplet'],
                        mutation["insertedNAs"],
                        mutation["isInsertion"],
                        mutation["isDeletion"],
                        mutation["isIndel"],
                        mutation["isAmbiguous"],
                        mutation["isApobecMutation"],
                        mutation["isApobecDRM"],
                        mutation["hasStop"],
                        mutation["isUnusual"],
                        mutation["primaryType"],
                        " ".join(comments),
                        sierra_version
                    ]

        # print(dfMutations.dtypes)
        for col in ['isInsertion','isDeletion','isIndel','isAmbiguous', 'isApobecMutation','isApobecDRM','hasStop','isUnusual']:
            dfMutations[col] = dfMutations[col].astype('bool')
        # print(dfMutations.dtypes)

        for dr in sequence['drugResistance']:
            for drugScore in dr['drugScores']:
                dr_mutations=[]
                for partialScore in drugScore['partialScores']:
                    for dr_mutation in partialScore['mutations']:
                        if dr_mutation['text'] not in dr_mutations:
                            dr_mutations.append(dr_mutation['text'])
                # print(", ".join(dr_mutations))

                dfDRMutations.loc[len(dfDRMutations)] = [
                        seqName,
                        drugScore['drugClass']['name'],
                        drugScore['drug']['fullName'],
                        drugScore['drug']['displayAbbr'],
                        drugScore['score'],
                        drugScore['level'],
                        drugScore['text'],
                        ", ".join(dr_mutations),
                        hivdb_version
                        ]

    return [dfProteinSequences,dfMutations,dfDRMutations]

############### mergeDB stuff ###############
def mergeWithDB(engine,dbTable,mergeColumns,df):
    # print("############## mergeWithDB  #############")
    # print(df.columns)
    # print(df.T)
    dbColumnsStr=",".join(mergeColumns)
    sql="select distinct "+dbColumnsStr+" from "+dbTable
    # print(sql)
    dfDB=pd.read_sql(sql,engine)
    # print(dfDB.columns)
    # print(dfDB)
    df = df.merge(dfDB, on=mergeColumns, how="left", indicator=True)
    df = df[df["_merge"] == "left_only"].drop(columns=["_merge"])

    return df


def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():

    parser = argparse.ArgumentParser(description="Create a summary report using various statistics files generated from HIV Genopipe")
    parser.add_argument('RunDir', help='path to Run level directory ex: /mnt/VISL/Instruments/MiSeq/240326_M07814_0037_000000000-DMR4L')
    parser.add_argument('InputFile', nargs='+', type=File, help='group of stats files')
    parser.add_argument('-s', '--samtools_stats', nargs='+', type=File, help='group of samtools stats files')
    parser.add_argument('-j', '--jsons', nargs='+', type=File, help='group of json files from Stanford DR')
    args = parser.parse_args()

    run_dir = args.RunDir
    # get seq_run_id part from run dir
    seq_run_id_long = os.path.basename(os.path.normpath(run_dir))
    seqRunIDParts=seq_run_id_long.split("_")
    seq_run_id = seqRunIDParts[0] + "_" + seqRunIDParts[1]

    interop_file = []
    pos_control_file = []
    negative_control_file = []
    test_sample_files = []

    metadata_file = []

    #parse through all input files and sort them by type
    for file in args.InputFile:
        # print(file)
        if file.startswith('control'):
            pos_control_file.append(file)
        elif file.startswith('interop'):
            interop_file.append(file)
        elif file.startswith('negative'):
            negative_control_file.append(file)
        elif file.startswith('metadata'):
            metadata_file.append(file)
        elif file.startswith('RTAComplete'):
            rta_file = file
        elif file.startswith('RunParameters'):
            runparam_file = file
        else:
            test_sample_files.append(file)
    test_sample_files= test_sample_files + negative_control_file

    # search for repeat samples in test samples and remove the prestats version if final stats are available
    prestats = []
    poststats = []

    for sample in test_sample_files:
        # prestats have a prefix in the name. post stats do not
        if sample.startswith('positive') | sample.startswith('test'):
            prestats.append(sample)
        else:
            poststats.append(sample)


    samples_in_pre_not_in_post = []
    for sample in prestats:
        sample_rootname = sample.removeprefix("positive_")
        sample_rootname = sample_rootname.removeprefix("test_")
        # print(sample_rootname)
        if sample_rootname not in poststats:
            samples_in_pre_not_in_post.append(sample)

    test_sample_files = samples_in_pre_not_in_post + poststats


        #need test samples first for interop input
    test_samples_df = get_test_sample_stats(test_sample_files)

    ##### build seqrun #####
    # set up path to run params files
    rta_file = run_dir + "RTAComplete.txt"
    ## change for old files
    runparam_file = run_dir + "runParameters.xml"

    runparams_df = get_runparams(rta_file,runparam_file,seq_run_id_long)

    if not interop_file:
        interop_stats_df = pd.DataFrame()
    else:
        interop_stats_df = get_interop_stats(interop_file, test_samples_df[1], test_samples_df[2])


    if not negative_control_file:
        # make empty df
        neg_control_df = pd.DataFrame()
    else:
        neg_control_df = get_neg_control_stats(negative_control_file)

    if not pos_control_file:
        pos_control_df = pd.DataFrame()
    else:
        pos_control_dfs = get_pos_control_stats(pos_control_file) # is a dict of dfs
        #currently handles more than one pos control but we only need one
        pos_control_df = next(iter(pos_control_dfs.values()))


    dfSeqRun = build_dfseqrun(interop_stats_df, neg_control_df,pos_control_df,runparams_df)

    print("THIS IS SEQRUN")
    print(dfSeqRun)

    ##### build hivseq #####
    dfHIVSequences = build_dfhivseq(test_samples_df[0], seq_run_id)
    print("THIS IS dfHIVSequences")
    print(dfHIVSequences)

    ##### build seqrunsample #####
    samtools_df = getSamtoolsStats(args.samtools_stats)

    dfSeqRunSample = build_dfseqrunsample(test_samples_df[0],samtools_df,seq_run_id_long)

    print("THIS IS SEQRUNSAMPLE COMPLETED")
    print(dfSeqRunSample)

    ##### build stanford dfs #####
    mutationFrames=[]
    drMutationFrames=[]
    proteinSequencesFrames=[]
    for json in args.jsons:
        sierraDFs = build_stanford_dfs(json,seq_run_id)
        proteinSequencesFrames.append(sierraDFs[0])
        mutationFrames.append(sierraDFs[1])
        drMutationFrames.append(sierraDFs[2])
    dfProteinSequences=pd.concat(proteinSequencesFrames)
    dfHIVMutations=pd.concat(mutationFrames)
    dfHIVDrugResistance=pd.concat(drMutationFrames)

    print("THIS IS STANFORD TABLES")
    print(dfProteinSequences)
    print(dfHIVMutations)
    print(dfHIVDrugResistance)



    ############### db engine #################
    engine = create_engine("mysql+pymysql://{user}:{pw}@{host}/{db}"
                        .format(host='00.0.0.00', db='DB_NAME_HERE', user='USER_NAME_HERE', pw='USER_PASSWORD_HERE'))

    ############### sample table check #################
    sample_name_list = dfSeqRunSample['sample_id'].tolist()

    queryString = 'SELECT * FROM sample WHERE sample_id in '+str(tuple(sample_name_list))
    dfValidSamples=pd.read_sql(queryString,engine)
    samples_in_db = dfValidSamples['sample_id'].tolist()
    print("THESE ARE MATCHING SAMPLES")
    print(samples_in_db)

    # get samples not in db
    samples_to_add = list(set(sample_name_list) - set(samples_in_db))

    if len(samples_to_add) == 0:
        loadDB=True
    else:
        print("samples to add: ", samples_to_add)

        #build df to add to Samples
        dfSamplesToAdd = test_samples_df[0][['sampleId','Sample Type']]
        print("samples to add DF from samples df")
        print(dfSamplesToAdd)
        dfSamplesToAdd.columns=['sample_id', 'sample_feature']
        print("rename columns")
        print(dfSamplesToAdd)
        print("get only samples to add")
        dfSamplesToAdd = dfSamplesToAdd[dfSamplesToAdd["sample_id"].isin(samples_to_add)]
        print(dfSamplesToAdd)

        dfSamplesToAdd.drop_duplicates(inplace=True)
        dfSamplesToAdd=dfSamplesToAdd.dropna(subset=['sample_id'])

        #merge to DB
        dfSamplesToAdd = mergeWithDB(engine,'sample',['sample_id'],dfSamplesToAdd)
        dfSamplesToAdd.to_sql('sample', engine, if_exists='append', index=False)
        loadDB=True

    ## see which ones are NOT in sample table and add them with info we have (sample type)

    ############### merge DB ###############

    if loadDB:
        try:
            # import each df into the db
            dfSeqRun=mergeWithDB(engine,'seq_run',['seq_run_id'],dfSeqRun)
            dfSeqRun.to_sql('seq_run', engine, if_exists='append', index=False)

            dfSeqRunSample = mergeWithDB(engine,'seq_run_sample',['seq_run_sample_id'],dfSeqRunSample)
            dfSeqRunSample.to_sql('seq_run_sample', engine, if_exists='append', index=False)

            dfHIVSequences = mergeWithDB(engine,'hiv_sequences',['seq_id','seq_run_sample_id'],dfHIVSequences)
            dfHIVSequences.to_sql('hiv_sequences', engine, if_exists='append', index=False)

            dfProteinSequences = mergeWithDB(engine,'hiv_protein_sequences',['seq_id','protein'],dfProteinSequences)
            dfProteinSequences.to_sql('hiv_protein_sequences', engine, if_exists='append', index=False)

            dfHIVMutations = mergeWithDB(engine,'hiv_mutations',['seq_id','protein','mutation_site','mutation'],dfHIVMutations)
            dfHIVMutations.to_sql('hiv_mutations', engine, if_exists='append', index=False)

            dfHIVDrugResistance = mergeWithDB(engine,'hiv_drug_resistance',['seq_id','drug_name'],dfHIVDrugResistance)
            dfHIVDrugResistance.to_sql('hiv_drug_resistance', engine, if_exists='append', index=False)
        except Exception as error:
            print("Data import has failed, run has not been loaded to DB. Please check the error and rerun")
            print(error)
            if error.__class__.__name__ == "DataError":
                #check mutations length
                print("Stanford may have included an extended mutation:")
                dfHIVMutations["mutations length"] = dfHIVMutations["mutation"].str.len()
                new_dfHIVMutations = dfHIVMutations[dfHIVMutations['mutations length'] >= 100]
                print(new_dfHIVMutations)







if __name__ == "__main__":
    main()
