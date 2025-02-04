#!/usr/bin/env python3

import os
import pandas as pd
from collections import OrderedDict
import sys
import re
from os.path import exists
import argparse

#function to format interop, return this section
    # interop db import is done directly from file so no setup needed here
    # has average insert size and average read length as mean for all samples
    # comes from test section?
    # add QC ranges for this section; dictionary of metric and QC range?

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



    # add params to change these, hard code for now
    yield_range = 200000000
    lower_cd_range = 300
    higher_cd_range = 1900
    cluster_passing_perc_range = 75
    Q30_score_range = 75
    average_insert_size_range = 200
    average_read_length_range = 150


    # use to make new formatted df
    interop_QC_dict = {
    "Yield (bp)": [interop_yield, ">=" + str(yield_range)],
    "Cluster Density (K/mm2)": [cluster_density, str(lower_cd_range) + "-" + str(higher_cd_range)],
    "Cluster Passing Filter (%)": [cluster_pf, ">=" + str(cluster_passing_perc_range) + "%"],
    "Q30 score (%)": [Q30_score, ">=" + str(Q30_score_range) + "%"],
    "Average Insert Size (bp)": [average_insert_size, ">=" + str(average_insert_size_range)],
    "Average Read Length (bp)": [average_read_length, ">=" + str(average_read_length_range)],
    }

    # compare columns to see if passing and add to df, need to use variable not df item because string
    interop_check = []
    if interop_yield >= yield_range:
        interop_check.append("PASSED")
    else:
        interop_check.append("FAILED")

    if lower_cd_range <= cluster_density <= higher_cd_range:
        interop_check.append("PASSED")
    else:
        interop_check.append("FAILED")

    if cluster_pf >= cluster_passing_perc_range:
        interop_check.append("PASSED")
    else:
        interop_check.append("FAILED")

    if Q30_score >= Q30_score_range:
        interop_check.append("PASSED")
    else:
        interop_check.append("FAILED")

    if average_insert_size >= average_insert_size_range:
        interop_check.append("PASSED")
    else:
        interop_check.append("FAILED")

    if average_read_length >= average_read_length_range:
        interop_check.append("PASSED")
    else:
        interop_check.append("FAILED")

    # print(interop_check)



    # print(interop_QC_dict)

    df_interop = pd.DataFrame(interop_QC_dict)

    df_interop_t = df_interop.transpose()
    # df_interop_t["passing"] =  interop_check
    df_interop_t.insert(0, 'passing', interop_check)

    return(df_interop_t)

#function to format positive control test sample(s), return this section
    # add QC ranges for this section
def get_pos_control_stats(pos_control_file):

    # add params to change these, hard code for now
    # Amb 1 will have 1%
    # The rest contig error will have 0.1
    filtered_reads_range = 5000 #gte
    HIV_reads_perc_low_range = 50
    HIV_reads_perc_high_range = 100
    med_depth_low_range = 200
    med_depth_high_range = 100000000
    depth_MAD_low_range = 200
    depth_MAD_high_range = 100000000
    contig_coverage_range = 100 #perc
    read_error_rate_range = 1 #lt perc
    contig_error_rate_1_range = 1 #lt perc
    contig_error_rate_other_range = 0.1 #lt perc



    # make QC dict or list here for adding to df with formatting

    dict_of_pos_control_dfs = {}
    # loop through in case of multiple positive controls
    for file in pos_control_file:
        sample_name = str(file).rsplit("_", 1)[0]
        # print("Sample name " + sample_name)
        dict_of_pos_control_dfs[sample_name] = pd.read_csv(file, delimiter= "\t")
        # remove extra stats
        dict_of_pos_control_dfs[sample_name] = dict_of_pos_control_dfs[sample_name].drop(columns=(['Total Raw Read1', 'Sample Type','Average Insert Size (bp)', 'Average Read Length (bp)',
                                                                                                'Bases Mapped Count', 'Average Quality Score' ]))
        # move error rate metric
        column_to_move = dict_of_pos_control_dfs[sample_name].pop("Error Rate (Read Variation %)")
        dict_of_pos_control_dfs[sample_name].insert(7, "Error Rate (Read Variation %)", column_to_move )


        # check column names, contig error rate number can change
        list_of_columns = list(dict_of_pos_control_dfs[sample_name].columns)
        list_of_amb = []
        for metric in list_of_columns:
            if metric.startswith('Contig Error'):
                list_of_amb.append(metric)
        # print(list_of_amb)

        # get metrics in list to compare passing
        row_list = dict_of_pos_control_dfs[sample_name].loc[0, 'Total HIV Read1':].values.flatten().tolist()
        # print(row_list)

        # deal with this separately
        contig_error_list = row_list[6:]


        pos_control_check = []
        # will need to add two empty indices to the front for sample name spaces
        pos_control_check.append("")
        pos_control_check.append("")

        for i in range(len(row_list)):
            if i == 0:
                if row_list[i] >= filtered_reads_range:
                    pos_control_check.append("PASSED")
                else:
                    pos_control_check.append("FAILED")
            elif i == 1:
                if HIV_reads_perc_low_range <= row_list[i] <= HIV_reads_perc_high_range:
                    pos_control_check.append("PASSED")
                else:
                    pos_control_check.append("FAILED")
            elif i == 2:
                if med_depth_low_range <= row_list[i] <= med_depth_high_range:
                    pos_control_check.append("PASSED")
                else:
                    pos_control_check.append("FAILED")
            elif i == 3:
                if depth_MAD_low_range <= row_list[i] <= depth_MAD_high_range:
                    pos_control_check.append("PASSED")
                else:
                    pos_control_check.append("FAILED")
            elif i == 4:
                if row_list[i] < contig_coverage_range:
                    pos_control_check.append("PASSED")
                else:
                    pos_control_check.append("FAILED")
            elif i == 5:
                if row_list[i] < read_error_rate_range:
                    pos_control_check.append("PASSED")
                else:
                    pos_control_check.append("FAILED")

        # set up pos control ranges to add to df, two space for sample info
        pos_control_ranges = ["", "", ">=" + str(filtered_reads_range), str(HIV_reads_perc_low_range) + "-" + str(HIV_reads_perc_high_range),
            str(med_depth_low_range) + "-" + str(med_depth_high_range), str(depth_MAD_low_range) + "-" + str(depth_MAD_high_range),
            str(contig_coverage_range) + "%", "<" + str(read_error_rate_range) + "%"]


        # for the rest of the metrics, check if it's Amb 1%, else compare to 0.1%
        for amb_file in range(len(list_of_amb)):
            # print(amb_file)
            # print(list_of_amb[amb_file])
            # print(contig_error_list[amb_file])

            if list_of_amb[amb_file] == "Contig Error Rate at 1% (%)":
                # add the qc range to the list to add to final df
                pos_control_ranges.append("<" + str(contig_error_rate_1_range) + "%")
                if contig_error_list[amb_file] < contig_error_rate_1_range:
                    pos_control_check.append("PASSED")
                else:
                    pos_control_check.append("FAILED")
            else:
                # add the qc range to the list to add to final df
                pos_control_ranges.append("<" + str(contig_error_rate_other_range) + "%")
                if contig_error_list[amb_file] < contig_error_rate_other_range:
                    pos_control_check.append("PASSED")

                else:
                    pos_control_check.append("FAILED")

        # print(pos_control_check)
        # print(pos_control_ranges)


        # compare metrics to passing values

        # flip the dataframe
        dict_of_pos_control_dfs[sample_name] = dict_of_pos_control_dfs[sample_name].transpose()
        # add new columns into df
        dict_of_pos_control_dfs[sample_name].insert(0, 'passing', pos_control_check)
        dict_of_pos_control_dfs[sample_name].insert(2, 'qc ranges', pos_control_ranges)

    return(dict_of_pos_control_dfs)


# how to parse dictionary of dfs example
    # for name, df in dict_of_pos_control_dfs.items():
    #     pos_stats_list = df[0].tolist()
        # print(pos_stats_list)

# function to format negative control info
    # only need filtered HIV reads and Contig coverage
def get_neg_control_stats(negative_control_file):
    for file in negative_control_file:
        df=pd.read_csv(file, delimiter= "\t")
        df = df.drop(columns=(['Total Raw Read1', 'Total HIV Read1(% of raw)']))

        row_list = df.loc[0, 'Total HIV Read1':].values.flatten().tolist()
        # print(row_list)

        #two spaces added for sample into rows
        neg_control_check = ["",""]
        neg_control_ranges = ["",""]

        # qc ranges
        filtered_hiv_reads_range = 1000 #lte
        contig_coverage_range = 0

        neg_control_ranges.append("<=" + str(filtered_hiv_reads_range))
        neg_control_ranges.append(contig_coverage_range)

        # check ranges
        if row_list[0] < filtered_hiv_reads_range:
            neg_control_check.append("PASSED")
        else:
            neg_control_check.append("FAILED")

        if row_list[1] == contig_coverage_range:
            neg_control_check.append("PASSED")
        else:
            neg_control_check.append("FAILED")

        df = df.transpose()
        df.insert(0, 'passing', neg_control_check)
        df.insert(2, 'qc ranges', neg_control_ranges)

        return(df)





#function to format each sample, need to get negative control stats from sample type here because negative samples dont have separate stats
def get_test_sample_stats(test_sample_files):
    dict_of_test_dfs = {}
    for file in test_sample_files:
        sample_name = str(file).rsplit("_", 1)[0]
        print("Sample name " + sample_name)
        dict_of_test_dfs[sample_name] = pd.read_csv(file, delimiter= "\t")

    # print(dict_of_test_dfs)

    build_df = pd.concat(dict_of_test_dfs.values(), ignore_index=True)
    # print(build_df)
    build_df = build_df.groupby(build_df["Sample"])

    final_df = build_df.first()

    print(final_df)

    average_insert_size = final_df['Average Insert Size (bp)'].mean()
    average_read_length = final_df['Average Read Length (bp)'].mean()

    return(final_df, average_insert_size, average_read_length)






def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():

    parser = argparse.ArgumentParser(description="Create a summary report using various statistics files generated from HIV Genopipe")
    parser.add_argument('InputFile', nargs='+', type=File, help='group of stats files')
    parser.add_argument('-o', '--output', type=str, help='<STR> Will be used as prefix for output')
    args = parser.parse_args()

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
        else:
            test_sample_files.append(file)

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

    # reconstruct for optional if no interop or control samples
    with open('final_report.tsv','a') as f:
        #need test samples first for interop input
        test_samples_df = get_test_sample_stats(test_sample_files)
        if not interop_file:
            pass
        else:
            interop_stats_df = get_interop_stats(interop_file, test_samples_df[1], test_samples_df[2] )
            f.write("Run Quality Control \n")
            f.write("QC Metric\t QC Status\t Value\t QC Range \n")
            interop_stats_df.to_csv(f, sep = "\t", header=False)
            f.write("\n\n")

        if not negative_control_file:
            pass
        else:
            neg_control_df = get_neg_control_stats(negative_control_file)
            f.write("Negative Control \n")
            f.write("QC Metric\t QC Status\t Value\t QC Range \n")
            neg_control_df.to_csv(f, sep = "\t", header=False)
            f.write("\n\n")
        if not pos_control_file:
            pass
        else:
            pos_control_dfs = get_pos_control_stats(pos_control_file) # is a dict of dfs
            f.write("Positive Control \n")
            for name, df in pos_control_dfs.items():
                df.to_csv(f, sep = "\t", header=False)
                f.write("\n")
            f.write("\n\n")
        f.write("Sample Quality \n")
        if metadata_file:
            # read csv meta data as new df and append to output test_samples_df
            metadata_df = pd.read_csv(metadata_file[0])
            metadata_df.rename(columns={"sample_id": "sampleId"}, inplace=True)
            joindf = test_samples_df[0].merge(metadata_df,on="sampleId", how="left")
            joindf.to_csv(f, sep = "\t",index=False)
        else:
            test_samples_df[0].to_csv(f, sep = "\t")



    # write new file, stitch all dfs together, remove header lines


if __name__ == "__main__":
    main()





        # with pd.ExcelWriter(summaryOut) as writer:
        # set_of_df[0].to_excel(writer, sheet_name='')
        # set_of_df[1].to_excel(writer, sheet_name='')
        # set_of_df[2].to_excel(writer, sheet_name='')
