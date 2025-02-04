#!/usr/bin/env python

#Author: Lynn Dotrang
#Email: Lynn.dotrang@nih.gov
#Date Dev start: 11-9-2023

import numpy as np
import sys
import pandas as pd

# Read in depth file

def stats_parse(stats_file, sample_name, sample_type):
    data=pd.read_csv(stats_file, delimiter='\t', header=None, names=['Metric','Value', 'Comment'])
    print(data)

    df_summary_stats = pd.DataFrame(columns=['Metric', 'Value'])
    df_summary_stats.set_index('Metric', inplace=True)

    for index, row in data.iterrows():
        if 'insert size average' in row['Metric']:
            average_insert = row['Value']
        elif 'average length' in row['Metric']:
            average_length = row['Value']
        elif 'error rate' in row['Metric']:
            error_rate = row['Value']
        elif 'bases mapped:' in row['Metric']:
            bases_mapped = row['Value']
        elif 'average quality:' in row['Metric']:
            average_quality = row['Value']

    # Info
    df_summary_stats.loc['Sample Name'] = sample_name
    df_summary_stats.loc['Sample Type'] = sample_type

    # Stats

    df_summary_stats.loc['Average Insert Size (bp)'] = average_insert
    df_summary_stats.loc['Average Read Length (bp)'] = average_length
    df_summary_stats.loc['Error Rate (Read Variation %)'] = error_rate
    df_summary_stats.loc['Bases Mapped Count'] = bases_mapped
    df_summary_stats.loc['Average Quality Score'] = average_quality

    print(df_summary_stats)

    # df_summary_stats.to_csv(sample_name + "_samtools_stats.csv", sep=",")
    return(df_summary_stats)

def depth_stats(depth_file, sample_name):
    data=pd.read_csv(depth_file, delimiter='\t', header=None, names=['Sequence','Position (bp)','Depth'])

    df_depth_stats = pd.DataFrame(columns=['Metric', 'Value'])
    df_depth_stats.set_index('Metric', inplace=True)

    # Median
    median_depth = data['Depth'].median()
    print("This is the median depth: ", median_depth)

    # Mean absolute deviation
    mad_depth = (data['Depth'] - data['Depth'].mean()).abs().mean()
    print("This is the new MAD depth: ", mad_depth)

    df_depth_stats.loc['Median Depth (reads)'] = median_depth
    df_depth_stats.loc['Depth MAD (reads)'] = mad_depth

    print(df_depth_stats)

    # df_depth_stats.to_csv(sample_name + "_depth_stats.csv", sep=",")
    return(df_depth_stats)

def main():
    # samtools stats file where grep ^SN | cut -f 2- has been performed already
    stats_file = open(sys.argv[1])
    depth_file = open(sys.argv[2])
    sample_name = sys.argv[3]
    sample_type = sys.argv[4]


    df_summary = stats_parse(stats_file, sample_name, sample_type)
    df_depth = depth_stats(depth_file, sample_name)

    df_all_data = pd.concat([df_summary, df_depth], axis=0)
    df_all_data.to_csv(sample_type + "_" + sample_name + "_combined_stats.csv", sep=",", header=False)

    print("Done!")

if __name__ == "__main__":
    main()
