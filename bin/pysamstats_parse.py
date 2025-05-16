#!/usr/bin/env python

import sys
import argparse
import os.path
import pandas as pd

# Author: Lynn Dotrang
# 12/2024
# Added more dynamic consensus calling

def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

parser = argparse.ArgumentParser(description="Parse through the tsv output file from pysamstats to determine base calls based on selected ambiguity thresholds and minimum read depth. \n \
                                If no ambiguity threshold is declared, script will perform maximal consensus call. Use the -c flag to include the consensus FASTA along with ambiguous FASTA. ")
parser.add_argument('InputFile', type=File, help='<TSV> formatted file from output of \'pysamstats --type variation_strand\'')
parser.add_argument('OutputFileName', type=str, help='<STR> Will be used as prefix for output tsv and/or FASTA')
parser.add_argument('-a','--ambiguity', type=int, nargs='+', metavar='', help='<INT> Threshold for ambiguous base call as a percentage (if not provided, produce consensus only)')
parser.add_argument('-d','--min-depth', type=int, metavar='', help='<INT> Minimum read depth permitted for base to be analyzed (default = 200), positions under threshold will be marked \'N\'', default=10)
parser.add_argument('-i','--info-output', help='Include TSV output of positions with possible InDels and overall quality metrics', action='store_true')
parser.add_argument('-t','--tsv-output', help='Include TSV formatted output file with base call columns for every position', action='store_true')
parser.add_argument('-c','--consensus-output', help='Include FASTA formatted consensus sequence (maximally likely base call at each position)', action='store_true')
args = parser.parse_args()

if args.ambiguity is None:
    args.consensus_output = True

df_pysamstats = pd.read_csv(args.InputFile, sep="\t")

# Get columns of interest
df_pos = df_pysamstats[['pos', 'reads_all', 'deletions','insertions','A','T','C','G','N']].copy()
# print(df_pos)

df_pos["A_frac"] = df_pos['A'] / df_pos['reads_all'] * 100
df_pos["C_frac"] = df_pos['C'] / df_pos['reads_all'] * 100
df_pos["G_frac"] = df_pos['G'] / df_pos['reads_all'] * 100
df_pos["T_frac"] = df_pos['T'] / df_pos['reads_all'] * 100
df_pos["del_frac"] = df_pos['deletions'] / df_pos['reads_all'] * 100
df_pos["ins_frac"] = df_pos['insertions'] / df_pos['reads_all'] * 100
df_pos = df_pos.fillna(0)

# Make IUPAC dictionary remove Bio dependency
# IUPAC_dict = IUPACData.ambiguous_dna_values
IUPAC_dict = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "ACGT",
}

# Small function for getting the key from a value in dict
def get_key(val):
    for key, value in IUPAC_dict.items():
        if val == value:
            return key

# amb = integer ambiguity threshold
# min_read_depth = integer value minimum quality

def base_amb_call(amb, min_read_depth):
    # initiate some counts and lists
    column_name = "Amb" + str(amb) + "_base"
    passing_depths = []
    base_list = []
    indel_list = []


    for ind in df_pos.index:
        base_detected = False
        base_list.append("")
        if df_pos['reads_all'][ind] < min_read_depth:
            base_list[ind] = base_list[ind] + "N"

        else:
            # # keep list of depths over the min, they can still be under the overall read depth needed per ambiguity
            passing_depths.append(df_pos['reads_all'][ind])

            if df_pos['A_frac'][ind] > amb and df_pos['A'][ind] >= min_read_depth:
                base_list[ind] = base_list[ind] + "A"
                base_detected = True
            if df_pos['C_frac'][ind] > amb and df_pos['C'][ind] >= min_read_depth:
                base_list[ind] = base_list[ind] + "C"
                base_detected = True
            if df_pos['G_frac'][ind] > amb and df_pos['G'][ind] >= min_read_depth:
                base_list[ind] = base_list[ind] + "G"
                base_detected = True
            if df_pos['T_frac'][ind] > amb and df_pos['T'][ind] >= min_read_depth:
                base_list[ind] = base_list[ind] + "T"
                base_detected = True

            # If Del is over amb and if there is no other base over threshold then add deletion
            if df_pos['del_frac'][ind] > amb and base_detected is False:
                base_list[ind] = ""
                base_detected = True # so as to not proceed with next logic check

            # If none of the above conditions are met then just put N
            ## example: is pos depth = 12, A = 5, G = 7, base_detected should be false based on above logic
            if base_detected is False:
                base_list[ind] = base_list[ind] + "N"


        # Search for InDels and flag
        indel_list.append("")
        if df_pos['del_frac'][ind] > amb:
            indel_list[ind] = indel_list[ind] + "Del"
        elif df_pos['ins_frac'][ind] > amb:
            indel_list[ind] = indel_list[ind] + "Ins"


    # Add the new base list to the data frame
    df_pos[column_name] = base_list

    IUPAC_col_name = "Amb"+ str(amb) + "_IUPAC"
    IUPAC_list = []
    for base in base_list:
        if base != 'N':
            IUPAC_list.append(get_key(base))
        else:
            IUPAC_list.append(base)
    # Add IUPAC translated column to df
    df_pos[IUPAC_col_name] = IUPAC_list

    # save subset df
    amb_base_df = df_pos[["pos","reads_all",IUPAC_col_name]]
    print(amb_base_df)

    # Add InDel Flag column; new column per amb
    InDel_col_name = "Amb" + str(amb) + "_InDel_Flag"
    df_pos[InDel_col_name] = indel_list

    # Need for FASTAs write out
    return IUPAC_list, passing_depths, amb_base_df

# Separate function for MAX (consensus) call
def consensus_call(min_read_depth):
    max_list = []
    max_indel_list = []
    for ind in df_pos.index:
        # keep track of maximal base read
        max_value = 0
        max_list.append("")
        if df_pos['reads_all'][ind] < min_read_depth:
            max_list[ind] = max_list[ind] + "N"
        else:
            # if the base > 0, this will be the new max to append
            if df_pos['A_frac'][ind] > max_value:
                max_list[ind] = max_list[ind] + "A"
                max_value = df_pos['A_frac'][ind]
            elif df_pos['A_frac'][ind] == max_value:
                max_list[ind] = max_list[ind] + "A"

            if df_pos['C_frac'][ind] > max_value:
                max_list[ind] = "C"
                max_value = df_pos['C_frac'][ind]
            elif df_pos['C_frac'][ind] == max_value:
                max_list[ind] = max_list[ind] + "C"

            if df_pos['G_frac'][ind] > max_value:
                max_list[ind] = "G"
                max_value = df_pos['G_frac'][ind]
            elif df_pos['G_frac'][ind] == max_value:
                max_list[ind] = max_list[ind] + "G"

            if df_pos['T_frac'][ind] > max_value:
                max_list[ind] = "T"
                max_value = df_pos['T_frac'][ind]
            elif df_pos['T_frac'][ind] == max_value:
                max_list[ind] = max_list[ind] + "T"

            if df_pos['del_frac'][ind] > max_value:
                max_list[ind] = ""
                max_value = df_pos['del_frac'][ind]

            ## Do we need to do something here for insertion?
            # Search for InDels and flag
        max_indel_list.append("")
        if df_pos['del_frac'][ind] >= max_value:
            max_indel_list[ind] = max_indel_list[ind] + "Del"
        elif df_pos['ins_frac'][ind] >= max_value:
            max_indel_list[ind] = max_indel_list[ind] + "Ins"


    # In the event there are ties, IUPAC ambiguity is needed
    max_IUPAC_list =[ ]
    for base in max_list:
        if base != 'N':
            max_IUPAC_list.append(get_key(base))
        else:
            max_IUPAC_list.append(base)
    # Add IUPAC translated column to df
    df_pos['Max_Base'] = max_IUPAC_list
    # Add max indel flag to df
    df_pos["Max_InDel_Flag"] = max_indel_list




    # Need for FASTAs write out
    return max_IUPAC_list

def calculate_quality(df, depth_threshold):
    amb_column_name = df.columns[2]
    # Filter out rows where the amb_base_column is a deletion
    filtered_df = df[df[amb_column_name] != ""]
    # Find the first and last index where the amb_base_column value is not "N"
    start_index = filtered_df[filtered_df[amb_column_name] != 'N'].index.min()
    end_index = filtered_df[filtered_df[amb_column_name] != 'N'].index.max()
    # Trim the filtered DataFrame to include only the middle portion
    trimmed_df = filtered_df.loc[start_index:end_index]
    # Count rows where 'reads_all' >= depth_threshold
    count_reads_above_threshold = (trimmed_df['reads_all'] >= depth_threshold).sum()
    # Calculate the total number of rows in trimmed_df
    total_rows = len(trimmed_df)
    # Calculate the proportion
    quality = count_reads_above_threshold / total_rows if total_rows > 0 else 0

    quality = quality * 100

    return quality

# Perform function, make sure input params are numbers
# Grab the output list for FASTA
if args.ambiguity is None:
    pass
else:
    collect_quality_metrics = {
        "Ambiguity": [],
        "Read depth threshold": [],
        "Quality": []
        }
    for amb in args.ambiguity:
        # added passing check in function return
        IUPAC_fasta_list, passing_depths, amb_base_df = base_amb_call(int(amb), int(args.min_depth))
        # Make list into string, trim leading and trailing N
        IUPAC_string = ''
        for i in IUPAC_fasta_list:
            # Handles internal bases that may not pass thresholds
            if i is None:
                IUPAC_string += ""
            else:
                IUPAC_string += i



        consensus_length_untrimmed = len(IUPAC_string)
        IUPAC_string = IUPAC_string.lstrip('N').rstrip('N')
        consensus_length_trimmed = len(IUPAC_string)

        bases_trimmed = consensus_length_untrimmed - consensus_length_trimmed

        # check overall quality
        overall_read_depth = (args.min_depth/amb) * 100

        quality_value = calculate_quality(amb_base_df,overall_read_depth)




        # add to metrics
        collect_quality_metrics['Ambiguity'].append(amb)
        collect_quality_metrics['Read depth threshold'].append(overall_read_depth)
        collect_quality_metrics['Quality'].append(quality_value)

        if quality_value >= 90: # percentage needed to pass overall
            quality_pass = True
        else:
            quality_pass = False

        if quality_pass is True:
            with open(args.OutputFileName + ".Amb" + str(amb) + ".fasta", 'w') as f:
                f.write('>' + args.OutputFileName + ".Amb" + str(amb) + '_IUPAC'+'\n' + IUPAC_string + '\n')
        else: # add warning in file and fasta name
            with open("QUALITY_FAIL_" + args.OutputFileName + ".Amb" + str(amb) + ".fasta", 'w') as f:
                f.write('>' + args.OutputFileName + ".Amb" + str(amb) + '_IUPAC_QUALITY_FAIL'+'\n' + IUPAC_string + '\n')

if args.consensus_output is True:
    max_fasta_list = consensus_call(args.min_depth)
    # Make list into string, trim leading and trailing N
    max_string = ''
    for i in max_fasta_list:
        if i is None:
            max_string += ""
        else:
            max_string += i
    max_string = max_string.lstrip('N').rstrip('N')
    with open(args.OutputFileName + ".consensus.fasta", 'w') as j:
        j.write('>' + args.OutputFileName + '_consensus'+'\n' + max_string + '\n')

if args.tsv_output == True:
    # Output the df
    df_pos.to_csv(args.OutputFileName + "_parsed.tsv", sep="\t")
else:
    pass

if args.info_output == True:
    output_filename = args.OutputFileName + "_Report.tsv"
    df_indel = df_pos.copy()
    amb_values = args.ambiguity
    amb_values.sort()
    min_amb = str(amb_values[0])
    lowest_indel_col_name = "Amb" + min_amb + "_InDel_Flag"
    #Get flags for lowest amb to get all of them
    df_indel = df_indel.loc[(df_indel[lowest_indel_col_name] == "Del") | (df_indel[lowest_indel_col_name] == "Ins")]
    df_indel = df_indel.drop(columns=['deletions', 'insertions', 'A', 'T', 'C', 'G', 'N'])
    df_indel.to_csv(output_filename, sep="\t")
    # add quality stuff here
    with open(output_filename, 'a') as k:
        k.write('\n\n\n' + 'Quality Metrics' + '\n')
        for key, value in collect_quality_metrics.items():
            k.write('%s\t' % key)
            for i in value:
                k.write('%s\t' % i)
            k.write('\n')

else:
    pass
