#!/usr/bin/env python3

import os
# import pysam
import Bio
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import identity_match
from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.Align import PairwiseAligner
from Bio import SeqIO
from collections import OrderedDict
import numpy as np
import pandas as pd
import sys
import re
from os.path import exists
import argparse

def coverage_check(consensus_files, ref_file, from_pos, to_pos, coverage_only):
    ##updated 01/2025

    # consensus_file = SeqIO.read(consensus_file, "fasta")
    #get sample name from first input file, add to data frame
    sample_name = str(consensus_files[0]).rsplit('.')[0]
    sample_name = str(sample_name).split('/')[-1]
    print("Sample name: ", sample_name)

    df_contig_stats = pd.DataFrame(columns=['Metric', 'Value'])
    df_contig_stats.set_index('Metric', inplace=True)

    df_contig_stats.loc['Sample Name'] = sample_name

    # ref only need once
    ref_file = SeqIO.read(ref_file, "fasta")

    trimmed_ref = ref_file.seq[from_pos:to_pos]
    trimmed_ref_length = len(ref_file.seq[from_pos:to_pos])
    print("Ref length = ", trimmed_ref_length)

    #### TODO: for file in amb_files list, do alignment and return stats for each, make one report for all files ####
    # contig coverage will only be needed once
    # loop through all amb fastas and get alignment stats
    for fasta in consensus_files:

        # get amb% info
        amb_percent=str(fasta).split(".")[-2]
        amb_percent = amb_percent.replace('Amb', '')
        print(amb_percent)

        consensus_file = SeqIO.read(fasta, "fasta")
        # alignment = pairwise2.align.localms(consensus_file.seq, trimmed_ref, 1, 0, -1, -1)
        alignment = pairwise2.align.localms(consensus_file.seq, trimmed_ref, 1, 0, -1, -1, one_alignment_only=True)


        ## manual match tracker
        # match = []
        score = 0
        N_count = 0

        for a, b in zip(alignment[0][0],alignment[0][1]):
            if a == b:
                # match.append('|')
                score += 1
            else:
                # match.append(' ')
                if a == "N":
                    N_count += 1

        contig_coverage = 1 - (N_count/trimmed_ref_length)
        contig_error_rate = 1 - (score/trimmed_ref_length)

        df_contig_stats.loc['Contig Coverage (%)'] = contig_coverage

        error_rate_name = 'Contig Error Rate at ' + str(amb_percent) + '% (%)'
        df_contig_stats.loc[error_rate_name] = contig_error_rate

    print(df_contig_stats)

    if coverage_only == True:
        df_coverage_only = df_contig_stats.loc[['Sample Name', 'Contig Coverage (%)']]
        return sample_name, df_coverage_only
    else:
        return sample_name, df_contig_stats



        ## print manual version of first alignment
        # print(alignment_score[0][0])
        # print("".join(match))
        # print(alignment_score[0][1])


def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():
    parser = argparse.ArgumentParser(description="Local pairwise alignment to get contig coverage and error rate")
    # add nargs + to allow any number of entry files (many ambiguity)
    parser.add_argument('-q','--query-sequence', type=File, nargs='+', help='<.fasta> fasta file to compare')
    parser.add_argument('-r','--ref-sequence', type=File, help='<.fasta> reference fasta file to align to')
    parser.add_argument('-f','--from-position', type=int, nargs='?', default=0, help='<int> starting base position for the reference, leave blank to stat from beginning')
    parser.add_argument('-t','--to-position', type=int, nargs='?', default=-1, help='<int> ending base position for the reference, leave blank for end of seq')
    parser.add_argument('-c','--coverage_only', help='skip calculation of error rate', action='store_true')
    args = parser.parse_args()


    #order the input fastas
    print(args.query_sequence)
    sorted_list = sorted(args.query_sequence, key=lambda y: int(y.rsplit('.', 2)[1].replace('Amb', '')))
    print(sorted_list)

    #make ref seq upper
    with open (args.ref_sequence, 'r') as old_file, open ('UPPER.fasta','w') as new_file:
        new_file.write(old_file.read().upper())

    with open ('UPPER.fasta', 'r') as new_file:
        coverage_check_output = coverage_check(sorted_list, new_file, args.from_position, args.to_position, args.coverage_only)

    sample_name = coverage_check_output[0]
    output_df = coverage_check_output[1]
    output_df.to_csv(sample_name + "_contig_stats.csv", sep=",", header=False)


    print("Done!")

if __name__ == "__main__":
    main()
