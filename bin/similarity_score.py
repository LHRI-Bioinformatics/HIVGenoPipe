#!/usr/bin/env python3
# Author: Lynn Dotrang
# Date: 02/2024
# Desc: performs pairwise alignment for group of files to generate similarity matrix

import os
# import pysam
import Bio
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import PairwiseAligner
from Bio import SeqIO
from os.path import exists
import argparse
import itertools
import csv


def matcher(fasta1, fasta2):

    #Gather info pieces to append to report
    # readx, ready, length, identity count, ambiguity count, id%, sim%
    info = []

    # consensus_file = SeqIO.read(consensus_file, "fasta")
    #get sample name from first input file, add to data frame
    fasta1_name = str(fasta1).split('.fasta')[0]
    info.append(fasta1_name)
    fasta2_name = str(fasta2).split('.fasta')[0]
    info.append(fasta2_name)



    # sets for base comparison (from original supermatcher)
    nt=set(['A','T','G','C'])
    ambiguity_set=[]
    ambiguity_set.append(set(['A','W','R','M','D','V','H']))
    ambiguity_set.append(set(['T','W','Y','K','D','B','H']))
    ambiguity_set.append(set(['G','S','R','K','D','V','B']))
    ambiguity_set.append(set(['C','S','Y','M','B','V','H']))

    # read in fastas
    fasta_1 = SeqIO.read(fasta1, "fasta")
    fasta_2 = SeqIO.read(fasta2, "fasta")
    # pairwise2 deprecated, use Align.PairwiseAligner()

    # align fastas using similar scoring from supermatcher
    aligner = Align.PairwiseAligner(match_score=1.0, open_gap_score=-10, extend_gap_score=-0.5, mode="local")
    # aligner = Align.PairwiseAligner(mode="local")
    aln = aligner.align(fasta_1.seq, fasta_2.seq)
    # print("alignment score from aligner: ", aln[0].score)
    # match tracker
    identity_score = 0
    N_identity = 0
    ambiguity_count = 0
    manual_length = 0
    no_N_length = 0
    # print(alignment)
    for a, b in zip(aln[0][0],aln[0][1]):
        manual_length += 1
        no_N_length += 1
        if a == "N" or b == "N":
            no_N_length -= 1
            if a == "N" and b == "N":
                N_identity += 1
                identity_score += 1
        elif a == b:
            # match.append('|')
            identity_score += 1
        else:
            # logic for similarity from old supermatcher_pair_similarity
            pair=set([a,b])
            similarity_found = False
            # print("non matching pair: ", pair)
            # print("pair intersect with nt: ", pair.intersection(nt))
            # if len(pair.intersection(nt)) > 0:
            for k in range(0,4):
                # print("checking sets:", k)
                # print("pair intersect with amb set: ", pair.intersection(ambiguity_set[k]))
                # if both bases in the mismatch pair are in the same amb set, it will count
                # if the mismatch pair is possible in two sets it will count twice
                if len(pair.intersection(ambiguity_set[k])) > 1:
                    similarity_found = True
            if similarity_found == True:
                ambiguity_count = ambiguity_count +1
                # print("amb count: ", ambiguity_count)
                # test if this flag is faster than a break, it was
                # break

    aln_length= len(aln[0][0])
    info.append(aln_length)
    info.append(identity_score)
    info.append(ambiguity_count)

    identity_perc = 100*(identity_score/aln_length)
    info.append(identity_perc)

    similarity_perc = 100*(identity_score + ambiguity_count)/aln_length
    info.append(similarity_perc)

    info.append(manual_length)
    info.append(no_N_length)
    no_N_identity = identity_score - N_identity
    info.append(N_identity)
    info.append(no_N_identity)
    no_N_identity_perc = 100*(no_N_identity/no_N_length)
    info.append(no_N_identity_perc)
    no_N_similarity_perc = 100*(no_N_identity + ambiguity_count)/no_N_length
    info.append(no_N_similarity_perc)


    # print(aln[0][0])
    # print("Length: ", aln_length)
    # print("identity score = ", identity_score)
    # print("identity % = ", 100*(identity_score/aln_length))
    # print("similarity count = ", ambiguity_count)
    # print("similarity % = ", similarity)

    return info



def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():
    parser = argparse.ArgumentParser(description="Local pairwise alignment to get contig coverage and error rate")
    # add nargs + to allow any number of entry files (many ambiguity)
    parser.add_argument('-f','--fastas', type=File, nargs='+', help='<.fasta> set of fastas to compare (list individual files)')
    args = parser.parse_args()

    # get list of all fastas
    list_of_fastas = args.fastas

    # sort fastas by amb/consensus
    #make a dictory, loop through file to sort the files into dictionary keys
    amb_files_dict={}
    for fasta in list_of_fastas:
        amb_percent=str(fasta).split(".")[-2]
        if amb_percent in amb_files_dict:
            amb_files_dict[amb_percent].append(fasta)
        else:
            amb_files_dict[amb_percent]= [fasta]
    # print(amb_files_dict)

    header = ["readX", "readY", "length of alignment", "identity count","ambiguity count", "identity%", "similarity%", "manual length", "no_N_length", "N matches", "no N identity", "no N identity%", "no N similarity%"]
    # iterate all combinations and put in matcher func



    for amb_percent,values in amb_files_dict.items():
        #open a new file for each amb_percent
        with open((amb_percent + "_similarity_stats.tsv"), "w+") as outfile:
            tsv_output = csv.writer(outfile, delimiter='\t')
            tsv_output.writerow(header)
            all_fasta_combos = list(itertools.combinations_with_replacement(values,2))
            for combo in all_fasta_combos:
                stats_list = matcher(combo[0],combo[1])
                # write out stats to tsv
                tsv_output.writerow(stats_list)

                # print(stats_list)


    print("Done!")

if __name__ == "__main__":
    main()
