#!/usr/bin/env python3

import os
import pysam
import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from collections import OrderedDict
import sys
import re
from os.path import exists
import argparse

#122123 revamp LTD
# dont need this any more
def total_read1_count(bamfile):
    """Gets the total number of read1s in an unaligned (raw or cleaned) BAM file.

	Args:
		bamfile (str): The absolute path of the unaligned BAM files.

	Returns:
		dictionary: Contains a key value map to the total number of read1s.
	"""

    read1_stats = OrderedDict([('total_read1_count', 0)])
    try:
        bamfile = pysam.AlignmentFile(str(bamfile), 'rb', check_sq=False)
        for read in bamfile.fetch(until_eof=True):
            if read.is_read1:
                read1_stats['total_read1_count'] += 1
    except:
        print("error {}".format(IOError))
        pass

    return read1_stats

# function to get read count from trimmomatic summary
def get_raw_reads(trimmomatic_summary):
    read1_stats = OrderedDict([('total_read1_count', 0)])
    try:
        with open(trimmomatic_summary, 'r') as f:
            read_pairs = f.readline()
            print(read_pairs)
            read_pairs = read_pairs.split()
            print(read_pairs)
            read_pairs = int(read_pairs[-1])
            # used to be total reads so we would * 2 but we want pairs 04/04/24
            # raw_reads = read_pairs * 2
            print(read_pairs)
            read1_stats['total_read1_count'] = read_pairs
            print("Raw reads:", read1_stats)
    except:
        print("error {}".format(IOError))
        pass
    return read1_stats

# function to get read count from bb tools seal
def get_cleaned_reads(bbtools_seal_summary):
    read1_stats = OrderedDict([('total_read1_count', 0)])
    try:
        with open(bbtools_seal_summary, 'r') as f:
            stats_lines = f.readlines()
            matched_line = stats_lines[2]
            matched_line = matched_line.split()
            matched_reads = int(matched_line[1])
            # we need read pairs 04/04/24
            read1_stats['total_read1_count'] = matched_reads / 2
            print("HIV reads:", read1_stats)

    except:
        print("error {}".format(IOError))
        pass

    return read1_stats

def mapped_bam_stats(mapped_bam, min_read_len):
    """Gets portion of read1 statistics. Done from aligned (mapped) BAM files.

	Args:
		mapped_bam (str): The absolute path of the aligned BAM files.
        min_read_len (int): The threshold length of the reads that will be processed.

	Returns:
		dictionary: Contains a key value map to multiple statistics.
	"""

    min_read_len = int(min_read_len)
    total_read1 = 0

    read1_stats = OrderedDict([
        ("filtered_HIV_read1", 0),
        ("filtered_HIV_read1%", 0),
        ("duplicate_HIV_read1", 0),
        ("duplicate_HIV_read1%", 0),
        ("short_HIV_read1", 0),
        ("short_HIV_read1%", 0),
        ("non_properly_paired_HIV_read1", 0),
        ("non_properly_paired_HIV_read1%", 0),
    ])
    try:
        bamfile = pysam.AlignmentFile(str(mapped_bam), "rb")
        mate_reads = get_mate_reads(mapped_bam)
        for read in bamfile.fetch(until_eof = True):
            if read.is_read1:
                total_read1 += 1
                if read.query_name in mate_reads:
                    if read.is_duplicate or mate_reads[read.query_name].is_duplicate:
                        read1_stats["duplicate_HIV_read1"] += 1
                    elif len(read.query_sequence) < min_read_len or len(mate_reads[read.query_name].query_sequence) < min_read_len:
                        read1_stats["short_HIV_read1"] += 1
                    elif not read.is_proper_pair or not mate_reads[read.query_name].is_proper_pair:
                        read1_stats["non_properly_paired_HIV_read1"] += 1
                    else:
                        read1_stats["filtered_HIV_read1"] += 1
                else:
                    if read.is_duplicate:
                        read1_stats["duplicate_HIV_read1"] += 1
                    elif len(read.query_sequence) < min_read_len:
                        read1_stats["short_HIV_read1"] += 1
                    elif not read.is_proper_pair:
                        read1_stats["non_properly_paired_HIV_read1"] += 1
                    else:
                        read1_stats["filtered_HIV_read1"] += 1
        read1_stats["filtered_HIV_read1%"] = round((read1_stats["filtered_HIV_read1"] / total_read1) * 100, 6)
        read1_stats["duplicate_HIV_read1%"] = round((read1_stats["duplicate_HIV_read1"] / total_read1) * 100, 6)
        read1_stats["short_HIV_read1%"] = round((read1_stats["short_HIV_read1"] / total_read1) * 100, 6)
        read1_stats["non_properly_paired_HIV_read1%"] = round((read1_stats["non_properly_paired_HIV_read1"] / total_read1) * 100, 6)

    except:
        print("error {}".format(IOError))
        pass

    return read1_stats

def get_mate_reads(mapped_bam):
    """Gets the corresponding mate (read2) for each read1.

	Args:
		mapped_bam (str): The absolute path of the aligned BAM files.

	Returns:
		dictionary: Contains a key value map to the mate reads.
	"""
    bamfile = pysam.AlignmentFile(str(mapped_bam), "rb")
    mate_reads = {}
    for read in bamfile.fetch(until_eof = True):
        if read.is_read2:
            mate_reads[read.query_name] = read
    return mate_reads

def fasta_stats(fasta_file, flag):
    """Gets portion of read1 statistics. Done from the FASTA (consensus and
    ambiguous consensus) files.

	Args:
		fasta_file (str): The absolute path of the FASTA files.
        flag (str): The identifier of the type of FASTA file (_Unamb, _Amb5%, etc.)

	Returns:
		dictionary: Contains a key value map to the oberved statistics from the FASTA
        file
	"""
    read1_stats = OrderedDict([
        ("length" + flag, 0),
        ("N_count" + flag, 0),
        ("N_count%" + flag, 0),
    ])
    try:
        for seq_record in SeqIO.parse(str(fasta_file), "fasta"):
            read1_stats["length" + flag] = len(seq_record)
            record = str(seq_record.seq)
            read1_stats["N_count" + flag] = record.count("N", 0, len(seq_record))
            if read1_stats["N_count" + flag] != 0:
                read1_stats["N_count%" + flag] = round((read1_stats["N_count" + flag] / read1_stats["length" + flag]) * 100, 6)
            elif read1_stats["length" + flag] == 0:
                read1_stats["N_count%" + flag] == "NaN"
            else:
                read1_stats["N_count%" + flag] = 0
    except:
        print("error {}".format(IOError))
        pass

    return read1_stats

def fasta_stats(fasta_file):
    """Gets portion of read1 statistics. Done from the FASTA (consensus and
    ambiguous consensus) files.

	Args:
		fasta_file (str): The absolute path of the FASTA files.
        flag (str): The identifier of the type of FASTA file (_Unamb, _Amb5%, etc.)

	Returns:
		dictionary: Contains a key value map to the oberved statistics from the FASTA
        file
	"""
    read1_stats = OrderedDict([
        ("length", 0),
        ("N_count", 0),
        ("N_count%", 0),
        ("Sequence", "")
    ])
    try:
        for seq_record in SeqIO.parse(str(fasta_file), "fasta"):
            read1_stats["length"] = len(seq_record)
            record = str(seq_record.seq)
            read1_stats["Sequence"] = str(record)
            read1_stats["N_count"] = record.count("N", 0, len(seq_record))
            if read1_stats["N_count"] != 0:
                read1_stats["N_count%"] = round((read1_stats["N_count"] / read1_stats["length"]) * 100, 6)
            elif read1_stats["length"] == 0:
                read1_stats["N_count%"] == "NaN"
            else:
                read1_stats["N_count%"] = 0
    except:
        print("error {}".format(IOError))
        pass

    return read1_stats

def get_samtools_stats(samtools_stats_and_depth_file):
    """Gets portion of read1 statistics. Done from the samtools depth and samtools stats parsed csv file files."""
    # read1_stats = OrderedDict([
    #     ("Average Insert Size", 0),
    #     ("Average Read Length", 0),
    #     ("Error Rate", 0),
    #     ("Bases Mapped Count", 0),
    #     ("Average Quality Score", 0),
    #     ("Median Depth", 0),
    #     ("Depth MAD", 0)
    # ])
    read1_stats = OrderedDict()
    try:
        with open(samtools_stats_and_depth_file, 'r') as f:
            stats_pairs = f.readlines()
            for line in stats_pairs:
                line = line.split(",")
                read1_stats[line[0]] = line[1].rstrip("\n")


    except:
        print("error {}".format(IOError))
        pass
    return read1_stats

def get_contig_stats(contig_coverage_file):
    """Gets portion of read1 statistics. Done from the contig coverage csv file files."""

    read1_stats = OrderedDict()
    list_of_error_rates = []
    try:
        with open(contig_coverage_file, 'r') as f:
            stats_pairs = f.readlines()
            for line in stats_pairs:
                print("CONTIG STATS CHECK - lines of stats: ", line)
                line = line.split(",")
                read1_stats[line[0]] = line[1].rstrip("\n")

            for line in stats_pairs[2:]:
                line = line.split(",")
                amb_percent = str(line[0])
                amb_percent = amb_percent.split(" ")[-2]
                amb_percent = amb_percent.rstrip("%")

                list_of_error_rates.append(amb_percent)



    except:
        print("error {}".format(IOError))
        pass
    return read1_stats, list_of_error_rates


# def get_stats(project_name, sample_name, min_read_len, miseqPath):
def get_stats(trimmomatic_summary='', bbtools_sealstats='', aligned_bam=None, min_read_len=None, pysamstats_fastas=None, samtools_stats_and_depth_file='', contig_coverage_file=None, output_prefix= ''):
    """
	Creates a read1 statistics (.txt) file for each selected sample.
    Implements other methods which return certain portions of the statistics.

    Args:
        project_name (str): The name of the selected Miseq project folder.
        sample_name (str): The name of the selected sample.
        min_read_len (int): The threshold length of the reads that will be processed.
        project_path (pathlib.Path): The absolute path where the Miseq projects are located.
        stats_path (pathlib.Path): The absolute path where the statistics files are located.
        raw_path (pathlib.Path): The absolute path where the unaligned (raw) BAM files are located.
        cleaned_path (pathlib.Path): The absolute path where the filtered (cleaned) BAM files are located.
        refGuided_consensus_path (pathlib.Path): The absolute path where the refGuided_consensus files are located.
        mapped_path (pathlib.Path): The absolute path where the aligned (mapped) BAM files are located.
        amb_consensus_path (pathlib.Path): The absolute path where the ambiguous consensus files are located.

    Returns:
        void
    """
    read1_stats = OrderedDict()
    # read1_stats['Project'] = project_name

    ## Basic read stats
    print("This is trim file: ", trimmomatic_summary)
    if trimmomatic_summary == []:
        sample_name = "failed read"
        read1_stats['Total Raw Read1'] = 0
    else:
        sample_name = str(trimmomatic_summary).split('.', 1)[0]
        sample_name = str(sample_name).split('/')[-1]
        print("Getting stats for "+sample_name)

        if len(re.findall(r'P\d{6}', sample_name)) > 0:
            sampleId = re.findall(r'P\d{6}', sample_name)[0]
        else:
            sampleId = sample_name
        read1_stats['Sample'] = sample_name
        read1_stats['sampleId'] = sampleId

        rawReadPairs = (get_raw_reads(trimmomatic_summary))['total_read1_count']
        print("Main raw read pairs:", rawReadPairs)
        read1_stats['Total Raw Read1'] = rawReadPairs

    print("This is bbtools file: ", bbtools_sealstats)
    if bbtools_sealstats == []:
        read1_stats['Total HIV Read1'] = 0
    else:
        hivReadPairs = (get_cleaned_reads(bbtools_sealstats))['total_read1_count']
        print("Main hiv read pairs:", hivReadPairs)
        read1_stats['Total HIV Read1'] = hivReadPairs
        read1_stats['Total HIV Read1(% of raw)'] = round((hivReadPairs/rawReadPairs) * 100, 6)


    ## Bam stats
    if aligned_bam == [] or aligned_bam == None:
        pass
    else:
        mappedBamStats = mapped_bam_stats(aligned_bam, min_read_len)
        # print(mappedBamStats)

        # add to data
        read1_stats['Filtered HIV Read1'] = mappedBamStats["filtered_HIV_read1"]
        read1_stats['Filtered HIV Read1(% of HIV)'] = mappedBamStats["filtered_HIV_read1%"]
        read1_stats['Duplicate HIV Read1'] = mappedBamStats["duplicate_HIV_read1"]
        read1_stats['Duplicate HIV Read1(% of HIV)'] = mappedBamStats["duplicate_HIV_read1%"]
        read1_stats['Short HIV Read1'] = mappedBamStats["short_HIV_read1"]
        read1_stats['Short HIV Read1(% of HIV)'] = mappedBamStats["short_HIV_read1%"]
        read1_stats['Non-properly paired HIV Read1'] = mappedBamStats["non_properly_paired_HIV_read1"]
        read1_stats['Non-properly paired HIV Read1(% of HIV)'] = mappedBamStats["non_properly_paired_HIV_read1%"]

    ## Samtools stats
    print ("This is samtools file: ", samtools_stats_and_depth_file)
    if samtools_stats_and_depth_file == [] or samtools_stats_and_depth_file == None:
        pass
    else:
        # get new columns from samtools depth and stats output
        samtools_stats = get_samtools_stats(samtools_stats_and_depth_file)
        # add to data
            # add samtools stats stuff
        read1_stats['Average Insert Size (bp)'] = samtools_stats['Average Insert Size (bp)']
        read1_stats['Average Read Length (bp)'] = samtools_stats['Average Read Length (bp)']
        read1_stats['Error Rate (Read Variation %)'] = samtools_stats['Error Rate (Read Variation %)']
        read1_stats['Bases Mapped Count'] = samtools_stats['Bases Mapped Count']
        read1_stats['Average Quality Score'] = samtools_stats['Average Quality Score']
        read1_stats['Median Depth (reads)'] = samtools_stats['Median Depth (reads)']
        read1_stats['Depth MAD (reads)'] = samtools_stats['Depth MAD (reads)']
        read1_stats['Sample Type'] = samtools_stats['Sample Type']

    ## contig coverage
    if contig_coverage_file == [] or contig_coverage_file == None:
        #automatically covers any failing files such as negative controls
        read1_stats['Contig Coverage (%)'] = 0
    else:
        contig_coverage_stats = get_contig_stats(contig_coverage_file)

        coverage_dict = contig_coverage_stats[0]
        amb_percents = contig_coverage_stats[1]

        read1_stats['Contig Coverage (%)'] = coverage_dict['Contig Coverage (%)']
        if amb_percents == []:
            pass
        else:
            for i in amb_percents:
                read1_stats['Contig Error Rate at ' + i + '% (%)'] = coverage_dict['Contig Error Rate at ' + i + '% (%)']





    # Need to make into for loop to handle any number of fasta files coming from pysamstats parser step
    # Create dictionary that dynamically creates variables for each file amb type
    print("This is the pysam file input: ", pysamstats_fastas)
    if pysamstats_fastas == [] or pysamstats_fastas == None:
        pass
    else:
        amb_files_dict={}
        for fasta in pysamstats_fastas:
            amb_percent=str(fasta).split(".")[-2]
            amb_files_dict[amb_percent]= fasta
        print(amb_files_dict)

        for amb_percent,value in amb_files_dict.items():
            stat_name = str(amb_percent) + "Stats"
            stat_name = fasta_stats(value)
            read1_stats[amb_percent + ' Length'] = stat_name['length']
            read1_stats['Num' + amb_percent + ' N\'s'] = stat_name['N_count']
            read1_stats['%' + amb_percent + ' N\'s'] = stat_name['N_count%']
            read1_stats[amb_percent + ' Seq'] = stat_name['Sequence']

    if output_prefix == '' or output_prefix == None:
        file_name = sample_name
    else:
        file_name = output_prefix + "_" + sample_name

    print("Writing out file...", sample_name)
    with open((file_name + "_stats.txt"), "w+") as outfile:
        outfile.write("\t".join(["%s" % key for key in read1_stats.keys()]) + "\n")
        outfile.write("\t".join(["%s" % val for val in read1_stats.values()]) + "\n")

def File(MyFile):
    if not os.path.isfile(MyFile):
        # raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
        MyFile = []
    return MyFile

def main():
# 1 Trimmomatic summary
# 2 BBtools sealstats
# 3 Self-aligned bam
# 4 Min read length
# 5 Pysamstats group of Amb/consensus fastas
# sample name will be the prefix in NF
# get_stats should have 5 arguments then

    parser = argparse.ArgumentParser(description="Gather read statistics from raw reads, cleaned HIV reads, ambiguous consensus call fastas, and mapped bam from throughout pipeline.")
    parser.add_argument('-t','--trimmomatic-summary', type=File, help='<.summary> formatted file from output of \'trimmomatic\' process', nargs='?', const='')
    parser.add_argument('-b','--bbtools-sealstats', type=File, help='<.sealstats.txt> formatted file from output of \'BBtools seal\' process', nargs='?', const='')
    parser.add_argument('-s','--self-align-bam', type=File, help='<.bam> BAM output after self-alignment', nargs='?', const='')
    parser.add_argument('-l','--min-len', type=int, metavar='', help='<INT> Minimum read len permitted for read to be analyzed (default = 100)', default=100)
    parser.add_argument('-p','--pysamstats-fastas', type=File, nargs='*', help='<*.fasta> Set of fasta files from pysamstats parser (can include any number of Amb% and consensus files)', default=[])
    parser.add_argument('-d','--samtools-stats-and-depth', type=File, help='<.csv> output from samtools_stats_and_depth.py', nargs='?', const='')
    parser.add_argument('-c','--contig-coverage', type=File, help='<.csv> output from check_coverage.py', nargs='?', const='')
    parser.add_argument('-o','--output-prefix', type=str, help='Prefix for output file, suffix is sample name')
    args = parser.parse_args()

    get_stats(args.trimmomatic_summary, args.bbtools_sealstats, args.self_align_bam, int(args.min_len), args.pysamstats_fastas, args.samtools_stats_and_depth, args.contig_coverage, args.output_prefix)

if __name__ == "__main__":
    main()
