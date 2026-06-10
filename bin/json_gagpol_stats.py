#!/usr/bin/env python3

# use stanford json result to get coordinates for gag and pol
# need to refer back to the parsed TSV to figure out where it starts
# can we get the coordinates for the submitted fasta after we make it in the parsed script?

import json
import pandas as pd
import argparse
import os
import Bio
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import PairwiseAligner
from Bio import SeqIO
from Bio.Align import substitution_matrices

def get_gene_coordinates(jsonFile, gene):
    gene_name = str(gene)
    with open(jsonFile, "r") as jsonIn:
        data = json.load(jsonIn)

    # Initialize before loop so we can detect a missing gene cleanly
    gene_coordinates = None
    gene_seq = None
    num_na_start_offset = 0
    num_na_end_offset = 0

    for sequence in data['data']['sequenceAnalysis']:
        seqName = sequence['inputSequence']['header']
        print("Parsing characterization data for {}".format(seqName))
        for gene in sequence['alignedGeneSequences']:
            if gene['gene']['name'] == gene_name:
                print(gene['gene']['name'])
                print("First", gene_name, "Nuc pos: ", gene['firstNA'])
                print("Last", gene_name, "Nuc pos: ", gene['lastNA'])
                gene_coordinates = [gene['firstNA'], gene['lastNA']]
                print("Stanford coordinates", gene_coordinates)

                gene_seq = gene['alignedNAs']
                print(gene_name, "sequence", gene_seq)

                if gene['firstAA'] == 1 and gene['lastAA'] == gene['gene']['length']:
                    pass
                elif gene['firstAA'] != 1 or gene['lastAA'] != gene['gene']['length']:
                    print("***ALERT: UNEXPECTED ALIGNMENT***")
                    print("First AA: ", gene['firstAA'])
                    print("Last AA: ", gene['lastAA'])
                    print("Ref gene length: ", gene['gene']['length'])

                    first_aa = gene['firstAA']
                    num_aa_start_offset = first_aa - 1
                    num_na_start_offset = num_aa_start_offset * 3

                    last_aa = gene['lastAA']
                    if last_aa > gene['gene']['length']:
                        num_na_end_offset = 0
                    elif last_aa < gene['gene']['length']:
                        num_aa_end_offset = gene['gene']['length'] - last_aa
                        num_na_end_offset = num_aa_end_offset * 3

    if gene_coordinates is None:
        raise ValueError(
            f"Gene '{gene_name}' not found in {jsonFile}. "
            f"Sierra may not have aligned this gene for this sample."
        )

    return gene_coordinates, gene_seq, num_na_start_offset, num_na_end_offset

def get_stats_by_coordinates(df_pos, gene, start, end, amb, min_depth=10):
    df_gene = df_pos.truncate(before=start, after=end - 1)

    amb_percent = amb / 100
    amb_reads_threshold = min_depth / amb_percent

    below_min_depth = (df_gene['reads_all'] < min_depth).sum()
    percent_below_min_depth = below_min_depth / len(df_gene)
    below_amb_reads_threshold = (df_gene['reads_all'] < amb_reads_threshold).sum()
    percent_below_amb_depth = below_amb_reads_threshold / len(df_gene)
    print("Positions below min depth", min_depth, "in", gene, "is", below_min_depth)
    print("Positions below amb depth", amb_reads_threshold, "in", gene, "is", below_amb_reads_threshold)
    if below_min_depth > 0:
        print("****************SAMPLE FAILED IN", gene, "********************")
        print(below_min_depth)
    else:
        print(below_min_depth)
    print("Majority bases in gene", gene)
    print(df_gene['Max_Base'])

    return percent_below_amb_depth, below_min_depth, df_gene

def align_consensus_biopython(consensus_fa, gene_seq):

    gene_seq = gene_seq.replace("-", "")
    consensus = str(SeqIO.read(consensus_fa, "fasta").seq).upper()
    consensus = consensus.replace("-", "")
    gene_seq = str(gene_seq).upper()

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"

    aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")

    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    alignment = aligner.align(consensus, gene_seq)[0]

    consensus_blocks = alignment.aligned[0]
    start = int(consensus_blocks[0][0])
    end = int(consensus_blocks[-1][1])

    print(f"Consensus coordinates: {start}–{end}")
    return start, end

def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile + ' does not exist or is not a file.')
    return MyFile

def main():

    parser = argparse.ArgumentParser(description="Parses Interop File for stats")
    parser.add_argument('-p', '--pysamstats-parsed', type=File, help='<Sample_parsed.tsv> parsed stats file from pysamstats_parse.py')
    parser.add_argument('-s', '--samtools-consensus', type=File, help='<sample_samtools_consensus.fasta> intermediate fasta sequence file from samtools consensus')
    parser.add_argument('-j', '--json', type=File, nargs="+", help='<Sample_sierra.json> json file returned from sierra hivdb')
    parser.add_argument('-d', '--min-depth', type=int, metavar='', help='<INT> Minimum read depth permitted for base to be analyzed (default = 10)', default=10)
    parser.add_argument(
        '-o', '--output',
        default='gagpol_quality.csv',
        help='Output file name'
    )
    args = parser.parse_args()

    df_pos = pd.read_csv(args.pysamstats_parsed, sep="\t")

    sample_name = os.path.basename(args.pysamstats_parsed).rsplit("_", 1)[0]
    output_dict = {
        "Sample_name": sample_name,
    }

    list_of_ambs = []
    percent_below_min_depth = None  # ensure defined even if loop is skipped
    truncated_out = f"{sample_name}_gagpol_truncated.csv"

    for json_file in args.json:
        amb_value = int(json_file.split("Amb")[1].split("_")[0])
        print("Ambiguous value: ", amb_value)
        list_of_ambs.append(amb_value)

        # --- GAG ---
        try:
            gag_coords_json, gag_seq, gag_num_na_start_offset, gag_num_na_end_offset = get_gene_coordinates(json_file, "gag")
            gag_start, gag_end = align_consensus_biopython(args.samtools_consensus, gag_seq)
            gag_ok = True
        except ValueError as e:
            print(f"WARNING: {e}")
            output_dict["Error"] = "gag not found in Sierra JSON"
            gag_ok = False

        # --- POL ---
        try:
            pol_coords_json, pol_seq, pol_num_na_start_offset, pol_num_na_end_offset = get_gene_coordinates(json_file, "pol")
            pol_start, pol_end = align_consensus_biopython(args.samtools_consensus, pol_seq)
            pol_ok = True
        except ValueError as e:
            print(f"WARNING: {e}")
            output_dict["Error"] = "pol not found in Sierra JSON — truncated CSV written from gag region only"
            pol_ok = False

        # --- Determine region to use for stats ---
        if gag_ok and pol_ok:
            adjusted_gag_start = gag_start - gag_num_na_start_offset
            adjusted_pol_end = pol_end + pol_num_na_end_offset
            region_start = adjusted_gag_start
            region_end = adjusted_pol_end
            region_label = "gagpol"
        elif gag_ok and not pol_ok:
            # Fall back to gag-only region so Nextflow always gets its file
            print("WARNING: Falling back to gag-only coordinates for truncated CSV.")
            region_start = gag_start - gag_num_na_start_offset
            region_end = gag_end
            region_label = "gag_only"
        else:
            # Neither gene found — write an empty stub so Nextflow doesn't stall
            print("WARNING: Neither gag nor pol found. Writing empty truncated CSV stub.")
            pd.DataFrame(columns=df_pos.columns).to_csv(truncated_out, index=False)
            output_dict["Error"] = "gag and pol both missing from Sierra JSON"
            output_dict[f"Gagpol_Amb_{amb_value}_Quality"] = None
            continue

        percent_below_amb_depth, below_min_depth, df_gagpol = get_stats_by_coordinates(
            df_pos, region_label, region_start, region_end, amb_value, args.min_depth
        )

        metric_name = f"Gagpol_Amb_{amb_value}_Quality"
        output_dict[metric_name] = 100 * (1 - percent_below_amb_depth)

        # Always write truncated CSV — Nextflow needs this regardless
        df_gagpol.to_csv(truncated_out, index=False)
        print(f"Wrote truncated {region_label} CSV: {truncated_out}")

        if len(df_gagpol) < (region_end - region_start):
            output_dict["Error"] = output_dict.get("Error", "") + " | Unexpected sequence length"

        percent_below_min_depth = below_min_depth

    min_depth_name = f"Gagpol Positions below min depth ({args.min_depth})"
    output_dict[min_depth_name] = percent_below_min_depth

    list_of_ambs.sort()
    ordered_keys = (
        ['Sample_name'] +
        [f'Gagpol_Amb_{amb}_Quality' for amb in list_of_ambs] +
        [f'Gagpol Positions below min depth ({args.min_depth})']
    )

    print(output_dict)
    df_output = pd.DataFrame(
        [[output_dict.get(k) for k in ordered_keys]],
        columns=ordered_keys
    )
    df_output.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()