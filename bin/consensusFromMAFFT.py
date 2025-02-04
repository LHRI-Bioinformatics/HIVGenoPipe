#!/usr/bin/env python3

import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import Counter

def read_sequences_from_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences


def write_sequence_to_fasta(id,sequence,file_path,description=""):
    seq_record=SeqRecord(Seq(sequence), id=id, description=description)
    with open(file_path, 'w') as output_handle:
        SeqIO.write(seq_record, output_handle, 'fasta')

def create_hybrid_consensus(sequences):
    ########This assumes a multiple sequence alignment where the HIV reference is the first sequence
    if not all(len(s) == len(sequences[0]) for s in sequences):
        raise ValueError("All sequences must be of the same length")

    ## check reference length after mafft align
    reference = str(sequences[0])
    right_strip_ref = reference.rstrip("-")
    left_strip_ref = reference.lstrip("-")
    number_right_strip_bases = len(reference) - len(right_strip_ref)
    print("right stripped bases= ", number_right_strip_bases)
    number_left_strip_bases = len(reference) - len(left_strip_ref)
    print("left strip bases= ", number_left_strip_bases)


    hybrid_consensus = ''

    for i in range(len(sequences[0])):
        # Exclude '-' characters and count occurrences in the second through last sequences
        char_count = Counter(s[i] for s in sequences[1:] if s[i] != '-')

        if char_count:
            #Second fix
            # if there is a tie in the most common, add the ref and check again
            if len(char_count.most_common()) > 1:
                for letter, count in char_count.most_common():
                    if letter == sequences[0][i]:
                        # this automatically means that base will break the tie
                        majority_char = sequences[0][i]
            else:
                majority_char = char_count.most_common(1)[0][0]

        else:
            majority_char = sequences[0][i]

        hybrid_consensus += majority_char

    if len(hybrid_consensus) != len(sequences[0]):
        raise ValueError("Hybrid sequence must be of the same length before trim")

    if number_left_strip_bases != 0:
        hybrid_consensus = hybrid_consensus[number_left_strip_bases:]
    else:
        pass
    if number_right_strip_bases != 0:
        hybrid_consensus = hybrid_consensus[:-number_right_strip_bases]
    else:
        pass

    print(">final hybrid\n", hybrid_consensus)
    print("length of hybrid= ", len(hybrid_consensus))

    return hybrid_consensus

def main():
    # Path to your MAFFT alignment FASTA file
    alignment_file=sys.argv[1]

    # name for the sequence (i.e. "sample.mafft.hybrid_consensus")
    hybrid_header=sys.argv[2]

    #output fasta path and name (i.e. "output_path/sample.mafft.hybrid_consensus.fasta")
    hybrid_output_file=sys.argv[3]

    #alignment_file = '/array/Members/dotrangtl/0123_BS_LD_sample19_pos_seal_norm_mindepth100_repeat3/mafft/19-P245626_S3_T1.fas'
    #hybrid_header = '19-P245626_S3.mafft.hybrid_consensus'
    #hybrid_output_file = '/array/Members/dotrangtl/0123_BS_LD_sample19_pos_seal_norm_mindepth100_repeat3/mafft/19-P245626_S3_T1.mafft.hybrid_consensus.fasta'

    # Read alignment fasta into an array of string sequences
    aligned_sequences = read_sequences_from_fasta(alignment_file)

    # Pass in the array of string sequences and create the hybrid
    hybrid_consensus = create_hybrid_consensus(aligned_sequences)
    #print(hybrid_consensus)

    # Write the hybrid consensus sequence to a file
    write_sequence_to_fasta(hybrid_header, hybrid_consensus, hybrid_output_file)

if __name__ == "__main__":
    main()
