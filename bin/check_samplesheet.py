#!/usr/bin/env python3

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/chipseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq_1,fastq_2,control
    SPT5_T0_REP1,SRR1822153_1.fastq.gz,SRR1822153_2.fastq.gz,negative
    SPT5_T0_REP2,SRR1822153_1.fastq.gz,SRR1822153_2.fastq.gz,positive
    SPT5_T0_REP3,SRR1822153_1.fastq.gz,SRR1822153_2.fastq.gz,test
    """

    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:

        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "fastq_1", "fastq_2", "sample_type"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, fastq_1, fastq_2, sample_type = lspl[: len(HEADER)]
            if sample.find(" ") != -1:
                print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
                sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Check sample_type columns have valid values
            ## make valid values negative, positive, test
            if sample_type:
                if sample_type.find(" ") != -1:
                    print(f"WARNING: Spaces have been replaced by underscores for sample_type: {sample_type}")
                    sample_type = sample_type.replace(" ", "_")
                if sample_type not in ('negative', 'positive', 'test'):
                    print_error(
                        "sample_type must be 'negative', 'positive', or 'test!",
                        "Line",
                        line,
                    )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2, antibody, sample_type]
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", fastq_1, fastq_2, sample_type]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", fastq_1, fastq_2, sample_type]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, antibody, sample_type ]]}
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(
                    [
                        "sample",
                        "single_end",
                        "fastq_1",
                        "fastq_2",
                        "sample_type",
                    ]
                )
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error(
                        f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                        "Sample",
                        sample,
                    )

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    control = val[-1]
                    # if control and control not in sample_mapping_dict.keys():
                    #     print_error(
                    #         f"Control identifier has to match does a provided sample identifier!",
                    #         "Control",
                    #         control,
                    #     )

                    fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
