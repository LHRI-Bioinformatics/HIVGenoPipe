#!/usr/bin/env python3

import os
import sys
import glob
import argparse
from pathlib import Path
from datetime import datetime

### Created by Nextflow nf-core community
### Edited by Lynn Dotrang 11-2024; public version

def parse_args(args=None):
    Description = "Generate HIVGenopipe samplesheet from a directory of FastQ files."
    Epilog = """Example usage: python fastq_dir_to_samplesheet.py <FASTQ_DIR> <OUT_DIR>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FASTQ_DIR", help="Folder containing raw FastQ files.")
    parser.add_argument(
        "-o",
        "--out_dir",
        default=".",
        dest="OUT_DIR",
        help="Output directory for samplesheet.")
    parser.add_argument(
        "-r1",
        "--read1_extension",
        type=str,
        dest="READ1_EXTENSION",
        default="_R1_001.fastq.gz",
        help="File extension for read 1.",
    )
    parser.add_argument(
        "-r2",
        "--read2_extension",
        type=str,
        dest="READ2_EXTENSION",
        default="_R2_001.fastq.gz",
        help="File extension for read 2.",
    )
    parser.add_argument(
        "-se",
        "--single_end",
        dest="SINGLE_END",
        action="store_true",
        help="Single-end information will be auto-detected but this option forces paired-end FastQ files to be treated as single-end so only read 1 information is included in the samplesheet.",
    )
    parser.add_argument(
        "-sn",
        "--sanitise_name",
        dest="SANITISE_NAME",
        default=True,
        action="store_true",
        help="Whether to further sanitise FastQ file name to get sample id. Used in conjunction with --sanitise_name_delimiter and --sanitise_name_index.",
    )
    parser.add_argument(
        "-sd",
        "--sanitise_name_delimiter",
        type=str,
        dest="SANITISE_NAME_DELIMITER",
        default="_",
        help="Delimiter to use to sanitise sample name.",
    )
    parser.add_argument(
        "-si",
        "--sanitise_name_index",
        type=int,
        dest="SANITISE_NAME_INDEX",
        default=2,
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    parser.add_argument(
        "-re",
        "--recursive",
        dest="RECURSIVE",
        action="store_true",
        help="Whether or not to search for FastQ files recursively in <FASTQ_DIR>.",
    )
    parser.add_argument(
        "-pc",
        "--positive_control_identifier",
        dest="POSITIVE_CONTROL_IDENTIFIER",
        default="Pos",
        help="String used to identify a positive control sample within name, sample_type will be marked as positive",
    )
    parser.add_argument(
        "-nc",
        "--negative_control_identifier",
        dest="NEGATIVE_CONTROL_IDENTIFIER",
        default="Neg",
        help="String used to identify a negative control sample within name, sample_type will be marked as negative",
    )
    return parser.parse_args(args)


def fastq_dir_to_samplesheet(
    fastq_dir,
    out_dir,
    sample_type="test",
    read1_extension="_R1_001.fastq.gz",
    read2_extension="_R2_001.fastq.gz",
    single_end=False,
    # default for our samples
    sanitise_name=True,
    sanitise_name_delimiter="_",
    sanitise_name_index=2,
    recursive=False,
    positive_control_identifier="Pos",
    negative_control_identifier="Neg"

):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[:sanitise_name_index]
            )
        return sample
    def get_sample_type(sample_name):
        """Retrieve sample type from name"""
        if sample.find(positive_control_identifier) != -1:
            sample_type="positive"
        elif sample.find(negative_control_identifier) != -1:
            sample_type="negative"
        else:
            sample_type="test"
        return sample_type

    def get_fastqs(extension, recursive=False):
        """
        Needs to be sorted to ensure R1 and R2 are in the same order
        when merging technical replicates. Glob is not guaranteed to produce
        sorted results.
        See also https://stackoverflow.com/questions/6773584/how-is-pythons-glob-glob-ordered
        """
        search_path = f"*{extension}"
        if recursive:
            search_path = f"**/*{extension}"
        return sorted(glob.glob(os.path.join(fastq_dir, search_path), recursive=recursive))

    read_dict = {}

    ## Get read 1 files
    for read1_file in get_fastqs(read1_extension, recursive):
        sample = sanitize_sample(read1_file, read1_extension)
        sample_type = get_sample_type(sample)
        if sample not in read_dict:
            read_dict[sample] = {"R1": [], "R2": [], "sample_type": []}
        read_dict[sample]["R1"].append(read1_file)
        read_dict[sample]["sample_type"].append(sample_type)

    ## Get read 2 files
    if not single_end:
        for read2_file in get_fastqs(read2_extension, recursive):
            sample = sanitize_sample(read2_file, read2_extension)
            read_dict[sample]["R2"].append(read2_file)

    ## Write to file
    timestamp = datetime.now().strftime("%Y%m%d")
    if len(read_dict) > 0:
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        #get name from fastqdir

        samplesheet_name = os.path.basename(os.path.normpath(fastq_dir))+"_" + timestamp + "_samplesheet.csv"

        output_path = os.path.join(out_dir,samplesheet_name)

        if os.path.isfile(output_path):
            print("file with name " + samplesheet_name + " already exists")
            sys.exit(1)
        else:
            with open(output_path, "w") as fout:
                header = ["sample", "fastq_1", "fastq_2", "sample_type"]
                fout.write(",".join(header) + "\n")
                for sample, reads in sorted(read_dict.items()):
                    for idx, read_1 in enumerate(reads["R1"]):
                        read_2 = ""
                        if idx < len(reads["R2"]):
                            read_2 = reads["R2"][idx]
                        sample_type = ""
                        sample_type = reads["sample_type"][idx]
                        sample_info = ",".join([sample, read_1, read_2, sample_type])
                        fout.write(f"{sample_info}\n")
    else:
        error_str = "\nWARNING: No FastQ files found so samplesheet has not been created!\n\n"
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastQ files\n"
        error_str += "  - '--read1_extension' parameter\n"
        error_str += "  - '--read2_extension' parameter\n"
        print(error_str)
        sys.exit(1)
    # use sheet to add metadata if option is selected
    return samplesheet_name



def main(args=None):

    args = parse_args(args)

    output_samplesheet_file = fastq_dir_to_samplesheet(
        fastq_dir=args.FASTQ_DIR,
        out_dir=args.OUT_DIR,
        read1_extension=args.READ1_EXTENSION,
        read2_extension=args.READ2_EXTENSION,
        single_end=args.SINGLE_END,
        sanitise_name=args.SANITISE_NAME,
        sanitise_name_delimiter=args.SANITISE_NAME_DELIMITER,
        sanitise_name_index=args.SANITISE_NAME_INDEX,
        recursive=args.RECURSIVE,
        positive_control_identifier=args.POSITIVE_CONTROL_IDENTIFIER,
        negative_control_identifier=args.NEGATIVE_CONTROL_IDENTIFIER
    )


if __name__ == "__main__":
    sys.exit(main())
