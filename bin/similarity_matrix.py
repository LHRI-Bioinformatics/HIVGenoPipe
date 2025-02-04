#!/usr/bin/env python3
# Author: Lynn Dotrang
# Date: 02/2024
# Desc: performs pairwise alignment for group of files to generate similarity matrix

import os
from collections import OrderedDict
import numpy as np
import pandas as pd
import sys
import re
from os.path import exists
import argparse
import seaborn as sns
import matplotlib.pyplot as plt


def matrix(file):
    df=pd.read_table(file, delimiter= "\t", usecols=["readX", "readY", "no N similarity%"])
    idx = sorted(set(df['readX']).union(df['readY']))
    print(idx)

    df = df.pivot(index='readX', columns='readY', values='no N similarity%').reindex(index=idx, columns=idx)
    # https://stackoverflow.com/questions/77242993/how-to-reshape-pandas-dataframe-into-a-symmetric-matrix-corr-like-square-matrix
    df = df.combine_first(df.T)
    # .reindex(index=idx, columns=idx).fillna(0, downcast='infer')
    # .pipe(lambda x: x+x.values.T)
    print(df)

    file_name = str(file).split(".tsv")[0]

    plt.figure(figsize=(20,15))
    sns.heatmap(df, annot=True, fmt=".2f", square=True, vmin=95, vmax=100, cmap="Spectral").set(xlabel="", ylabel="")

    plt.savefig(file_name+"_matrix.png", bbox_inches='tight')


    return



def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():
    parser = argparse.ArgumentParser(description="Local pairwise alignment to get contig coverage and error rate")
    # add nargs + to allow any number of entry files (many ambiguity)
    parser.add_argument('-i','--input', type=File, nargs='+', help='<.tsv> set of one or more stats tables from similarity_score.py')
    args = parser.parse_args()

    for file in args.input:
        matrix(file)



    print("Done!")

if __name__ == "__main__":
    main()
