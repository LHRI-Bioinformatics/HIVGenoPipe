#!/usr/bin/env python

#Author: Lynn Dotrang
#Email: Lynn.dotrang@nih.gov
#Date Dev start: 08-30-2023

import seaborn as sns
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Read in depth file

def plot_depth(depth_file, sample_name):
    data=pd.read_csv(depth_file, delimiter='\t', header=None, names=['Sequence','Position (bp)','Depth'])

    sns.lineplot(data=data,x='Position (bp)', y='Depth')

    plt.title(sample_name+" Depth")
    plt.savefig(sample_name+".depth.png", bbox_inches='tight', dpi=400)
    plt.close()

def main():
    depth_file = open(sys.argv[1])
    sample_name = sys.argv[2]
    plot_depth(depth_file, sample_name)
    print("Done!")

if __name__ == "__main__":
    main()
