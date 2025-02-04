#!/usr/bin/env python3

import pandas as pd
import argparse
import os

# from interop import py_interop_run_metrics, py_interop_run, py_interop_summary

def interop_info(input):
    RunDF = pd.DataFrame(columns=['Metric', 'Value'])
    RunDF.set_index('Metric', inplace=True)

    InteropSummary = pd.read_csv(input, sep=",", nrows=6, skiprows=2)
    print("New stats", InteropSummary,"\n\n")
    Q30 = InteropSummary.iloc[5]['%>=Q30']
    print(Q30)
    Yield = InteropSummary.iloc[5]['Yield'] * 1000000000
    print("Total yield: ", Yield)
    InteropDetails = pd.read_csv("interop.csv", sep=",", nrows=1, skiprows=12)
    print(InteropDetails)

    ClusterDensity = InteropDetails.iloc[0]['Density']
    ClusterDensity = ClusterDensity.split(' ', 1)[0]
    ClustersPF = InteropDetails.iloc[0]['Cluster PF']
    ClustersPF = ClustersPF.split(' ', 1)[0]
    # print("Cluster Density: {}\nClusters Passing Filter: {}\nQ30: {}".format(ClusterDensity, ClustersPF, Q30))
    # RunDF.loc['ClusterDensity'] = float(ClusterDensity)
    RunDF.loc['ClusterDensity'] = float(ClusterDensity)
    RunDF.loc['Q30'] = float(Q30)
    RunDF.loc['ClustersPF'] = float(ClustersPF)
    RunDF.loc['Yield'] = float(Yield)

    print(RunDF, "\n\n")

    RunDF.to_csv("interop_stats.csv", sep=",")

def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():

    parser = argparse.ArgumentParser(description="Parses Interop File for stats")
    parser.add_argument('-i','--interop-file', type=File, help='<interop.csv> interop file from Illumina MiSeq')
    args = parser.parse_args()


    interop_info(args.interop_file)





if __name__ == '__main__':
    main()
