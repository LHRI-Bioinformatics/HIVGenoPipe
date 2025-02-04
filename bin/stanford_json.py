#!/usr/bin/env python3

import sys
import os
from textwrap import TextWrapper
from datetime import datetime
from sierrapy import fastareader, SierraClient
from sierrapy.sierraclient import ResponseError
from gql import gql
import json
import re
import argparse

def get_seq(fastas):

    amb_files_dict={}


    for fasta in fastas:
        amb_percent=str(fasta).split(".")[-2]
        amb_files_dict[amb_percent]= fasta
        sample_name = str(fasta).rsplit('.', 1)[0]
        sample_name = str(sample_name).split('/')[-1]
        currSeq = {}

        with open(fasta, 'r') as f:
            seq = f.readlines()[1]
            currSeq['header']='mySeq'
            currSeq['sequence']=seq.rstrip("\n")
            # print(sequences)
            # sequences[sample_name]=seq.rstrip("\n")

        # print(amb_files_dict)
        # print(sequences)
        dir_path = os.path.dirname(os.path.realpath(__file__))
        client = SierraClient('https://hivdb.stanford.edu/graphql')
        # try:



        print("Getting characterization data from Stanford for {}".format(sample_name))
        data = client.execute(gql(open(os.path.join(dir_path,'sierra.gql')).read()),variable_values={"sequences": [currSeq]})
        for sequence in data['sequenceAnalysis']:
            sequence['inputSequence']['header']=sample_name

        # dfHIVSeqs['pipeline']=data['currentVersion']['text']

        # Serializing json
        # json_object = json.dumps(data, indent=4)

        json.dump(data, open((sample_name + ".json"),'w'),separators=(',', ': '), indent = 4)
        # with open(jsonFile, "w") as jsonOut:
        #     # Writing data to a file
        #     jsonOut.write(json_object)
        #     jsonOut.close()  # to change file access modes



def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():

    parser = argparse.ArgumentParser(description="Take fasta inputs and get Stanford drug report for each file.")
    parser.add_argument('-f','--fastas', type=File, nargs='+', help='<*.fasta> Set of fasta files from pysamstats parser (can include any number of Amb% and consensus files)')
    args = parser.parse_args()

    get_seq(args.fastas)




if __name__ == '__main__':
    main()
