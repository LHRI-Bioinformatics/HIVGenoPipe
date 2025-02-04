#!/usr/bin/env python3

import json
import pandas as pd
from inspect import getsourcefile
import os
import sys
from os import stat
from pwd import getpwuid
import re
import os.path
import argparse

def parseSierraCharacterization(jsonFile):
    analysis_date=os.path.getmtime(jsonFile)
    # print("analysis_date: {}".format(analysis_date))
    # print("Reading characterization data from Stanford")
    with open(jsonFile, "r") as jsonIn:
        # Writing data to a file
        data = json.load(jsonIn)
        jsonIn.close()  # to change file access modes

    hivdb_version=data['currentVersion']['text']
    sierra_version=data['currentProgramVersion']['text']

    dfMutations=pd.DataFrame(columns = ['seq_id','protein','reference','mutation','mutation_site','codon','insertedNAs','isInsertion','isDeletion','isIndel','isAmbiguous','isApobecMutation','isApobecDRM','hasStop','isUnusual','mutation_type','remarks','sierra_version'])

    dfDRMutations=pd.DataFrame(columns = ['seq_id', 'drug_class', 'drug_name', 'drug_abbr', 'score', 'level', 'interpretation', 'mutations','hivdb_version'])

    dfProteinSequences=pd.DataFrame(columns = ['seq_id', 'protein', 'alignedNT','alignedAA', 'clade', 'matchPcnt','sierra_version'])

    for sequence in data['sequenceAnalysis']:
        seqName=sequence['inputSequence']['header']
        clade=sequence['bestMatchingSubtype']['displayWithoutDistance']
        print("Parsing characterization data for {}".format(seqName))
        for gene in sequence['alignedGeneSequences']:
            # print(gene)
            dfProteinSequences.loc[len(dfProteinSequences)] = [
                        seqName,
                        gene['gene']['name'],
                        gene['alignedNAs'],
                        gene['alignedAAs'],
                        clade,
                        gene['matchPcnt'],
                        sierra_version]

            for mutation in gene['mutations']:
                # print("{}\t{}({})\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                comments=[]
                for comment in mutation["comments"]:
                    comments.append(comment["text"])
                dfMutations.loc[len(dfMutations)] = [
                        seqName,
                        # "{}({})".format(gene['gene']['name'],gene['gene']['length']),
                        gene['gene']['name'],
                        # gene['firstAA'],
                        # gene['lastAA'],
                        # "{}{}{}".format(mutation['reference'],mutation['position'],mutation['AAs']),
                        mutation['reference'],
                        mutation['AAs'],
                        mutation['position'],
                        mutation['triplet'],
                        mutation["insertedNAs"],
                        mutation["isInsertion"],
                        mutation["isDeletion"],
                        mutation["isIndel"],
                        mutation["isAmbiguous"],
                        mutation["isApobecMutation"],
                        mutation["isApobecDRM"],
                        mutation["hasStop"],
                        mutation["isUnusual"],
                        mutation["primaryType"],
                        " ".join(comments),
                        sierra_version
                    ]

        # print(dfMutations.dtypes)
        for col in ['isInsertion','isDeletion','isIndel','isAmbiguous', 'isApobecMutation','isApobecDRM','hasStop','isUnusual']:
            dfMutations[col] = dfMutations[col].astype('bool')
        # print(dfMutations.dtypes)

        for dr in sequence['drugResistance']:
            for drugScore in dr['drugScores']:
                dr_mutations=[]
                for partialScore in drugScore['partialScores']:
                    for dr_mutation in partialScore['mutations']:
                        if dr_mutation['text'] not in dr_mutations:
                            dr_mutations.append(dr_mutation['text'])
                # print(", ".join(dr_mutations))

                dfDRMutations.loc[len(dfDRMutations)] = [
                        seqName,
                        drugScore['drugClass']['name'],
                        drugScore['drug']['fullName'],
                        drugScore['drug']['displayAbbr'],
                        drugScore['score'],
                        drugScore['level'],
                        drugScore['text'],
                        ", ".join(dr_mutations),
                        hivdb_version
                        ]

    return [dfProteinSequences,dfMutations,dfDRMutations]




def File(MyFile):
    if not os.path.isfile(MyFile):
        raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile

def main():

    parser = argparse.ArgumentParser(description="Take fasta inputs and get Stanford drug report for each file.")
    parser.add_argument('-j','--jsons', type=File, nargs='+', help='<*.json> Set of fasta files from Stanford json script (can include any number of Amb% and consensus files)')
    args = parser.parse_args()

    for json in args.jsons:
        sample_name = str(json).rsplit('.', 1)[0]
        sample_name = str(sample_name).split('/')[-1]
        # summaryOut = sample_name+".summary.xlsx"

        # print(parseSierraCharacterization(json))

        set_of_df = parseSierraCharacterization(json)
        # sheet_names = ["Protein_Sequences","Mutations","Drug_Resistance"]

        set_of_df[0].to_csv(sample_name+".Protein_Sequences.tsv", sep = "\t", index=False)
        set_of_df[1].to_csv(sample_name+".Mutations.tsv", sep = "\t", index=False)
        set_of_df[2].to_csv(sample_name+".Drug_Resistance.tsv", sep = "\t", index=False)

        # with pd.ExcelWriter(summaryOut) as writer:
        # set_of_df[0].to_excel(writer, sheet_name='')
        # set_of_df[1].to_excel(writer, sheet_name='')
        # set_of_df[2].to_excel(writer, sheet_name='')





if __name__ == '__main__':
    main()
