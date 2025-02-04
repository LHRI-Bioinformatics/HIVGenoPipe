try:
    import json
except ImportError:
    import simplejson as json
import sys
from deepdiff import DeepDiff
import os

__author__ = "bsherman"
__version__ = "1.0"

def Diff(li1, li2):
    return (list(set(li1) - set(li2)))

def getMutationList(jsonObj):
    aligned = jsonObj['viewer']['sequenceAnalysis'][0]['drugResistance']
    mutationList = []
    for segment in aligned:
        gene=segment['gene']
        genename=gene['name']
        mutationTypes=segment['mutationsByTypes']
        for mutationType in mutationTypes:
            type=mutationType['mutationType']
            if type != "Other":
                mutations=mutationType['mutations']
                for mutation in mutations:
                    name=mutation['text']
                    #position=mutation['position']
                    mutationList.append(type+"~"+genename+"~"+name)
                    #print(name+"\t"+genename+"\t"+type)
    return mutationList


def main():
    args = sys.argv[1:]


    a = json.load(open(args[0]))

    b = json.load(open(args[1]))

    mutations_a=getMutationList(a)
    mutations_b=getMutationList(b)

    sampleA=os.path.splitext(os.path.basename(args[0]))[0]
    sampleB=os.path.splitext(os.path.basename(args[1]))[0]

    print(sampleA+"\t"+sampleB+"\t"+str(Diff(mutations_a,mutations_b))+"\t"+str(Diff(mutations_b,mutations_a)))

if __name__ == "__main__":
    main()
