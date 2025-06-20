### Inigo Banales - 21.04.2021
### Combine Blastn HSPs that are distant by less than maxGAP and overlap by less than maxOVERLAP

# Set directories and input
import argparse,os
parser = argparse.ArgumentParser(description="Combine Blastn HSPs that are distant by less than maxDist and overlap by less than maxOverlap")
parser.add_argument('-i','--input',help='input file with the Blastn matches to combine HSPs ',required=True)
parser.add_argument('-g','--guide',help='(optional) txt file with the Blastn txt files to filter ',required=False)
parser.add_argument('-o','--out',help='output folder',required=True)
parser.add_argument('-mg','--maxGAP',help='maxGAP value',required=True)
parser.add_argument('-mo','--maxOVERLAP',help='maxOVERLAP value',required=True)


args=parser.parse_args()

file = args.input
output_dir = args.out

# Import libraries
import re, pandas as pd
import statistics as st
import math as mt
import numpy as np
import bottleneck
import numexpr
from itertools import groupby
from operator import itemgetter



# FUNCTIONS 
def removeNestedHits(dataframe):

    # Open a list to store indexes of rows to drop because they are nested within any other row
    drop = []

    # Check for every row if the alignment in the query or in the subject is "nested" within any other
    for index, row in dataframe.iterrows():
        qstart = getattr(row,"qstart")
        qend = getattr(row,"qend")
        sstart = getattr(row,"sstart")
        send = getattr(row,"send")

        # Compare this row against all other rows to find conflicts
        for i in range(0, len(dataframe)):
            if (qstart > dataframe["qstart"][i] and qend < dataframe["qend"][i]) or (sstart > dataframe["sstart"][i] and send < dataframe["send"][i]):
                  drop.append(index)

    dataframe.drop(drop, inplace = True)
    dataframe.reset_index(drop=True,inplace=True)
    return dataframe


def combineHSPs (dataframe, maxGAP, maxOVERLAP):

    #calculate gap between query sequences
    qdistance = (dataframe["qstart"].shift(-1) - dataframe["qend"])[0:len(dataframe)-1]
    sdistance = (dataframe["sstart"].shift(-1) - dataframe["send"])[0:len(dataframe)-1]

    # See which distances meet our criteria
    which = np.where((qdistance < maxGAP) & (sdistance < maxGAP) & (qdistance > -maxOVERLAP) & (sdistance > -maxOVERLAP))[0]

    # If no distances meet criteria, or if different "groups" of HSPs could be merged
    # (so distances are not consecutive -very very rare case-) return a dataframe only with the top HSP
    if len(which) == 0:
        dataframe.sort_values(["length","bitscore","evalue"], inplace=True, ascending = False, ignore_index = True)
        dataframe.drop(range(1,len(dataframe)),inplace=True)

    # If one or more distances meet criteria, try to merge HSPs
    else:
       #find consecutive numbers in list of distances that can be matched (so
        which_con = [[group[1] for group in group] for _, group in groupby(enumerate(which), lambda x: x[0] - x[1])]

        #for each group of consecutive numbers, merge the HSP's.
        for i in which_con:
            merge = list(range(i[0],i[-1]+2))

            newhit = {
            "query": list(set(dataframe[dataframe.index.isin(merge)]["query"]))[0],
            "subject": list(set(dataframe[dataframe.index.isin(merge)]["subject"]))[0],
            "pident": st.mean(dataframe[dataframe.index.isin(merge)]["pident"]),
            "length": max(dataframe[dataframe.index.isin(merge)]["qend"]) - min(dataframe[dataframe.index.isin(merge)]["qstart"]),
            "mismatch": sum(dataframe[dataframe.index.isin(merge)]["mismatch"]),
            "gapopen": sum(dataframe[dataframe.index.isin(merge)]["gapopen"]),
            "qstart": min(dataframe[dataframe.index.isin(merge)]["qstart"]),
            "qend": max(dataframe[dataframe.index.isin(merge)]["qend"]),
            "sstart": min(min(subdf[subdf.index.isin(merge)]["sstart"]), min(subdf[subdf.index.isin(merge)]["send"])),
            "send": max(max(subdf[subdf.index.isin(merge)]["sstart"]), max(subdf[subdf.index.isin(merge)]["send"])),
            "evalue": mt.prod(dataframe[dataframe.index.isin(merge)]["evalue"]),
            "bitscore": sum(dataframe[dataframe.index.isin(merge)]["bitscore"]),
            }

            dataframe.drop(i, inplace = True)
            dataframe.loc[len(dataframe)] = newhit

    #in case of multiple HSP being kept (one or multiple large gaps prohibit HSPs from being merged), keep longest HSP
    dataframe = dataframe.loc[[dataframe["length"].idxmax()]]


    #reorder the dataframe based on group, qstart and sstart
    dataframe.sort_values(["qstart", "sstart"], inplace=True, ascending = True, ignore_index = True)

    return dataframe

# Option to tell Pandas we are modifying copies of other DataFrames so that it doesn't give warnings
pd.options.mode.chained_assignment = None

# List with blast-tab column names
fields = ["query", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

# Functions required later
blastn = pd.read_table(file,
                   sep = "\t",  #Setting  separator allows for using "c" engine
                   engine= "c", #"c" engine is faster than "python"
                   names = fields)


# Open new data frame!
newblastn = pd.DataFrame(data = None,
                         columns = fields)

# 1) Add column with number of directions (so hit is forward-forward: 1 or forward-reverse: 2).
blastn["#directions"] = [1 if sstart > send else 2 for sstart, send in zip(blastn["sstart"], blastn["send"])]

# 2) Add a new column with querysubject as "pair"
blastn["pair"] = blastn["query"] + blastn["subject"] + str(blastn["#directions"])

# 3) Add a new Numerical column based on the query-subject pair group ("group")
blastn["group"] = pd.factorize(blastn["pair"])[0]

# 4) Order rows by group, qstart and sstart
blastn.sort_values(["group","qstart","sstart"], inplace=True, ascending = True, ignore_index = True)

# 5) Delete "pair" column
blastn.drop("pair", axis=1, inplace=True)

# 6) Choose the rows that do not belong to a group with multiple HSPs (groups which only have one row) and adds this rows to the new data frame
onehspgroups = np.unique(blastn.group)[np.unique(blastn.group,return_counts = True)[1] == 1].tolist()
onehsprows = blastn[blastn["group"].isin(onehspgroups)]

newblastn = pd.concat([newblastn,onehsprows])

newblastn.reset_index(drop=True, inplace=True)
newblastn.drop("group", axis=1, inplace=True)

# 7) Procceed to try to merge HSPs for the groups with more than one "HSPs"
multiplehspgroups = np.unique(blastn.group)[np.unique(blastn.group,return_counts = True)[1] > 1].tolist()

    # Loop through each group
for j in list(set(multiplehspgroups)):

    # Subselect rows for each group
    subdf = blastn[blastn["group"]==j]
    subdf.reset_index(drop=True,inplace=True)
    subdf.drop("group", axis=1, inplace=True)

    # If there are more than 6 HSPs just skip this group and add the top HSP (it will increase performance in big files with many HSPs)
    if len(subdf) > 6:

        subdf.sort_values(["length","bitscore","evalue"], inplace=True, ascending = False, ignore_index = True)
        newblastn.loc[len(newblastn)] = subdf.loc[0]

    # Remove "nested" hits (lower quality because they are shorter)
    else:
        subdf = removeNestedHits(subdf)

        # If only one HSP is left after removing "nested" hits, add it to the new data frame
        if len(subdf) == 1:

            newblastn = pd.concat([newblastn,subdf])
            newblastn.reset_index(drop=True,inplace=True)

        # If 1-6 HSPs are left, merge them (or at least try...)
        elif len(subdf) > 1 and len(subdf) < 7:

            subdf = combineHSPs(subdf,int(args.maxGAP),int(args.maxOVERLAP))   # maxGAP = 300 and maxOVERLAP = 300

            newblastn = pd.concat([newblastn,subdf])
            newblastn.reset_index(drop=True,inplace=True)

#8) Export new blastn file
if not os.path.exists(output_dir.rsplit("/",1)[0]):
    os.makedirs(output_dir.rsplit("/",1)[0])
newblastn.to_csv(output_dir, index = None, header=None, float_format='%.2f',sep="\t")
