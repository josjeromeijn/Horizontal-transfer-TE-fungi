#script to carry out initial clustering of candidate HTTs
#by: Josje Romeijn, Oct '21

import argparse
import os, itertools
import pandas as pd
from collections import Counter
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from threading import Lock
import numpy as np

parser = argparse.ArgumentParser(description="Apply criterion 1")
parser.add_argument('-e','--edges',help='Output file to store edges between TE-pairs',required=True)
parser.add_argument('-bw','--blast_file_within',help='blast file with blast hits WITHIN clades',required=True)
parser.add_argument('-c','--clade_file',help='file that on each line has genome identifiers that are together in collapsed clade',required=True)
parser.add_argument('-ks','--ks_filtered_file',help='file with TE pairs after ks filtering - so TEs that are candidate HTTs',required=True)
parser.add_argument('-t','--threads', help='number of threads', required=True)
parser.add_argument('-n','--nodes',help='file with nodes', required=True)
parser.add_argument('-no','--nodes_output',help='file location for nodes output', required=True)
args = parser.parse_args()



def find_clade(value, clade_dict):
    for clade, objects in clade_dict.items():
        if value.split("__")[0] in objects:
            return clade
    return value.split("__")[0]

def get_combination(row, clade_dict):
    query = find_clade(row['q'], clade_dict)
    subject = find_clade(row['s'], clade_dict)
    return f"{min(query, subject)}--{max(query, subject)}"

def criterion_1(pair, can_htt, blast_within):
    """
    This function takes a batch of hits, the hits DataFrame, a pIDmatrix representing
    percentage identity between clades, and nHits as the total number of hits.

    It returns pairs of hits where the highest within-clade identity (intra) is
    greater than the between-clade identity (inter) of the hits, based on criterion 1.
    """

    #Step 1: create all possible pairs of hits between this pair of clades
    clade1, clade2 = pair.split("--")

    hits = can_htt.loc[lambda can_htt: can_htt['combination'] == pair]
    if len(hits) == 1:
        print(f"{pair} completed! Only 1 hit found.")
        return pd.DataFrame(columns=['hit1', 'hit2', 'status'])



    #filter blast_within file so that it only contains hits of these clades
    TE_COPIES = list(set(list(hits['q']) + list(hits['s'])))
    blast_filtered = blast_within[(blast_within['qseqid'].isin(TE_COPIES)) & (blast_within['sseqid'].isin(TE_COPIES))]



    grouped_blast = (blast_filtered.sort_values('length', ascending=False)
                 .drop_duplicates(subset=['qseqid', 'sseqid'])
                 .loc[:, ['qseqid', 'sseqid', 'pident', 'length']]
                 .reset_index(drop=True))


    pid_within = {(row['qseqid'], row['sseqid']): row['pident'] for _, row in grouped_blast.iterrows()}

    #create dataframe of all possible hit pairs
    hit_pairs = pd.DataFrame(list(itertools.combinations(hits['hit_id'], 2)), columns=['hit1', 'hit2'])

    #split hit IDs into TE columns, sorted for intra-clade comparisons
    hit_pairs[['h1_te1', 'h1_te2']] = hit_pairs['hit1'].str.split("___", expand=True).apply(np.sort, axis=1, result_type="expand")
    hit_pairs[['h2_te1', 'h2_te2']] = hit_pairs['hit2'].str.split("___", expand=True).apply(np.sort, axis=1, result_type="expand")

    #check if hits involve the same TEs
    hit_pairs['connected'] = (hit_pairs.apply(lambda row: (row['h1_te1'] in row['hit2']) or (row['h1_te2'] in row['hit2']), axis=1))

    connections = hit_pairs[hit_pairs['connected']].copy()
    connections = connections.assign(status='connected')
    hit_pairs = hit_pairs[~hit_pairs['connected']].copy()

    #retrieve pID values between and within genomes
    def retrieve_pid(te1, te2):
        return pid_within.get((te1, te2), pid_within.get((te2, te1), 0))

    hit_pairs['pID_within1'] = hit_pairs.apply(lambda x: float(retrieve_pid(x['h1_te1'], x['h2_te1'])), axis=1)
    hit_pairs['pID_within2'] = hit_pairs.apply(lambda x: float(retrieve_pid(x['h1_te2'], x['h2_te2'])), axis=1)

    hit_pairs = hit_pairs.merge(hits[['hit_id', 'pid']], left_on='hit1', right_on='hit_id').rename(columns={'pid': 'pID_between1'})
    hit_pairs = hit_pairs.merge(hits[['hit_id', 'pid']], left_on='hit2', right_on='hit_id').rename(columns={'pid': 'pID_between2'})

    #apply criterion 1
    #connect TE pairs that have higher within-species blast pID's than between-species blast pID's
    hit_pairs['criterion1'] = ((hit_pairs[['pID_within1', 'pID_within2']].max(axis=1) >
                                hit_pairs[['pID_between1', 'pID_between2']].max(axis=1)) &
                               ~hit_pairs['connected'])

    #collect results
    criterion1_connections = hit_pairs[hit_pairs['criterion1']].assign(status='criterion1')

    #combine results
    final_connections = pd.concat([connections[['hit1', 'hit2', 'status']],
                                   criterion1_connections[['hit1', 'hit2', 'status']]], ignore_index=True)

    #free up memory
    del hit_pairs
    del criterion1_connections
    del connections
    return final_connections

def write_to_file(pair, nodes):
    #get edges [connected or criterion1]
    if os.path.exists(args.edges + "_" + pair + ".csv"):
        pass
    else:
        result = criterion_1(pair, can_htt, blast_within)

        #write edges to file
        result.to_csv(args.edges + "_" + pair + ".csv", mode='w', header=["Target","Source","connection"], index=False, sep = '\t')
        del result

        #get nodes from node file and write to outputfile
        sp1, sp2 = pair.split("--")
        cl1 = clade_dict.get(sp1, [sp1])
        cl2 = clade_dict.get(sp2, [sp2])

        #get nodes that have a species of cl1 or cl2
        result = nodes[nodes['Node'].str.contains('|'.join([item + "__" for item in cl1])) & nodes['Node'].str.contains('|'.join([item + "__" for item in cl2]))]
        result.to_csv(args.nodes_output + pair + ".csv", mode = 'w', index=False, sep='\t')




#-------------------------------------------------------------------------------


#remove output file if it already exists
if os.path.exists(args.edges):
    os.remove(args.edges)


#read in candidate hits and blast file and nodes file
can_htt = pd.read_csv(args.ks_filtered_file, sep = '\t')
can_htt.columns = ['q','s','pid']

print("reading in candidate HTTs complete")
colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
blast_within = pd.concat([chunk for chunk in tqdm(pd.read_csv(args.blast_file_within, chunksize=1000, sep = '\t', names=colnames, header=None), \
                                                  total=49073, desc='Loading data')])
print("reading in within blast hits complete")

nodes = pd.read_csv(args.nodes, sep = '\t')

print("reading in nodes file complete")

uniq_tes = list(set(list(can_htt['q']) + list(can_htt['s'])))
print(f"number of unique TEs is {len(uniq_tes)}")

#STEP 2: determine all genome-genome combinations that have candidate HTTs
#load in collapsed nodes file
clade_dict = {}
with open(args.clade_file, 'r') as file:
    for index, line in enumerate(file, start =1):
        clade_dict[f"clade{index}"] = line.strip().split()


#make column with the combination of clades, add to candidate hits dataframe
tqdm.pandas()
can_htt['combination'] = can_htt.progress_apply(get_combination, axis=1, clade_dict=clade_dict)

#make column that will become hit_id
can_htt['hit_id'] = can_htt.apply(lambda row: str(row['q']) + "___" + str(row['s']), axis=1)

print(f"number of clade-clade comparisons is {len(set(can_htt['combination']))}")



#STEP 3: PERFORM CRITERION 1
print("Now perform criterion 1")

freq = Counter(can_htt['combination'])



with ProcessPoolExecutor(int(args.threads)) as executor:
    futures = {executor.submit(write_to_file, pair, nodes): pair for pair, _ in freq.items()}

    for future in as_completed(futures):
        pair = futures[future]
        try:
            result = future.result()
            print(f"{pair} completed!")
        except Exception as e:
            print(f"{pair} generated an issue: {e}")
