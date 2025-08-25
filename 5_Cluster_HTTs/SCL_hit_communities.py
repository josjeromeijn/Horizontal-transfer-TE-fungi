#single linkage clustering of hit communities
#Josje Romeijn, 13-12-2024
import os, itertools, argparse
import pandas as pd
from tqdm import tqdm
import networkx as nx

#1------------------------GET INPUT FROM COMMAND--------------------------------
parser = argparse.ArgumentParser(description="Single linkage clustering of hit communities")
parser.add_argument('-i','--input',help='Directory containing hit community files with TE information',required=True)
parser.add_argument('-o','--output',help='Output directory ', required=True)

args=parser.parse_args()

def single_linkage_clustering(G):
    clusters = list(nx.connected_components(G))
    return clusters

#2)-----------------------read in all community TE files------------------------
print("read in files")
genome_pairs = []
for file in os.listdir(args.input):
    if file.endswith("_TEs.csv"):
        genome_pair = pd.read_csv(args.input + file, sep = "\t")
        pair = os.path.basename(file).replace("comm_","").replace("_TEs.csv","")
        genome_pair["community"] = [pair + "__" + str(value) for value in list(genome_pair["community"])]
        genome_pairs.append(genome_pair)

print("concatenate files")
all_genome_pairs = pd.concat(genome_pairs, ignore_index=True)

#3)-----------create nodes file (all hit communities in all genome pairs)-------
with open(args.output + "nodes_hit_communities.csv", 'w') as nodes:
    nodes.write("ID\tsize_of_community\n")
    for comm in tqdm(set(all_genome_pairs["community"])):
        nodes.write(comm + "\t" + str((all_genome_pairs["community"] == comm).sum()) + "\n")

#4)---------------create edges and store in edges file--------------------------
#create empty dataframe to store the number of TEs that are shared between communities
shared_TEs = {}

#go over each TE to see in which communities it's present
print("find in which communities TEs occur")
for TE in tqdm(set(all_genome_pairs["TE"])):
    occurences = [value for value in all_genome_pairs[all_genome_pairs["TE"] == TE]["community"]]
    if len(occurences) > 1:
        for comb in list(itertools.combinations(occurences, 2)):
            if comb in shared_TEs.keys():
                shared_TEs[comb] += 1
            else:
                shared_TEs[comb] = 1



#create edges file
print("create edges file")
with open(args.output + "edges_hit_communities.csv", 'w') as edges:
    edges.write("Source\tTarget\tnumber_shared_TEs\tpercentage_shared\n")
    for key in tqdm(shared_TEs.keys()):
        percentage_shared = int(shared_TEs[key]) / int(min((all_genome_pairs["community"] == key[0]).sum(), (all_genome_pairs["community"] == key[1]).sum()))
        edges.write(key[0] + "\t" + key[1] + "\t" + str(shared_TEs[key]) + "\t" + str(percentage_shared) + "\n")

#5)---------------perform single linkage clustering-----------------------------

#load in node + edge file and create graph
edges_df = pd.read_csv(args.output + "edges_hit_communities.csv", sep = '\t')
nodes_df = pd.read_csv(args.output + "nodes_hit_communities.csv", sep = '\t')

g = nx.Graph()
g.add_nodes_from(nodes_df['ID'])
for _, row in edges_df.iterrows():
    g.add_edge(row['Source'], row['Target'])

#perform single linkage clustering
clusters = single_linkage_clustering(g)

#save clusters
with open(args.output + "clusters_hit_communities.csv", "w") as f:
    f.write("ID\tcluster\n")
    for i, cluster in enumerate(clusters, 1):
        for node in cluster:
            f.write(node + '\t' + str(i) + '\n')
