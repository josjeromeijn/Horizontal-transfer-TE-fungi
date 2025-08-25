#written by: Josje Romeijn Nov '24

#script to find communities in a graph using greedy modularity maximization

import igraph as ig
import pandas as pd
from itertools import combinations
import argparse, re, os

#1------------------------GET INPUT FROM COMMAND--------------------------------
parser = argparse.ArgumentParser(description="Creates list of unique nodes from multiple files containing edges")
parser.add_argument('-e','--edges',help='Directory of files with edges (files must end with .csv)',required=True)
parser.add_argument('-o','--output',help='Output file containing clusters', required=True)
parser.add_argument('-n','--nodes',help='Directory of files with classifications', required=True)

args=parser.parse_args()



#function to find the final community of a given community
def find(community):
    while community in community_map:
        community = community_map[community]
    return community



#2-------------------------------LOAD IN FILES----------------------------------
#read  edge and node files
edges_df = pd.read_csv(args.edges, sep='\t') if os.path.getsize(args.edges) > 0 else pd.DataFrame(columns=["Target", "Source", "connection"])
#nodes_df = pd.read_csv(args.nodes, sep='\t')
nodes_df = pd.read_csv(args.nodes, sep='\t')

print("Done loading in files")

#create graph
graph = ig.Graph(directed=False)
graph.add_vertices(nodes_df['Node'].tolist())
graph.add_edges(list(zip(edges_df['Source'], edges_df['Target'])))

print("Graph created")

del edges_df, nodes_df

#3-------------------------------CLUSTERING-------------------------------------
g = graph.community_fastgreedy(weights=None).as_clustering()
print("Clustering done")


#4-------------CALCULATE WITHIN AND BETWEEN COMMUNITY CONNECTIVITY--------------
#get unique number of communities and get commnity assignments for each node
membership = g.membership
comm_idx = list(set(g.membership))

#calculate within cluster connectivity
within_connectivity = {}

for comm in comm_idx:
    #get nodes that are in this community
    nodes_in_comm = [i for i, m in enumerate(membership) if m == comm]
    if len(nodes_in_comm) == 1:
        within_connectivity[comm] = 1
    else:
        subgraph = graph.subgraph(nodes_in_comm)

        #within connectivity is number of realized edges divided by total edges
        within_connec = subgraph.ecount() / ( len(nodes_in_comm) * (len(nodes_in_comm) - 1) // 2 )

        #add to dict
        within_connectivity[comm] = within_connec

print("Within connectivity calculated")

#calculate between cluster connectivity
between_connectivity = {}

for comm1, comm2 in combinations(comm_idx, 2):
    nodes_in_comm1 = [i for i, m in enumerate(membership) if m == comm1]
    nodes_in_comm2 = [i for i, m in enumerate(membership) if m == comm2]

    #between connectivity is number of realized edges between communities
    #divided by total edges possible between communities.
    #first, create subgraph and create list of indices for subgraph
    subgraph = graph.subgraph(nodes_in_comm1 + nodes_in_comm2)
    subgraph_to_original = [graph.vs.find(name=subgraph.vs[node]["name"]).index for node in range(subgraph.vcount())]

    #count number of edges between comm1 and comm2
    edge_count = sum(1 for edge in subgraph.es \
        if (membership[subgraph_to_original[edge.source]] == comm1 and membership[subgraph_to_original[edge.target]] == comm2) or \
        (membership[subgraph_to_original[edge.source]] == comm2 and membership[subgraph_to_original[edge.target]] == comm1))

    #total number of possible edges between communities
    total_edges = len(nodes_in_comm1) * len(nodes_in_comm2)

    #calculate between connectivity of comm1 and comm2
    between_connec = edge_count / total_edges

    #add to dict
    between_connectivity[(comm1, comm2)] = [between_connec, edge_count, total_edges]

print("Between connectivity calculated")

#5-----------------COLLAPSE COMMUNITIES BASED ON CONNECTIVITY-------------------
#if between connectivity > 0.05, merge communities in one big community.

#get community combinations that show "high" connectivity
#to_collapse = [k for k, v in between_connectivity.items() if v >= 0.05]
community_map = {}

for (comm1, comm2), connectivity in between_connectivity.items():
    if connectivity[0] >= 0.05:
        #find the final communities
        final_comm1 = find(comm1)
        final_comm2 = find(comm2)

        #merge by mapping one community to the other
        if final_comm1 != final_comm2:
            community_map[final_comm1] = final_comm2

final_communities = {comm: find(comm) for comm in set(community_map.keys())}
print("Communities with high connectivity are checked")



#6---------------------------CREATE OUTPUT FILES--------------------------------

with open(args.output + ".csv", "w") as f:
    f.write("ID\tpid\tclass_o_s_f\tagreement_o_s_f\tclass_o_s\tagreement_o_s\tclass_o\tagreement_o\told_cluster\tnew_cluster\twithin_connec\n")
    with open(args.nodes, "r") as node_file:
        next(node_file)
        for line in node_file:
            node = line.split("\t")[0]
            line = line.replace("\n", "")
            #get index of node
            idx = graph.vs["name"].index(node)

            #get community
            old_comm = membership[idx]
            new_comm = final_communities.get(old_comm, old_comm)
            f.write(f"{line}\t{old_comm}\t{new_comm}\t{within_connectivity.get(old_comm)}\n")


with open(args.output + "_between_connectivity.csv", "w") as f:
    if between_connectivity:
        f.write("cluster1\tcluster2\tbetween_connec\tobserved_#_edges\tpossible_#_edges\n")
        for key in between_connectivity.keys():
            #write pair to file
            between_connec, edge_count, total_edges = between_connectivity.get(key)
            f.write(f"{key[0]}\t{key[1]}\t{between_connec}\t{edge_count}\t{total_edges}\n")




#6---------------------GET STATS OF COMMUNITIES---------------------------------
#per cluster, we want:
#   -   the number of different classifications (based on order,
#       superfam and family level) [TE level]
#   -   ^same, but on TE pair level
#   -   number of TEs involved
#   -   number of TE pairs involved
#   -   number of unknown TEs involved

#get nodes in graph
nodes = pd.read_csv(args.output + ".csv", sep = '\t')

#get the unique TEs
TE = [item for i in nodes["ID"] for item in i.split("___")]
#and their classifications
classif = [item for s in nodes["class_o_s_f"] for item in (
            s[:[i for i, char in enumerate(s) if char == "_"][4]], #split string on 5th underscore
            s[([i for i, char in enumerate(s) if char == "_"][4] + 3 ):] ) #retrieve 2nd classification from 7th underscore onwards
            ]


#and their community
comm = [item for item in nodes["new_cluster"] for _ in range(2)]

#make dataframe and drop duplicate rows
TEs = pd.DataFrame({"TE": TE, "classif_o_s_f": classif, "community": comm})
TEs = TEs.drop_duplicates()

TEs["classif_o_s"] = [item.rsplit("__",1)[0] for item in TEs["classif_o_s_f"]]
TEs["classif_o"] = [item.rsplit("__",1)[0] for item in TEs["classif_o_s"]]


#write TE information to output file
TEs.to_csv(args.output + "_TEs.csv", mode = "w", index=False, sep = '\t')


#make output file with communities and within community connectivity
with open(args.output + "_communities.csv", "w") as f:
    f.write("Community\tWithin_connectivity\tNumber_of_nodes\tunique_o_s_f\tunique_o_s\tunique_o\tunknowns\n")
    for comm in within_connectivity.keys():
        #how many nodes (=TE pairs) in community
        nodes_in_comm = len([i for i, m in enumerate(membership) if m == comm])

        #how many TEs in community
        TEs_in_comm = TEs[TEs["community"] == comm]

        #how many different classifications, without unknowns (TE level)
        unique_o_s_f = len({x for x in TEs_in_comm["classif_o_s_f"].tolist() if x != "Unknown"})
        unique_o_s = len({x for x in TEs_in_comm["classif_o_s"].tolist() if x != "Unknown"})
        unique_o = len({x for x in TEs_in_comm["classif_o"].tolist() if x != "Unknown"})

        #how many unknown TEs
        unknowns = len([x for x in TEs_in_comm["classif_o"].tolist() if x == "Unknown"])

        f.write(f"{comm}\t{within_connectivity[comm]}\t{nodes_in_comm}\t{unique_o_s_f}\t{unique_o_s}\t{unique_o}\t{unknowns}\n")
