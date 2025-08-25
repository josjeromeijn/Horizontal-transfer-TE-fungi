import leidenalg
import igraph as ig
import pandas as pd
from itertools import combinations
import argparse, re, os

def find(community, community_map):
    while community in community_map:
        community = community_map[community]
    return community

#1------------------------GET INPUT FROM COMMAND--------------------------------
parser = argparse.ArgumentParser(description="Creates list of unique nodes from multiple files containing edges")
parser.add_argument('-e','--edges',help='Directory of files with edges (files must end with .csv)',required=True)
parser.add_argument('-o','--output',help='Output file containing clusters', required=True)

args=parser.parse_args()

# Build Networkit graph
graph = ig.Graph.Read_Ncol(args.edges, directed=False)


print("Graph created")

# Apply Louvain clustering
clustering = leidenalg.find_partition(graph, leidenalg.ModularityVertexPartition)  #leidenalg


print("Clustering done")

#write output to file
with open(args.output, 'w') as f:
    for vertex_id in range(graph.vcount()):
        f.write(f"{vertex_id}\t{clustering.membership[vertex_id]}\n")

print("Done")
