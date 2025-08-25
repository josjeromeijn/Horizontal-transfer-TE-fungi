#a script to count minimal number of HTTs from hit communities.
#By Josje Romeijn, 18-12-2024

from ete3 import Tree
from collections import defaultdict, Counter
import pandas as pd
import argparse
from tqdm import tqdm


parser = argparse.ArgumentParser(description="Single linkage clustering of hit communities")
parser.add_argument('-c','--clusters',help='File with clusters of hit communities',required=True)
parser.add_argument('-o','--output',help='Output directory', required=True)
parser.add_argument('-t', '--tree', help='Tree file', required=True)
parser.add_argument('-cl','--clade_file',help='file that on each line has genome identifiers that are together in collapsed clade',required=True)


args=parser.parse_args()


#-------------------------DEFINE FUNCTIONS--------------------------------------
def prune_tree_and_count_edges(treefile, tips, edges, clade_dict):
    #load tree
    tree = Tree(treefile, format=1)

    #prune tree to keep only the specified tips
    tree.prune(tips, preserve_branch_length=True)

    #convert edges to a lookup dictionary for quick access
    edge_counts = defaultdict(int)
    for a, b in edges:
        edge_counts[frozenset([a, b])] += 1

    #helper function to collect all tip names under a node
    def get_tips(node):
        return set(leaf.name for leaf in node.get_leaves())

    #counts variable to track valid edges
    counts = 0

    #traverse nodes from tips to root
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            #get children nodes
            child1, child2 = node.get_children()

            #get tip names for each child clade
            group1_tips = get_tips(child1)
            group2_tips = get_tips(child2)

            #check edges between group1_tips and group2_tips
            group_edges = set()
            for tip1 in group1_tips:
                for tip2 in group2_tips:
                    edge_key = frozenset([tip1, tip2])
                    if edge_key in edge_counts:
                        group_edges.add(edge_key)

            #add counts for unique edges
            if group_edges:
                counts += max([edge_counts[edge_key] for edge_key in group_edges])


    return counts

#1)-----------------------------LOAD FILES--------------------------------------
treefile = args.tree
clusters = pd.read_csv(args.clusters, sep = '\t')
#load in collapsed nodes file
clade_dict = {}
with open(args.clade_file, 'r') as file:
    for index, line in enumerate(file, start =1):
        clade_dict[f"clade{index}"] = line.strip().split()



#2)---------------COUNT HTT PER CLUSTER OF HIT COMMUNITIES----------------------
dupl = 0
min_htt = {}
for cluster in tqdm(set(clusters['cluster'])):
    #get list of all hit communities in this cluster
    involved_hc = [hc.split("__")[0] for hc in clusters[clusters['cluster'] == cluster]['ID']]
    edges = [(item.split("--")[0],item.split("--")[1]) for item in involved_hc]
    #since collapsed clades are involved, if "clade" is in tips, check which tips it involves and take the first tip
    edges = [tuple(clade_dict[tip][0] if tip.startswith("clade") else tip for tip in tup) for tup in edges]

    tips = {item for tup in edges for item in tup}

    #count separate hit communities of same pair in same cluster (shouldn't happen)
    edge_counts = Counter(tuple(sorted(edge)) for edge in edges)
    duplicates = {edge for edge, count in edge_counts.items() if count > 1}
    dupl += len(duplicates)

    #get minimal number of HTTs in this cluster
    min_htt[cluster] = prune_tree_and_count_edges(treefile, tips, edges, clade_dict)


print(f"number of multiple hit communities of same clade pair in single cluster is {dupl}")

with open(args.output + "minimal_number_HTTs.txt", "w") as f:
    f.write("Cluster\tminimal_number_of_HTTs\n")
    for cluster in min_htt:
        f.write(str(cluster) + '\t' + str(min_htt[cluster]) + '\n')
