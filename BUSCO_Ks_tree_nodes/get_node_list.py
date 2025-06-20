#Script written by Josje Romeijn, May 28 2024

################################################################################
#------------------------SCRIPT DESCRIPTION-------------------------------------
#This script first goes over a tree to obtain a list of all its nodes.
#Next, per node it provides the two "children clades" (the next adjecent 2 clades,
#note: these can be leaves). For these two clades, it obtains the busco set
#these two children clades share, and for each gene in this busco set it retrieves
#(for each child clade) 1 species that has the longest "version" of this busco gene.

#In the end, it returns a list in this format:
#"{node_name}__{species1}__{species2}__{busco_set}__{busco_gene}"
#e.g.:
#"node1__Zymtr1__Phaal1__fungi__100957at4751"
################################################################################

#-------------------------------STEP 1: GET DATA--------------------------------
#load modules
import argparse, os, glob, random, threading
import numpy as np
from ete3 import Tree
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

#get input from command
parser = argparse.ArgumentParser(description="Lists nodes and children")
parser.add_argument('-i','--input',help='Input tree',required=True)
parser.add_argument('-b','--busco',help='Directory of BUSCO output', required=True)

args=parser.parse_args()

#--------------------------STEP 2: DEFINE FUNCTIONS------------------------------

#define busco sets for order, class and phylum level
orders = ["agaricales", "boletales", "capnodiales","chaetothyriales",
            "eurotiales","glomerellales","helotiales","hypocreales",
            "mucorales","onygenales","pleosporales","polyporales"]
classes = ["agaricomycetes","dothideomycetes","eurotiomycetes",
            "leotiomycetes", "saccharomycetes","sordariomycetes","tremellomycetes"]
phyla = ["ascomycota", "basidiomycota", "microsporidia", "mucoromycota"]

#function to get the BUSCO level of 2 species
def get_busco_level(species_a,species_b):
    busco_a = [i.replace(species_a + "_","") for i in os.listdir(busco_loc) if i.startswith(species_a)]
    busco_b = [i.replace(species_b + "_","") for i in os.listdir(busco_loc) if i.startswith(species_b)]
    overlap = [i for i in busco_a if i in busco_b and not i.endswith(".check") and not i.endswith(".fasta")]


    #if overlap = 1: fungi, 2: phylum, 3: class, 4: order
    if len(overlap) == 1:
        return "fungi"
    elif len(overlap) in range(2, 5):
        level = [d for d in orders if d in overlap] or \
                [d for d in classes if d in overlap] or \
                [d for d in phyla if d in overlap]
        return level[0] if level else "fungi"
    else:
        raise ValueError(f"Unexpected number of overlaps: {len(overlap)}")

        return busco








#function to get the busco sequences that overlap between 2 clades
#for each seq, it also returns 1 species from clade A and 1 from clade B with
#the LONGEST busco seq
def get_busco_seq_overlap(list_species_clade_a,list_species_clade_b, busco_level, node):
    #get list of fasta seqs of both species
    fasta_a_loc = [glob.glob(busco_loc + species_a + "_" + busco_level + "/run_" + busco_level + "_odb10/busco_sequences/single_copy_busco_sequences/*.faa")
                for species_a in list_species_clade_a]
    #this returns a nested list, make it a 1D list
    fasta_a_loc = [item for sublist in fasta_a_loc for item in sublist]
    #remove the file path, just get the busco seq (last part of file name)
    fasta_a = [i.split("/")[-1].split(".")[0] for i in fasta_a_loc]

    #do the same for clade b
    fasta_b_loc = [glob.glob(busco_loc + species_b + "_" + busco_level + "/run_" + busco_level + "_odb10/busco_sequences/single_copy_busco_sequences/*.faa")
                for species_b in list_species_clade_b]
    fasta_b_loc = [item for sublist in fasta_b_loc for item in sublist]
    fasta_b = [i.split("/")[-1].split(".")[0] for i in fasta_b_loc]


    #determine overlap
    overlap = set(fasta_a) & set(fasta_b)






    #now loop over the overlap sequences and only store 1 species per clade with
    #longest busco sequence.
    final_seqs = {"node":[],"species_a":[], "species_b":[], "busco_level":[], "seq":[]}
    for seq in overlap:
        #add busco seq to dict to store species that have longest busco seqs
        final_seqs["seq"].append(seq)
        final_seqs["node"].append(node)
        final_seqs["busco_level"].append(busco_level)

        #create temp dict to store all species that have this busco seq (and their lengths)
        existing_seqs = {"species_a":[], "loc_a":[], "length_a":[], "species_b":[], "loc_b":[], "length_b":[]}





        #loop over file_locs of species in clade A that have this busco seq
        for file_loc in fasta_a_loc:
            if seq in file_loc:
                for species in list_species_clade_a:
                    if species + "_" + busco_level in file_loc:
                        if seq == file_loc.split("/")[-1].split(".")[0]:
                            existing_seqs["species_a"].append(species)
                            existing_seqs["loc_a"].append(file_loc)
                            #get length of busco seq of this species and store it
                            fasta = SeqIO.read(file_loc, "fasta")
                            existing_seqs["length_a"].append(len(fasta.seq))
                            break #to exit inner for loop once species is found

        #again for clade B
        for file_loc in fasta_b_loc:
            if seq in file_loc:
                for species in list_species_clade_b:
                    if species + "_" + busco_level in file_loc:
                        if seq == file_loc.split("/")[-1].split(".")[0]:
                            existing_seqs["species_b"].append(species)
                            existing_seqs["loc_b"].append(file_loc)
                            #get length of busco seq of this species and store it
                            fasta = SeqIO.read(file_loc, "fasta")
                            existing_seqs["length_b"].append(len(fasta.seq))
                            break #to exit inner loop once species is found


        #find the longest busco seq
        #first: find the idx of where the longest busco seq is located in dict
        idx_longest_seq_a = [i for i, x in enumerate(existing_seqs["length_a"]) if x == max(existing_seqs["length_a"])]
        idx_longest_seq_b = [i for i, x in enumerate(existing_seqs["length_b"]) if x == max(existing_seqs["length_b"])]

        #first the case where there is one seq the longest
        if len(idx_longest_seq_a) == 1:
            final_seqs["species_a"].append(existing_seqs["species_a"][idx_longest_seq_a[0]])
        #now the case that multiple busco seqs are the longest (pick a random one)
        else:
            final_seqs["species_a"].append(existing_seqs["species_a"][random.choice(idx_longest_seq_a)])

        #repeat for clade b.
        if len(idx_longest_seq_b) == 1:
            final_seqs["species_b"].append(existing_seqs["species_b"][idx_longest_seq_b[0]])
        #now the case that multiple busco seqs are the longest (pick a random one)
        else:
            final_seqs["species_b"].append(existing_seqs["species_b"][random.choice(idx_longest_seq_b)])

    print(f"done with node {node}")
    return final_seqs


def process_node(node_name):
    #get children for each node
    node = tree.search_nodes(name=node_name)
    print(node)
    children = node[0].get_children()

    #for both children: get list of genomes in this clade)
    ch1 = [leaf.replace("_fungi","") for leaf in children[0].iter_leaf_names()]
    ch2 = [leaf.replace("_fungi","") for leaf in children[1].iter_leaf_names()]

    #get busco set for this node (do this by using the first species in both child_
    #clades)
    busco_level = get_busco_level(ch1[0], ch2[0])
    print(f"busco level is {busco_level}")

    #obtain a dictionary that per busco seq of this clade has stored the longest
    #species of clade A and clade B
    busco_seqs = get_busco_seq_overlap(ch1, ch2, busco_level, node_name)

    return busco_seqs



def write_to_file(node):
    result = process_node(node)
    with file_lock:
        with open("list_busco_seqs_species.txt","a") as f:
            for node, sp_a, sp_b, bus, seq in zip(result["node"],result["species_a"], result["species_b"], result["busco_level"], result["seq"]):
                f.write(f"{node}__{sp_a}__{sp_b}__{bus}__{seq}\n")









#load busco folder
busco_loc = args.busco

#load treefile
tree = Tree(args.input, format=1, quoted_node_names=True)


#-----------------STEP 3: FOR EACH NODE, GET LIST OF BUSCO SEQS-----------------

#get list of all nodes in the tree (which are NOT tips (=leaves))
nodes = [node.name for node in tree.traverse() if not node.is_leaf()]

#Initialize lock for multithreading
file_lock = Lock()


f = open("list_busco_seqs_species.txt","w")
f.close()

with ThreadPoolExecutor() as executor:
    [executor.submit(write_to_file, node) for node in nodes]

