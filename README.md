# Horizontal transfer TE fungi
Scripts used to detect horizontal transfer of transposable elements (HTT) in 1,348 publicly available fungal genomes (JGI), see our [preprint](https://www.biorxiv.org/content/10.1101/2025.06.16.659975v1). Folders include Snakemake pipelines, scripts, and a file containing software versions used for each specific pipeline. All scripts included in these pipelines are created by Josje Romeijn (MGE lab, Utrecht University). 

1. Extraction of TEs from fungal genomes 
2. BUSCO Ks estimations of genome comparisons
3. TE vs TE blast
4. Calculation of Ks of TE pairs to obtain candidate HTTs 
5. Clustering of candidate HTTs and counting of minimal TE transfer events 
6. Retrieving HTT-associated TEs

## General software versions used for all pipelines: 
- Snakemake v8.12.0
- Python package pandas v2.2.3
- Python package tqdm v4.67.1

### 1. Extraction of TEs from fungal genomes
Make sure to only store .fasta files of fungal genomes included in this analysis in `genome_folder`, as the pipeline scans this directory for .fasta files for input. 
- BEDtools v2.31.1
- Transeq v6.6.0.0
- HMMER v3.4
- RPSTBLASTN v2.12.0+

### 2. BUSCO Ks estimations of genome comparisons
This pipeline starts by creating a file that identifies the best BUSCO gene comparisons for each comparison of 2 nodes in a phylogenetic tree (`list_busco_seqs_species.txt`). How it works is explained in the tree below. 

<img width="400" height="350" alt="image" src="https://github.com/user-attachments/assets/e550adc2-45b7-4a71-8936-166649a909c4" />


Let's take node 20 as an example. First, the most suitable BUSCO set is determined (priority: class-level, order-level, phylum-level, kingdom-level BUSCO set). In this case that's the kingdom-level BUSCO set, as the clades involves two different phyla. Next, for each BUSCO gene in the determined BUSCO set, the two longest copies of each BUSCO gene are searched in the two children clades (node 21 and node 22). This information is then added to a list (`list_busco_seqs_species.txt`). In case of equal lengths, a random genome is chosen from each child clade. This procedure is carried out for each node in the tree (20, 21, 22, 23, 24).

Format of `list_busco_seqs_species.txt`:

node_name__species1__species2__busco_set__busco_gene

e.g.:

node1__Zymtr1__Phaal1__fungi__100957at4751

Using this list, for each node, all the synonymous mutation rates of BUSCO genes are calculated. The amount of jobs in this pipeline are too much to handle for Snakemake (>4 million BUSCO gene comparisons x 3 jobs). Therefore, we've found a work-around by giving a subset of BUSCO comparisons to Snakemake using a for-loop in bash (`run_snakemake_in_batches.sh`). It's not pretty, but it works. 

Software versions:
- python package ete3 v.3.1.3
- mafft v7.505
- PAL2NAL v14
- R package seqinr v.4.2-36
- python package numpy v1.26.3
- python package biopython v1.78

After running this pipeline, the R-script `Ks_vs_ABL.R` (R v4.5.1, ggtree v4.4.3, patchwork v1.3.1, ggplot2 v3.5.2, treedataverse v0.0.1, stringr v1.5.1) can be used to obtain the 0.5% quantile of Ks distributions for each node in the tree. Note that this script is not part of the Snakefile. At the end, this script saves a file listing the node numbers (according to python `ete3`, the 0.5% quantiles, the nodeids (according to R `ggtree`), and all combinations possible within that treenode. An example of this file is included (`all_sp_sp_combinations.txt`). This file is needed for **4. Calculation of Ks of TE pairs to obtain candidate HTTs**.  

### 3. TE vs TE blast
Make sure to only store .fasta files of fungal genomes included in this analysis in `genome_folder`, as the pipeline scans this directory for .fasta files for input. Also, it is important to note for further analyses that we removed completely identical TE sequences per genome to reduce the workload. For later analyses, these "duplicate" TEs will have to be manually added. 

From the BUSCO comparisons, we estimated which nodes in the tree contain species that are too recently diverged to detect horizontal transfer (referred to as "collapsed nodes"). We've compiled a list of all species comparisons that are for that reason excluded from further analysis (`collapsed_sp_sp_combinations.txt`). TE blast hits involving TEs of a collapsed node are therefore filtered out, as well as "self"-hits (TE hits between TEs residing in the same genome). TE blast hits are then chained if they involve the same two TEs, their gap is <600 bp and overlap <600 bp. Next, only TE hits are retained where the alignment covers >60% of the length of both TEs. Finally, all TE hits are combined into a single file, after which only one copy of bidirectional TE hits is kept (only keep TE1_sp1 -- TE2_sp2 instead of also keeping TE2_sp2 -- TE1_sp1). This file is then split in multiple files containing max 4000 lines for further analysis. 

Software versions: 
- BLAST 2.12.0+
- python package collections v3.13.5
- python package numpy v1.26.3
  

### 4. Calculation of Ks of TE pairs to obtain candidate HTTs 
A custom python script (`calculate_Ks_TE_hit.py`) determines per TE hit whether there are any TE-related protein domains that overlap and determines the synonymous mutation rate of these TE-related protein domains (must be at minimum 100bp alignment). It takes into account that TE blast hits can be forward-forward or forward-reverse and makes sure TE-related protein domains are at the **exact** same location in the alignment and that the protein domains are in the right ORF compared to the alignment. Frameshift mutations are ruled out due to how HMMER and RPSTBLASTN perform protein calling. Be aware that the `calculate_Ks_TE_hit.py` script can take quite a while, depending on the size of your dataset. For 7 million TE pairs it ran a couple of weeks (using 90 threads simulatenously). 

Software version:
- mafft v7.505
- PAL2NAL v14
- R package seqinr v.4.2-36
- python package numpy v1.26.3
- python package biopython v1.78



### 5. Clustering of candidate HTTs and counting of minimal TE transfer events 
In order to group TE-TE hits that were likely derived from the same transfer event, we used this pipeline to cluster TE-TE hits and count the minimal number of HTT events needed to explain the data. Per clade-clade comparison (either single genomes or collapsed clades involving multiple genomes), we connected TE-TE hits if the percent identity of a within-genome within a clade was higher than between clades. For example, take TE-TE hits A1-B1 and A2-B2 between clade A and clade B, both have a percent identity of 80%. If the percent identity of the blast hit within a single clade (in this case: A1-A2 or B1-B2) has a higher percent identity (>80%), the two TE-TE hits are connected (`first_clustering.py`). 

With TE-TE hits as nodes and TE-TE hit connections as edges per clade-clade comparison, a fast greedy clustering approach is implemented to form hit communities in these graphs (`cluster_fast_greedy.py`). Hit communities having more than 5% of possible edges between them realized were merged (communities with a connectivity greater than 0.05), to prioritize minimizing false positives over false negatives. For computational efficiency, community detection for three clade pairs sharing more than 80,000 candidate HTTs and 3 billion edges was performed using the Leiden algorithm implemented in C++ in the leidenalg package, which relies on the igraph package in Python (`cluster_louvain.py `). Note that connecting TE-TE hits and subsequent clustering will use a lot of memory (>1.5Tb RAM). 

Next, using single linkage clustering, hit communities were connected if they shared at least 1 TE (`SCL_hit_communities.py`). These clusters of hit communities were then projected onto the tree and used to count the minimal number of HTT events needed to explain the data (`count_min_HTTs.py`).  

Note that classifications for TEs are needed in this pipeline. We've constructed these files for de novo TEs per genome, with TE family name in column 1 and TE classifications on order, superfamily and family level in column 7, 8, and 9, respectively. We've used the Repeatmasker classification as well as classifications from PASTEC for this. In case of conflicting classifications, Repeatmasker was favored. For the TEs annotated from the Repeatmasker library (`2759.RepeatMasker.lib`), we've constructed a similar file for the annotations of the TE consensus sequences of the library. 

Software versions: 
- BLAST 2.12.0+
- python package collections v3.13.5
- python package numpy v1.26.3
- python package ete3 v.3.1.3
- python package igraph v0.11.8
- python package leidenalg v0.10.2
- python package networkX v2.5.1

### 6. Retrieving HTT-associated TEs
To be able to infer the proportion of TEs likely associated with horizontal transfer, we clustered all TEs per genome using MMseqs2. We labelled clusters containing a TE involved in transfer as "associated with horizontal transfer". Next, we extract all the TEs per genome that are associated with horizontal transfer. Nested or overlapping TEs are counted only once.

Software versions: 
- BEDtools v2.31.1
- MMseqs2 15.6f452
