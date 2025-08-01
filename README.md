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

### 3. TE vs TE blast
Make sure to only store .fasta files of fungal genomes included in this analysis in `genome_folder`, as the pipeline scans this directory for .fasta files for input. Also, it is important to note for further analyses that we removed completely identical TE sequences per genome to reduce the workload. For later analyses, these "duplicate" TEs will have to be manually added. 

From the BUSCO comparisons, we estimated which nodes in the tree contain species that are too recently diverged to detect horizontal transfer (referred to as "collapsed nodes"). We've compiled a list of all species comparisons that are for that reason excluded from further analysis (`collapsed_sp_sp_combinations.txt`). TE blast hits involving TEs of a collapsed node are therefore filtered out, as well as "self"-hits (TE hits between TEs residing in the same genome). TE blast hits are then chained if they involve the same two TEs, their gap is <600 bp and overlap <600 bp. Next, only TE hits are retained where the alignment covers >60% of the length of both TEs. 

Software versions: 
- BLAST 2.12.0+
- python package collections v3.13.5
- python package numpy v1.26.3

### 4. Calculation of Ks of TE pairs to obtain candidate HTTs 

