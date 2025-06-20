# Horizontal transfer TE fungi
Scripts used to detect horizontal transfer of transposable elements (HTT) in 1,348 publicly available fungal genomes (JGI), see [reference]. Folders include Snakemake pipelines, scripts, and a file containing software versions used for each specific pipeline. All scripts included in these pipelines are created by Josje Romeijn (MGE lab, Utrecht University). 

1. Extraction of TEs from fungal genomes 
2. TE classification 
3. BUSCO Ks estimations of genome comparisons
4. TE vs TE blast
5. Calculation of Ks of TE pairs to obtain candidate HTTs 
6. Clustering of candidate HTTs and counting of minimal TE transfer events 
7. Retrieving HTT-associated TEs

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
This pipeline starts by creating a file that identifies the best BUSCO gene comparisons for each comparison of 2 nodes in a phylogenetic tree. How it works is explained in the tree below. 

First, for node 20 is determined which fungal BUSCO set is most suitable (priority: class-level, order-level, phylum-level, kingdom-level BUSCO set). Next, for each BUSCO gene in the determined BUSCO set the two longest copies of this BUSCO gene are searched in the two children clades. In this case, genome A and D both have the longest copy, and are added for this BUSCO gene to the list. 

Format of `list_busco_seqs_species.txt`:
{node_name}__{species1}__{species2}__{busco_set}__{busco_gene}
e.g.:
node1__Zymtr1__Phaal1__fungi__100957at4751
