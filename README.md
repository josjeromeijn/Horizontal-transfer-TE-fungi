# Horizontal transfer TE fungi
Scripts used to detect horizontal transfer of transposable elements (HTT) in 1,348 publicly available fungal genomes (JGI), see [reference]. Folders include Snakemake pipelines, scripts, and a file containing software versions used for each specific pipeline. 

1. Extraction of TEs from fungal genomes 
2. TE classification 
3. BUSCO Ks estimations of genome comparisons
4. TE vs TE blast
5. Calculation of Ks of TE pairs to obtain candidate HTTs 
6. Clustering of candidate HTTs and counting of minimal TE transfer events 
7. Retrieving HTT-associated TEs

## General software versions used for all pipelines: 
Snakemake v8.12.0
Python package pandas v2.2.3
Python package tqdm v4.67.1

## 1. Extraction of TEs from fungal genomes
