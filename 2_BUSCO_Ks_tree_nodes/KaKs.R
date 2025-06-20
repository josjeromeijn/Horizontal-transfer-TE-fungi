#calculate Ks using seqinr package
library("seqinr")

#get input and output from command
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
species_a <- args[3]
species_b <- args[4]
busco <- args[5]

aln <- read.alignment(file = input_file, format = "fasta")
data <- kaks(aln, verbose = FALSE)



cat(file = output_file, c(data$ka, data$ks, species_a, species_b, busco), sep = "\t")
