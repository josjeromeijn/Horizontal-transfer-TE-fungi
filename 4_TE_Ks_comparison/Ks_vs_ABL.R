#Made by: Josje Romeijn 13-05-24
#script to visualize Ks values of BUSCO genes 1411 genome-genome comparisons
library(stringr)
library(ggplot2)
library(patchwork)
library(treedataverse)
library(ggtree)

#----------create matrix of RAD values of comparisons---------------------------
args <- commandArgs(trailingOnly = TRUE)
working_dir <- args[1]
output_dir <- args[2]



working_dir <- "/your/path/"
output_dir <- "/your/path/"
setwd(working_dir)

tree <- read.tree("/your/path/to/tree") #rooted, RED-scaled tree
original_tree <- read.tree("/your/path/to/tree") #rooted, unscaled tree (only used for visualization later on)

distMat <- cophenetic(tree)

#----------------------load in nodes and Ks values------------------------------
nodes <- list.files(paste0(output_dir,"/Ks/"))

#remove Ks temp directory
nodes <- nodes[-which(nodes == "Ks_temp")]

#create dataframe to store results
df <- data.frame(nodes = str_remove(nodes, "_KaKs.txt"))

#--------------get Ks value of all busco seqs per node and filter---------------
df[, c("Ks", "Ks_nr", "Ks_nr_filtering", "quantile")] <- t(sapply(1:nrow(df), function(i) {
  data <- read.table(paste0(output_dir,"/Ks/",nodes[i]), sep = "\t", fill = TRUE)
  c(Ks = mean(data$V2[data$V2 < 9]),
    Ks_nr = length(data$V2),
    Ks_nr_filtering = length(data$V2[data$V2 < 9]),
    quantile = quantile(data$V2[data$V2 < 9], 0.005, na.rm = TRUE) )
}))

#-------------------add column with the used BUSCO set--------------------------
#first read in file with list of all nodes and their busco sets
system("awk -F '__' '!seen[$1]++' list_busco_seqs_species.txt > list_busco_seqs_species_filtered.txt")
list_nodes_busco <- read.table("list_busco_seqs_species_filtered.txt", sep = "\t")

#for each node, find a line starting with this node in the file and filter out the busco set
df$busco_set <- sapply(1:nrow(df), function(i) {
  id <- which(grepl(paste0("^", df$nodes[i], "__"), unlist(list_nodes_busco)))
  str_split(list_nodes_busco[id,], "__")[[1]][4]
})

#------------------------get average branch length------------------------------
df$RAD <- sapply(1:nrow(df), function(i) {
  #get node id for this node name
  id = nodeid(tree,df$nodes[i])
  
  #get two children clades
  children <- child(tree, id)
  
  #get random child from each clade (and if clade is tip, keep the tip)
  results <- lapply(children, function(child) {
    first_offspring <- offspring(tree, child, type = "tips")[1]
    if (length(first_offspring) == 0 || is.na(first_offspring)) {
      return(child)
    }
    return(first_offspring)
  })
  
  #find distance between these 2 random children (in dist matrix)
  #(tree is normalized, all children should give same phylogenetic distance)
  distMat[which(colnames(distMat) == tree$tip.label[results[[1]]]), which(colnames(distMat) == tree$tip.label[results[[2]]])] / 2
})





#-----------------add column with "original number of busco seqs"---------------
#use awk because it's quicker
#(to check that everything worked in snakemake)
orig_num_busco_seq <- system("awk -F '__' '{count[$1]++} END {for (num in count) print num, count[num]}' list_busco_seqs_species.txt", intern = TRUE)
orig_num_busco_seq <- as.data.frame(do.call(rbind, strsplit(orig_num_busco_seq, " ")), stringsAsFactors = FALSE)

df$orig_num_bsc <- sapply(1:nrow(df), function(i) {
  orig_num_busco_seq$V2[which(orig_num_busco_seq$V1 == df$nodes[i])]
})






#-----------remove nodes that only have oversaturated Ks values-----------------
#(Ks > 9, so no reliable Ks estimation)
print(paste("number of comparisons with saturation of all BUSCO seqs:", length(which(is.na(df$Ks)))))
#df <- df[-c(which(is.na(df$Ks))),]

#give this node the same Ks value as the highest Ks of the two children nodes
children <- child(tree, nodeid(tree, df$nodes[which(is.na(df$Ks))]))
df$quantile[which(is.na(df$Ks))] <- min(df$quantile[which(df$nodes == nodelab(tree, children[1]))], df$quantile[which(df$nodes == nodelab(tree, children[2]))])

#---add column with busco level and default number of busco seqs per busco set--
busco_sets <- data.frame(sets= c(
  "fungi","ascomycota","basidiomycota","microsporidia","mucoromycota",
  "agaricomycetes", "boletales", "capnodiales","chaetothyriales",
  "agaricales", "dothideomycetes", "eurotiomycetes", "leotiomycetes",
  "saccharomycetes", "sordariomycetes", "tremellomycetes", "polyporales",
  "eurotiales", "glomerellales", "helotiales", "hypocreales", "mucorales",
  "onygenales", "pleosporales"), 
  nr_busco_seqs= c(
    758, 1706, 1764, 600, 1614, 2898, 4878, 3578, 6265, 3870, 3786, 3546, 3234, 
    2137, 3817, 4284, 4464, 4191, 6841, 5177, 4494, 2449, 4862, 6641),
  level=c(
    "kingdom", "phylum", "phylum", "phylum", "phylum", "class", "order", "order",
    "order", "order", "class", "class", "class", "class", "class", "class",
    "order","order","order","order","order","order","order","order"))

df[,c("num_busco_of_busco_set", "busco_level")] <- t(sapply(1:nrow(df), function(i) {
  set_index <- match(df$busco_set[i], busco_sets$sets)
  c(num_busco_of_busco_set = as.numeric(busco_sets$nr_busco_seqs[set_index]),
    busco_level = busco_sets$level[set_index])
}))

#-------------add column with percentages of kept busco seqs--------------------
df$percentage_kept <- df$Ks_nr_filtering / df$Ks_nr
df$percantage_kept_of_busco_set <- df$Ks_nr_filtering / as.numeric(df$num_busco_of_busco_set)


#--------------------------------plots------------------------------------------
#plot Ks ~ RAD with circle sizes indicating how many busco seqs were kept after 
#filtering for oversaturation
ggplot(data = df, aes(x = RAD, y = Ks, color = "black", fill = busco_level, size = percentage_kept)) + 
  geom_point(shape = 21) + scale_size_continuous(name = "Percentage kept", labels = scales::percent_format(), guide = "legend") + 
  scale_color_identity() + scale_fill_discrete() + theme_bw()

#plot Ks ~ RAD with circle sizes indicating how many busco seqs were present in
#node AND were kept after filtering for oversaturation
ggplot(data = df, aes(x = RAD, y = Ks, color = "black", fill = busco_level, size = percantage_kept_of_busco_set)) + 
  geom_point(shape = 21) + scale_size_continuous(name = "Percentage kept", labels = scales::percent_format(), guide = "legend") + 
  scale_color_identity() + scale_fill_discrete() + theme_bw()

#plot Ks ~ RAD, colored whether quantile value is > or < 0.01 (cut-off)
ggplot(data = df, aes(x = RAD, y = Ks, color = "black", fill = quantile > 0.01)) + 
   geom_point(shape = 21, size = 2) + 
   scale_color_identity() + scale_fill_discrete() + theme_bw()

#check to make sure all busco seqs did run in snakemake (you should see diagonal)
plot(df$Ks_nr,df$orig_num_bsc)


#reorder df based on RAD values
df <- df[order(df$RAD),]

#make a 1 plot containing the 96 lowest RAD's and their individual Ks distribution
plots <- list()
for (j in 1:1) {
  for (i in 1:nrow(df))  {
    #find index in dataframe
    idx <- i+24*(j-1)
    
    #load in data
    data <- as.data.frame(read.table(paste0(output_dir,"/Ks/",df$nodes[idx],"_KaKs.txt"), sep = "\t", fill = TRUE)) 
    #only keep Ks values below 9
    data <- data[data$V2 < 9,]
    
    #create plot
    p <- ggplot(data, aes(x=V2, ..scaled..)) + geom_density(alpha=.5, fill="#FF6666") + 
      theme_classic() + theme(plot.title = element_text(size=11)) + 
      labs(title = paste("Node number = ", df$nodes[idx], "\nRAD =", round(df$RAD[idx], digits = 4), 
                         "| busco set =", df$busco_set[idx],
                         "\n[#seqs in set|#seqs present|#seqs Ks<9] : [", df$num_busco_of_busco_set[idx],
                         "|", df$Ks_nr[idx], "|", df$Ks_nr_filtering[idx],"]") , x= "Ks") +
      geom_vline(xintercept = df$quantile[idx], color = "darkred", linewidth = 1) +
      annotate("text", y=0.9, x =df$quantile[idx] + 0.8, label = round(df$quantile[idx], digits = 3)) + lims(x=c(0,9), y=c(0,1)) +
      coord_cartesian(expand = FALSE)
    plots[[i]] <- p #add to list
  }
  ggsave(wrap_plots(plots,ncol =4), filename = paste("Ks_distr_lowest_RADs_", (j-1)*nrow(df),"_", j*nrow(df),".pdf", sep = ""), device = "pdf",
         height = 1200, width = 20, units = "in", limitsize = FALSE)
}


#---------------------find missing values---------------------------------------
#find missing values
for (i in c(1:nrow(df))) {
  if (!(i %in% as.numeric(df$nodes))) {
    print(paste("missing node is node", i))
  } 
}




#-----------------visualize collapsed nodes in tree-----------------------------


#get nodes that are going to be collapsed
nodes = sapply(as.numeric(df$nodes[which(df$quantile < 0.01)]), function (i) nodeid(tree, i))

#get phylum data
data <- read.table("earlgray_stats_1348_genomes.csv", sep = "\t", header = TRUE, fill = TRUE)
lineages <- strsplit(data$lineage, ",")
data$phylum <- sapply(1:1348, function(i) lineages[[i]][2])
tree$phylum <- sapply(1:length(tree$tip.label), function(i) data$phylum[which(data$genome.name == str_remove(tree$tip.label[i], "_fungi"))])
tree$phylum <- as.factor(tree$phylum)

#create dataframe containing MRCA nodes of phyla
y <- as_tibble(tree)
indiv_per_phyla <- lapply(unique(tree$phylum), function(i) which(tree$phylum == i))
indiv_per_phyla <- lapply(indiv_per_phyla, function(i) toString(i))
#second, find MRCA of tip numbers per phyla
phyla_node_number <- sapply(1:length(indiv_per_phyla), function(i) eval(parse(text = paste("MRCA(y, c(", indiv_per_phyla[i], "), root.edge=FALSE)"))))
#put it in dataframe with corresponding phylum name
d <- data.frame(node=c(unlist(phyla_node_number[2,])), Phyla=unique(tree$phylum))
d$color <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
             "#44AA99", "#999933", "#882255")



#plot tree
gg_norm <- ggtree(tree, size=0.2) + geom_point2(aes(subset = (node %in% nodes)), color= 'red', size = 1, alpha=1/4) +
  geom_cladelab(node = d$node, label = d$Phyla, barcolour = d$color, barsize = 3, angle = 45, fontsize = 3)
gg_orig <- ggtree(original_tree, size=0.2) + geom_point2(aes(subset = (node %in% nodes)), color= 'red', size = 1, alpha=1/4)

#visualize both trees and save
ggsave(gg_norm + gg_orig, filename = "collapsed_nodes_tree.pdf", device = "pdf",
       height = 20, width = 30, units = "in")


#------create files with collapsed nodes species-species comparisons------------

#how many species-species comparisons are lost? 
#find nodes to collapse
collapse_nodes <- as.numeric(df$nodes[which(df$quantile < 0.01)])

#now turn into node_ids 
nodes_to_collapse = data.frame(node_names = collapse_nodes, node_ids = c(sapply(collapse_nodes, function(i) nodeid(tree, i))))

internal_collapsed_nodes <- list()
for (node in nodes_to_collapse$node_ids) {
  parent_id <- parent(tree, node)
  if (parent_id %in% nodes_to_collapse$node_ids) {
    print("parent node present!")
    print(node)
    internal_collapsed_nodes <- append(internal_collapsed_nodes, node)
  }
}

#remove these "internal" nodes to collapse
nodes_to_collapse <- nodes_to_collapse[-c(which(nodes_to_collapse$node_ids %in% internal_collapsed_nodes)),]

#create a file that per line puts all the genomes in a collapsed node.
nodes_to_collapse$sp <- sapply(nodes_to_collapse$node_ids, function(i) str_remove(na.omit(tree$tip.label[offspring(tree, i)]), "_fungi"))
lines <- sapply(nodes_to_collapse$sp, function(x) paste(x, collapse = ' '))
writeLines(lines, con = "collapsed_nodes_with_species.txt")
write.table(nodes_to_collapse, file="collapsed_nodes_with_species.txt", row.names=FALSE, quote=FALSE, sep = "\t")

#for each node, find all the tips (=offspring), and calculate number of possible
#species-species comparisons in this clade (ncol(combn)), and take sum of this
#to get all species-species comparisons that are collapsed
collapsed_sp_sp_comp <- sum(sapply(nodes_to_collapse$node_ids, function(i) ncol(combn(length(offspring(tree, i)), 2))))

#!!!!!!!!! GOES WRONG BECAUSE IT ALSO COUNTS INTERNAL NODES!!!!!!!!!!


#get list of sp-sp comparisons that will be collapsed
comb<- sapply(nodes_to_collapse$node_ids, function(i) {
  #of all offspring, find the name
  offspring_ids <- str_remove(tree$tip.label[offspring(tree, i)], "_fungi")
  
  #remove NA's -> internal nodes
  offspring_ids <- offspring_ids[!is.na(offspring_ids)]
  
  #get all possible combinations
  combinations <- combn(offspring_ids, 2)
  t(combinations)
})

#make this one big list
sp_sp_comb <- do.call(rbind, comb)

#save lists of collapsed nodes and collapsed sp-sp combinations
write.table(sp_sp_comb, file="collapsed_sp_sp_combinations.txt",row.names=FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
write.table(nodes_to_collapse, file="collapsed_nodes.txt",row.names=FALSE, quote = FALSE, sep = "\t")

#------create table with species-species comparisons of not-collapsed nodes-----


#make a table containing all nodes of tree and their species-species comparisons
all_sp_sp_comb <- data.frame(nodes = df$nodes, quantile = df$quantile)
all_sp_sp_comb$nodeid <- sapply(all_sp_sp_comb$nodes, function(i) nodeid(tree, i))

all_sp_sp_comb$comb <- lapply(all_sp_sp_comb$nodeid, function(i) {
  children_names <- list()
  #look per children clade
  for (j in c(1,2)) {
    #if it is a tip, assign name of this tip
    if (isTip(tree, child(tree, i)[j])) {
      children_names[[j]] <- str_remove(tree$tip.label[child(tree, i)[j]], "_fungi")
    } else { #if it is not a tip, get all names of the tips of this branch
      offspring <- offspring(tree,child(tree,i)[j])
      offspring_tips <- offspring[isTip(tree, offspring)]
      children_names[[j]] <- str_remove(tree$tip.label[offspring_tips], "_fungi")
    }
  }
  unlist(lapply(children_names[[1]], function(k) lapply(children_names[[2]], function(l) paste0(k,"__",l))))
})

all_sp_sp_comb <- unnest(all_sp_sp_comb, comb)

#save this table
write.table(all_sp_sp_comb, file="all_sp_sp_combinations.txt", row.names=FALSE, quote=FALSE, sep = "\t")



