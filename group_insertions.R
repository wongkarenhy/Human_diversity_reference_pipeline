#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: group_insertions.R
##
## Purpose of script: Group insertions if they have overlapping reference coordinates. 
##
## Author: Karen Wong
##
## Date Created: 06/07/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------
##
## Notes: Insertions with overlapping reference coordinates are grouped together into components. 
##        Components with only a single sample (singleton) are removed to speed up steps downstream.
## ---------------------------------------------------------------------------------
##

# remove all previous variables
rm(list=ls())

# load up packages  
suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(igraph))
suppressMessages(library(data.table))
#suppressMessages(library(plyr))

# command line options
## ---------------------------------------------------------------------------------

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory"),
  make_option(c("-t", "--threads"), action="store", default=NA, type='numeric',
              help="number of threads")
  
)
opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
dir = opt$dir
threads = opt$threads

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

# read the combined assemblytics text file
assemblytics_combined = fread(paste0(dir, "/discovery/assemblytics_combined_results.txt"), stringsAsFactors = F, header = T)

# read metadata
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source", "project")

# add 0.5 to each end before creating grange objects
assemblytics_combined$ref_start = assemblytics_combined$ref_start-0.5
assemblytics_combined$ref_end = assemblytics_combined$ref_end+0.5

# igraph can't handle long vector
# break assemblytics_combined dataframe by chr 
#for (i in unique(assemblytics_combined$ref_chr)){
assemblytics_combined_added = NULL
assemblytics_combined_added = mclapply(unique(assemblytics_combined$ref_chr), assignComponent, mc.cores = threads)
comp_per_chr = sapply(assemblytics_combined_added, function(x) length(unique(x$component)))
comp_per_chr_cumsum = c(0,cumsum(comp_per_chr)[-length(cumsum(comp_per_chr))])

for (i in 1:length(assemblytics_combined_added)){
  
    # adjust the component number
    assemblytics_combined_added[[i]]$component = assemblytics_combined_added[[i]]$component+comp_per_chr_cumsum[i]
    assemblytics_combined_added_df_per_chr = as.data.frame(assemblytics_combined_added[[i]], stringsAsFactors = F)
    
    # sort by component number before writing an output file
    assemblytics_combined_added_df_per_chr = assemblytics_combined_added_df_per_chr[order(assemblytics_combined_added_df_per_chr$component),]
    
    # correct ref_start and ref_end values
    assemblytics_combined_added_df_per_chr$ref_start = assemblytics_combined_added_df_per_chr$ref_start+0.5
    assemblytics_combined_added_df_per_chr$ref_end = assemblytics_combined_added_df_per_chr$ref_end-0.5
    
    # partition dataframe by size and chr
    assemblytics_combined_added_df_per_chr_small = assemblytics_combined_added_df_per_chr[assemblytics_combined_added_df_per_chr$insert_size<50,]
    assemblytics_combined_added_df_per_chr_big = assemblytics_combined_added_df_per_chr[assemblytics_combined_added_df_per_chr$insert_size>=50,]
    
    # define chr
    chr = assemblytics_combined_added_df_per_chr_small$ref_chr[1]
    
    # write output
    options(scipen=999)
    write.table(assemblytics_combined_added_df_per_chr_small, paste0(dir, "/discovery/assemblytics_combined_results_with_component_group_small_", chr,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')
    write.table(assemblytics_combined_added_df_per_chr_big, paste0(dir, "/discovery/assemblytics_combined_results_with_component_group_big_", chr,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')
    
}

# make list into individual dataframes 



# convert a list of dataframes into a single dataframe
#assemblytics_combined_added <- ldply(assemblytics_combined_added, data.frame)

# # find overlapping reference ranges
# assemblytics.gr = makeGRangesFromDataFrame(assemblytics_combined, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
# assemblytics_assemblytics_ov_self = findOverlaps(assemblytics.gr, assemblytics.gr, type = "any", ignore.strand=T)
# 
# # group insertions based on overlapping reference ranges
# # long vector error if edges are too large 
# # in this case we need to partition the data into different chrs before grouping them into components
# edges = as.data.frame(assemblytics_assemblytics_ov_self)
# 
# 
# net = graph_from_data_frame(d=edges, directed = F)
# membership = components(net)$membership
# components_count = count_components(net)
# assemblytics_combined_added$component = membership # add the component number to each INS

# # remove singleton (components with just one sample)
# singleton = names(which(table(assemblytics_combined_added[,c("component")])==1))
# doubleton = names(which(table(assemblytics_combined_added[,c("component")])==2))
# doubleton_df = assemblytics_combined_added[assemblytics_combined_added$component %in% doubleton,]
# doubleton_df = doubleton_df[order(doubleton_df$component),]
# doubleton_df_odd = doubleton_df[seq(1,nrow(doubleton_df),2),]
# doubleton_df_even = doubleton_df[seq(2,nrow(doubleton_df),2),]
# doubleton_df_odd$even = doubleton_df_even$sample
# # here doubleton is defined as two INS but from the same sample
# doubleton = doubleton_df_odd$component[doubleton_df_odd$sample==doubleton_df_odd$even]
# assemblytics_combined_added = assemblytics_combined_added[!(assemblytics_combined_added$component %in% c(singleton, doubleton)), ]

# # sort by component number before writing an output file
# assemblytics_combined_added = assemblytics_combined_added[order(assemblytics_combined_added$component),]
# 
# # correct ref_start and ref_end values
# assemblytics_combined_added$ref_start = assemblytics_combined_added$ref_start+0.5
# assemblytics_combined_added$ref_end = assemblytics_combined_added$ref_end-0.5
# 
# # partition dataframe by size and chr
# options(scipen=999)
# assemblytics_combined_added_small = assemblytics_combined_added[assemblytics_combined_added$insert_size<50,]
# assemblytics_combined_added_big = assemblytics_combined_added[assemblytics_combined_added$insert_size>=50,]
# chr_list = paste0("chr", c(1:22,"X","Y"))
# for (chr in chr_list){
#   
#   assemblytics_combined_added_small_chr = assemblytics_combined_added_small[assemblytics_combined_added_small$ref_chr==chr,]
#   write.table(assemblytics_combined_added_small_chr, paste0(dir, "/discovery/assemblytics_combined_results_with_component_group_small_", chr,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')
# 
#   assemblytics_combined_added_big_chr = assemblytics_combined_added_big[assemblytics_combined_added_big$ref_chr==chr,]
#   write.table(assemblytics_combined_added_big_chr, paste0(dir, "/discovery/assemblytics_combined_results_with_component_group_big_", chr,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')
#   
# }
# 
# 
# 
# 
