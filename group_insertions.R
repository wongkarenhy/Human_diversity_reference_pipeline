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

# command line options
## ---------------------------------------------------------------------------------

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory"),
  make_option(c("-v", "--vcf"), action="store", default=NA, type='character',
              help="EEE's vcf file"),
  make_option(c("-t", "--threads"), action="store", default=NA, type='numeric',
              help="number of threads")
  
)
opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
dir = opt$dir

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

# find overlapping reference ranges
assemblytics.gr = makeGRangesFromDataFrame(assemblytics_combined, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
assemblytics_assemblytics_ov_self = findOverlaps(assemblytics.gr, assemblytics.gr, type = "any", ignore.strand=T)

# group insertions based on overlapping reference ranges
edges = as.data.frame(assemblytics_assemblytics_ov_self)
net = graph_from_data_frame(d=edges, directed = F)
membership = components(net)$membership
components_count = count_components(net)
assemblytics_combined$component = membership # add the component number to each INS

# remove singleton (components with just one sample)
singleton = names(which(table(assemblytics_combined[,c("component")])==1))
doubleton = names(which(table(assemblytics_combined[,c("component")])==2))
doubleton_df = assemblytics_combined[assemblytics_combined$component %in% doubleton,]
doubleton_df = doubleton_df[order(doubleton_df$component),]
doubleton_df_odd = doubleton_df[seq(1,nrow(doubleton_df),2),]
doubleton_df_even = doubleton_df[seq(2,nrow(doubleton_df),2),]
doubleton_df_odd$even = doubleton_df_even$sample
# here doubleton is defined as two INS but from the same sample
doubleton = doubleton_df_odd$component[doubleton_df_odd$sample==doubleton_df_odd$even]
assemblytics_combined = assemblytics_combined[!(assemblytics_combined$component %in% c(singleton, doubleton)), ]

# sort by component number before writing an output file
assemblytics_combined = assemblytics_combined[order(assemblytics_combined$component),]

# correct ref_start and ref_end values
assemblytics_combined$ref_start = assemblytics_combined$ref_start+0.5
assemblytics_combined$ref_end = assemblytics_combined$ref_end-0.5

options(scipen=999)
write.table(assemblytics_combined, paste0(dir, "/discovery/assemblytics_combined_results_with_component_group.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

