#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: process_multi_results.R
##
## Purpose of script: Cleans up multi-alignment results
##
## Author: Karen Wong
##
## Date Created: 06/09/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------
##
## Notes: This script picks the best sequence in each component to represent that sequence
##
## ---------------------------------------------------------------------------------
##

# remove all previous variables
rm(list=ls())

# load up packages  
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(plyr))

# command line options
## ---------------------------------------------------------------------------------

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory")

)
opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
dir = opt$dir

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

# read multi-alignment results and sort by component number
res = fread(paste0(dir, "/multi_results.csv"), stringsAsFactors = F, drop = 3)
colnames(res) = c("component", "group")
res = res[order(res$component),]

# extract cluster count
cluster_count = strsplit(res$group, ",")
res$cluster = sapply(cluster_count,length)

# discard res entries if 0 cluster
res = res[res$cluster!=0,]

# extract item count
# item_count = strsplit(res$group, ',|;')
# res$item = sapply(item_count, length)

# read edge result 
#edge = fread(paste0(dir, "/discovery/assemblytics_component_edge.txt"), stringsAsFactors = F)

# read the input to multi-alginments
assemblytics_path = list.files(path=paste0(dir,"/discovery"), pattern = "^assemblytics_combined_results_with_component_group_*")
assemblytics = NULL
for (i in assemblytics_path){
    print(i)
    assemblytics_per_chr = fread(paste0(dir, "/discovery/", i), stringsAsFactors = F, select = c("ref_chr", "ref_start", "ref_end", "ngap_boundaries_size_left", "ngap_boundaries_size_right", "sample", "INS_id", "component"))
    assemblytics = rbind.data.frame(assemblytics, assemblytics_per_chr)
}

# read sample metadata
metadata = read.table(paste0(dir, '/TMP_sample_metadata.txt'), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source", "project")

# make sure all components in res are also in assemblytics, and vice versa
assemblytics = checkResAssemblyticsComponentConcord(assemblytics, res)
stopifnot(length(which(unique(assemblytics$component) == res$component)) == length(unique(assemblytics$component)))

# pick representative seq per component
assm_tbl = c(table(assemblytics$component))

fileline_start = c(1, cumsum(assm_tbl[1:(length(assm_tbl)-1)])+1)
fileline_end = c(cumsum(assm_tbl[1:(length(assm_tbl))]))

representative_df = NULL
for (i in 1:length(assm_tbl)){
#for (i in 1:2){
  
  print(i)
  start_line = as.numeric(fileline_start[i])
  end_line = as.numeric(fileline_end[i])

  # subset assemblytics dataframe
  assm_sub_df = assemblytics[start_line:end_line,]
  
  representative_df[[i]] = findRepresentativeSeq(assm_sub_df, res)
  
  #representative_df = rbind.data.frame(representative_df, rep_seq, stringsAsFactors = F)
  
}

representative_df = ldply(representative_df, data.frame)

# representative_df = as.data.frame(t(as.data.frame(representative_df, stringsAsFactors = F)), stringsAsFactors = F)
# colnames(representative_df) = c(colnames(assemblytics), "sample_count", "sample_perct", "sample_record", "sample_AFR", "sample_AMR", 'sample_EAS', "sample_EUR", "sample_SAS", "sample_nonAFR", "percent_AFR", "percent_AMR", "percent_EAS", "percent_EUR", "percent_SAS", "percent_nonAFR", "edge_start", "edge_end")
# row.names(representative_df) = 1:nrow(representative_df)

# write output
write.table(representative_df[representative_df$type==0,], paste0(dir, "/discovery/assemblytics_representative_seq.txt"), row.names = F, col.names = T, quote = F, sep = '\t')
write.table(representative_df[representative_df$type==1,], paste0(dir, "/discovery/assemblytics_representative_seq_supplementary.txt"), row.names = F, col.names = T, quote = F, sep = '\t')




















