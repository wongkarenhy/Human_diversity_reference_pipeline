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
cluster_count = (strsplit(res$group, ";"))
res$cluster = sapply(cluster_count,length)

# read edge result 
edge = fread(paste0(dir, "/discovery/assemblytics_component_edge.txt"), stringsAsFactors = F)

# read the input to multi-alginments
assemblytics_path = list.files(path=paste0(dir,"/discovery"), pattern = "^assemblytics_combined_results_with_component_group_*")
assemblytics = NULL
for (i in assemblytics_path){
    
    assemblytics_per_chr = fread(paste0(dir, "/discovery/", i), stringsAsFactors = F, select = "component")
    assemblytics = rbind.data.frame(assemblytics, assemblytics_per_chr)
}

# read sample metadata
metadata = read.table(paste0(dir, '/TMP_sample_metadata.txt'), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source", "project")

# make sure all components in res are also in assemblytics
assemblytics = checkResAssemblyticsComponentConcord(assemblytics, res)

# pick representative seq per component
assm_tbl = c(table(assemblytics$component))

start_counter = 1
representative_df = NULL
for (i in 1:length(assm_tbl)){
#for (i in 1:7006){
  end = start_counter + assm_tbl[i] - 1

  representative_df = findRepresentativeSeq(start_counter, end)
  
  start_counter = end + 1
  
}

representative_df = as.data.frame(t(as.data.frame(representative_df, stringsAsFactors = F)), stringsAsFactors = F)
colnames(representative_df) = c(colnames(assemblytics), "sample_count", "sample_perct", "sample_record", "sample_AFR", "sample_AMR", 'sample_EAS', "sample_EUR", "sample_SAS", "sample_nonAFR", "percent_AFR", "percent_AMR", "percent_EAS", "percent_EUR", "percent_SAS", "percent_nonAFR", "edge_start", "edge_end")
row.names(representative_df) = 1:nrow(representative_df)

# write output
write.table(representative_df, paste0(dir, "/discovery/assemblytics_representative_seq.txt"), row.names = F, col.names = T, quote = F, sep = '\t')




















