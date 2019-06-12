#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: calc_edge.R
##
## Purpose of script: calculate edge score for each component
##
## Author: Karen Wong
##
## Date Created: 06/09/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------
##
## Notes: This script calculate the breakpoint consistency for all insertions in each component
##
## ---------------------------------------------------------------------------------
##

# remove all previous variables
rm(list=ls())

# load up packages  
suppressMessages(library(optparse))
suppressMessages(library(data.table))

# command line options
## ---------------------------------------------------------------------------------

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory")

)

opt = parse_args(OptionParser(option_list=option_list))

## ---------------------------------------------------------------------------------

# parse user input
dir = opt$dir

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))


# read input 
assemblytics = fread(paste0(dir, "/discovery/assemblytics_combined_results_with_component_group.txt"), stringsAsFactors = F, header = T)

# assemblytics = assemblytics[order(assemblytics$component),]
assm_tbl = c(table(assemblytics$component))

start_counter = 1
edge_start = NULL
edge_end = NULL
for (i in 1:length(assm_tbl)){

    end = start_counter + assm_tbl[i] - 1
    assm_sub_df = assemblytics[start_counter:end,]
    
    assm_sub_df = assm_sub_df[assm_sub_df$ngap_boundaries!="yes",]
    
    # Duplicate the lines if haplo=="unphased"
    dup = assm_sub_df[assm_sub_df$haplo=="unphased", ]
    assm_sub_df = rbind.data.frame(assm_sub_df, dup, stringsAsFactors = F)
    
    # calculate edge score for breakpoint start and breakpoint end
    edge_start[[i]] = cal_edge_size(assm_sub_df$ref_start)
    edge_end[[i]] = cal_edge_size(assm_sub_df$ref_end)
    
    start_counter = end + 1
}

assemblytics_component_edge = data.frame(names(assm_tbl), edge_start, edge_end, stringsAsFactors = F)
colnames(assemblytics_component_edge) = c("component", "edge_start", "edge_end")

write.table(assemblytics_component_edge, paste0(dir, "/discovery/assemblytics_component_edge.txt"), row.names = F, col.names = T, quote = F)




