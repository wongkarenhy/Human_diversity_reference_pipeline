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
suppressMessages(library(GenomicRanges))

# # command line options
## ---------------------------------------------------------------------------------

option_list = list(

  # make_option(c("-d", "--dir"), action="store", default=NA, type='character',
  #             help="work directory"),
  # make_option(c("-c", "--chr"), action="store", default=NA, type='character',
  #             help="list of chromosomes")
  make_option(c("-c", "--chr_list"), action="store", default=NA, type='character',
                help="chr to run")
  

)
opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
#dir = opt$dir
#chr_list = opt$chr
dir = "../"
#chr_list = c(1:22, "X", "Y")
#threads = opt$threads
chr_list = opt$chr_list

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

# read multi-alignment results and sort by component number
# res = NULL
# for (i in 1:length(chr_list)){
#   
#   chr = chr_list[i]
#   # Read individual processed assemblytics file
#   res[[i]] = fread(paste0(dir, "/discovery/multi_results_chr",chr,".csv"), stringsAsFactors = F, drop = 3)
# 
# }
# 
# res = ldply(res, data.frame)
for (chr in chr_list){
  
    print(chr)
  
    res = NULL
    for (size in c("small", "big")){
        res_df = fread(paste0("../discovery_08132019/ma/multi_results_",size,"_chr",chr,".txt"), stringsAsFactors = F, drop = 3)
        res = rbind.data.frame(res, res_df, stringsAsFactors = F)
    }

    colnames(res) = c("component", "group")
    res = res[order(res$component),]
    
    # extract cluster count
    cluster_count = strsplit(res$group, ",")
    res$cluster = sapply(cluster_count,length)
    
    # discard res entries if 0 
    res = res[res$cluster!=0,]
    
    # extract item count
    item_count = strsplit(res$group, ',|;')
    res$item = sapply(item_count, length)
    
    # get size of first cluster
    first_cluster_size = strsplit(str_split_fixed(res$group, ",",2)[,1], ';')
    res$first_cluster_size = sapply(first_cluster_size, length)
    
    # get first cluster percent
    res$first_cluster_perct = res$first_cluster_size/res$item
    
    # read the input to multi-alginments
    assemblytics_path = list.files(path=paste0(dir,"/discovery_08132019/"), pattern = "^assemblytics_combined_results_with_component_group_*")
    keep = grep(paste0("chr",chr,".txt"), assemblytics_path)
    assemblytics_path = assemblytics_path[keep]
    assemblytics = NULL
    for (i in assemblytics_path){
        print(i)
        #assemblytics_per_chr = fread(paste0(dir, "/discovery/", i), stringsAsFactors = F, select = c("ref_chr", "ref_start", "ref_end", "ngap_boundaries_size_left", "ngap_boundaries_size_right", "sample", "INS_id", "component"))
        assemblytics_per_chr = fread(paste0(dir, "/discovery_08132019/", i), stringsAsFactors = F)
        assemblytics = rbind.data.frame(assemblytics, assemblytics_per_chr)
    }
    
    # read sample metadata
    metadata = read.table(paste0(dir, '/TMP_sample_metadata.txt'), stringsAsFactors = F)
    colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source", "project")
    
      disc = which(! assemblytics$component %in% res$component)

      if (length(disc)!=0){
        assemblytics_missing = assemblytics[disc,]
        assemblytics = assemblytics[-disc,]
      }
      
      # sort assemblytics by component
      assemblytics = assemblytics[order(assemblytics$component),]
      
    # assm_tbl = c(table(assemblytics$component))
    # disc = names(assm_tbl)[which(!names(assm_tbl) %in% res$component)]
    # # make sure all components in res are also in assemblytics, and vice versa
    # assemblytics = checkResAssemblyticsComponentConcord(assemblytics, res)
    
    # remove duplicate in res (only if there's error in running)
    disc = which(duplicated(res$component))
    if (length(disc) !=0){
       res = res[-disc,]
    }
    
    stopifnot(length(which(unique(assemblytics$component) == res$component)) == length(unique(assemblytics$component)))
    
    # pick representative seq per component
    assm_tbl = c(table(assemblytics$component))
    
    # make sure the number of items in each component is identical before and after multiple alignments
    
    #assemblytics = assemblytics[!(assemblytics$component %in% disc),]
    
    #assm_tbl = c(table(assemblytics$component))
    
    # now check if the number of items in results is identical the line count per component
    disc = which(assm_tbl != res$item)
    
    if (length(disc)!=0){
        assemblytics_wrong_count = assemblytics[assemblytics$component %in% names(disc),]
        assemblytics = assemblytics[! assemblytics$component %in% names(disc),]
        res = res[!(res$component %in% names(disc)),]
    }
    
    assm_tbl = c(table(assemblytics$component))
    
    
    fileline_start = c(1, cumsum(assm_tbl[1:(length(assm_tbl)-1)])+1)
    fileline_end = c(cumsum(assm_tbl[1:(length(assm_tbl))]))
    
    representative_df = NULL
    singleton_clus_df = NULL
    for (i in 1:length(assm_tbl)){
    #for (i in 1:37){
      
      print(i)
      start_line = as.numeric(fileline_start[i])
      end_line = as.numeric(fileline_end[i])
    
      # subset assemblytics dataframe
      assm_sub_df = assemblytics[start_line:end_line,]
      
      singleton_clus_df[[i]] = findSingletonClus(assm_sub_df, res, i, metadata)
      representative_df[[i]] = findRepresentativeSeq(assm_sub_df, res, i, metadata)
      
      #representative_df = rbind.data.frame(representative_df, rep_seq, stringsAsFactors = F)
      
    }
    
    representative_df = ldply(representative_df, data.frame)
    singleton_clus_df = ldply(singleton_clus_df, data.frame)

    # write output
    write.table(representative_df, paste0(dir, "/discovery/assemblytics_representative_seq_raw_chr",chr,".txt"), row.names = F, col.names = T, quote = F, sep = '\t')
    
}

# res = mclapply(chr_list, processChromosomes, mc.cores = threads)
# job_err = res[1:length(chr_list)]
# stopifnot(length(which(grepl("Error", job_err)))==0)
# 




df2 = fread("../discovery/assemblytics_representative_seq_raw_chr10.txt", stringsAsFactors = F)











