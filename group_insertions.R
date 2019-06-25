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
    
    # for each component, make sure each entry comes from a different haplotype
    comp_tbl = table(assemblytics_combined_added_df_per_chr$component)
    
    comp_start = c(1,cumsum(comp_tbl[1:(length(comp_tbl)-1)])+1)
    comp_end = c(cumsum(comp_tbl[1:length(comp_tbl)]))
    
    discard = NULL
    for (comp in 1:length(comp_tbl)){
      
      if (comp_tbl[comp]==1) next
      
      comp_df = assemblytics_combined_added_df_per_chr[comp_start[comp]:comp_end[comp],]
      
      if (all(data.frame(table(comp_df[,c("sample", "haplo")]))$Freq < 2)) next
      
      else {
        
        count_df = data.frame(table(comp_df[,c("sample", "haplo")]))
        problematic_vec = which(count_df$Freq>=2)
        
        # identify which one we want to discard
        for (p in problematic_vec){
          
            problematic_df = comp_df[comp_df$sample==count_df$sample[p] & comp_df$haplo==count_df$haplo[p],]
            problematic_df = problematic_df[order(problematic_df$ref_gap_size),]
            
            # don't do anything if the ref coords are not overlapping
            problematic_df.gr = makeGRangesFromDataFrame(problematic_df, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end", ignore.strand = T)
            ov = findOverlaps(problematic_df.gr, problematic_df.gr, type = "any")
            
            if (length(ov[ov@from!=ov@to])==0) next
            
            if (length(unique(problematic_df$insert_size))!=1) { # if the two scaffolds report different insert size, we will discard them both
            
                discard = c(discard, as.vector(problematic_df$INS_id))
              
            } else { # discard the one with shorter scaffold length
                
              keep = problematic_df$INS_id[which.max(problematic_df$scaffold_length)] # we keep the one with the longest scaffold length
              
              discard = c(discard, as.vector(problematic_df$INS_id[!problematic_df$INS_id %in% keep])) # rm anything not the longest scaffold
              
            }
            
        }
        
      }
      
    }
    
    # discard redundant coords
    assemblytics_combined_added_df_per_chr = assemblytics_combined_added_df_per_chr[!assemblytics_combined_added_df_per_chr$INS_id %in% discard, ]
    
    # write singleton to a separate file (components with just one sample)
    singleton = names(which(table(assemblytics_combined_added_df_per_chr[,c("component")])==1))
    doubleton = names(which(table(assemblytics_combined_added_df_per_chr[,c("component")])==2))
    doubleton_df = assemblytics_combined_added_df_per_chr[assemblytics_combined_added_df_per_chr$component %in% doubleton,]
    doubleton_df = doubleton_df[order(doubleton_df$component),]
    doubleton_df_odd = doubleton_df[seq(1,nrow(doubleton_df),2),]
    doubleton_df_even = doubleton_df[seq(2,nrow(doubleton_df),2),]
    doubleton_df_odd$even = doubleton_df_even$sample
    # here doubleton is defined as two INS but from the same sample
    doubleton = doubleton_df_odd$component[doubleton_df_odd$sample==doubleton_df_odd$even]
    singleton_per_chr_df = assemblytics_combined_added_df_per_chr[(assemblytics_combined_added_df_per_chr$component %in% c(singleton, doubleton)), ]
    # When output singleton, make sure every component is output exactly once (not twice as seen in both haplotypes)
    singleton_per_chr_df = singleton_per_chr_df[-which(duplicated(singleton_per_chr_df$component)),]
    
    # Remove singleton entries from the original dataframe
    assemblytics_combined_added_df_per_chr = assemblytics_combined_added_df_per_chr[!(assemblytics_combined_added_df_per_chr$component %in% c(singleton, doubleton)), ]
    
    # get component number if any insert_size is greater than 50
    large_comp = unique(assemblytics_combined_added_df_per_chr$component[assemblytics_combined_added_df_per_chr$insert_size>=50])
    
    # partition dataframe by size and chr
    assemblytics_combined_added_df_per_chr_small = assemblytics_combined_added_df_per_chr[!(assemblytics_combined_added_df_per_chr$component %in% large_comp), ]
    assemblytics_combined_added_df_per_chr_big = assemblytics_combined_added_df_per_chr[assemblytics_combined_added_df_per_chr$component %in% large_comp, ]
    
    # define chr
    chr = assemblytics_combined_added_df_per_chr_small$ref_chr[1]
    
    # write output
    options(scipen=999)
    write.table(assemblytics_combined_added_df_per_chr_small, paste0(dir, "/discovery/assemblytics_combined_results_with_component_group_small_", chr,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')
    write.table(assemblytics_combined_added_df_per_chr_big, paste0(dir, "/discovery/assemblytics_combined_results_with_component_group_big_", chr,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')
    write.table(singleton_per_chr_df, paste0(dir, "/discovery/singleton_", chr,".txt"), col.names = T, row.names = F, quote = F, sep = '\t')
    
}

