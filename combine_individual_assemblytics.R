#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: combine_individual_assemblytics.R
##
## Purpose of script: Combine individual raw assemblytics files
##
## Author: Karen Wong
##
## Date Created: 06/07/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------
##
## Notes: This script takes the raw files in the WORKDIR/raw/* and combine into a single tab-delimited text file
##   
## ---------------------------------------------------------------------------------
##

# remove all previous variables
rm(list=ls())

# load up packages  
suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))

# command line options
## ---------------------------------------------------------------------------------

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory")

)
opt = parse_args(OptionParser(option_list=option_list))

## ---------------------------------------------------------------------------------

# prase user input
dir = opt$dir

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

# read metadata 
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("SAMPLE", "SEX", "FASTQ_DIR", "LONGRANGER_DIR", "ASSM_DIR", "BN_DIR", "ENZYME", "SUPERNOVA_VER", "ALT_NAME", "NUCMER_DIR", "POPULATION", "SRC")

# declare assemblytics path
assemblytics_path = paste0(dir, "/discovery/raw/")
assemblytics_list = list.files(assemblytics_path, pattern = "raw_results.txt")

# keep samples in the metadata sample list
file_list_sample = str_split_fixed(assemblytics_list, "_",2)[,1]
assemblytics_list = assemblytics_list[(file_list_sample %in% metadata$SAMPLE) | (file_list_sample %in% metadata$ALT_NAME)]

# read all the processed assemblytics files and append to a new variable
assemblytics_ALL = NULL
for (i in assemblytics_list){
  
  print(i)
  # Read individual processed assemblytics file
  assemblytics = read.table(paste0(assemblytics_path,i), stringsAsFactors = F, header = T)
  
  # Combine
  assemblytics_ALL = rbind.data.frame(assemblytics_ALL, assemblytics, stringsAsFactors = F)
  
}

# turn off scientific notation
options(scipen=999)
write.table(assemblytics_ALL, paste0(dir, "/discovery/assemblytics_combined_results.txt"), col.names = T, row.names = F, quote = F, sep = '\t')



