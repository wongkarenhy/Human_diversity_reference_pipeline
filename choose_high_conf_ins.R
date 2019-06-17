#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: make_annotation.R
##
## Purpose of script: Annotate representative insertion sequences
##
## Author: Karen Wong
##
## Date Created: 06/11/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------
##
## Notes: This script filters for a set of high confident insertions
##
## ---------------------------------------------------------------------------------
##

# remove all previous variables
rm(list=ls())

# load up packages  
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

# ------------------------------------------------  Choose confident insertions -----------------------------------------------------
# read all annotated insertions
rep_seq = fread(paste0(dir, "/discovery/assemblytics_representative_seq_annotated.txt"), stringsAsFactors = F)

conf = rep_seq[rep_seq$ngap_boundaries!="yes" & rep_seq$edge_start<1 & rep_seq$edge_end<1 & rep_seq$ngap<=10,]
TR_discard = which(conf$TRF_N_perct>0.9 & conf$sample_perct<0.9)
small_variant_discard = which(conf$insert_size<50 & conf$sample_perct<0.5)
disc = c(unique(TR_discard, small_variant_discard))

if (length(disc)>=1){
  conf = conf[-disc,]
}

# Make sure variants are at least 50bp apart
#conf = conf[order(conf$ref_chr, conf$ref_start),]
#which(diff(conf$ref_start)<=50)

write.table(conf, paste0(dir, "/discovery/assemblytics_representative_seq_conf_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

