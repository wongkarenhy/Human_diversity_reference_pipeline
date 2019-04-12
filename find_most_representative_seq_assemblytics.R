#!/usr/lib64/R/bin/Rscript
# Remove all previous variables
rm(list=ls())
######################################################################################################################
#command line options
library(optparse)

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory"),
  make_option(c("-v", "--vcf"), action="store", default=NA, type='character',
              help="EEE's vcf file")
  
)
opt = parse_args(OptionParser(option_list=option_list))
################################################################################################################
suppressMessages(library(GenomicRanges))
suppressMessages(library(igraph))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

# Parse user input
dir = opt$dir
vcf = opt$vcf

#source("insertion_filtering_functions.R")
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

#lastz_combined = fread("../discovery/assemblytics_pseudohap_combined_results.txt", stringsAsFactors = F, header = T)
#metadata = read.table("../ALL_sample_metadata.txt", stringsAsFactors = F)
lastz_combined = fread(paste0(dir, "/discovery/assemblytics_pseudohap_combined_results.txt"), stringsAsFactors = F, header = T)
metadata = read.table(paste0(dir, "/ALL_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population")

# Read EEE's vcf file
vcf = read.table(vcf, stringsAsFactors = F)
colnames(vcf) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "CHM1", "CHM13", "HG00514", "HG00733", "NA12878", "HG02818", "NA19434", "HG01352", "HG02059", "NA12878", "HG04217", "HG02106", "HG00268", "HX1")

# Convert the validated column to a numeric factor
lastz_combined$validated[is.na(lastz_combined$validated)] = -1
lastz_combined$validated[lastz_combined$validated=="yes"] = 1
lastz_combined$validated[lastz_combined$validated=="no"] = 0

# Find total non-redundant sequences 
lastz.gr = makeGRangesFromDataFrame(lastz_combined, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
lastz_lastz_ov_self = findOverlaps(lastz.gr, lastz.gr, type = "any", ignore.strand=T)

# Find a representative sequence in each group
## Representative is defined by the sequence with the min difference between insert_size and BN_size, and possibly no ngap
## Here we're only using alignments that are also identified by BN and at least 500+/- bp

edges = as.data.frame(lastz_lastz_ov_self)
net = graph_from_data_frame(d=edges, directed = F)
membership = components(net)$membership
components_count = count_components(net)

# Search for the representative sequence for this non-redundant dataset
representative_seq = findRepresentativeSeq(membership, components_count, lastz_combined, metadata)

# In the representative set, insertions >2kb has to be confirmed by BN (size doesn't matter due to Ngap)
# representative_seq = representative_seq[!(representative_seq$insert_size>=2000 & representative_seq$BN_size==0),]

# Need to adjust a few things to make downstream analysis easier
representative_seq$adjusted_assm_start = ifelse(representative_seq$ref_gap_size <0,representative_seq$assm_start+representative_seq$ref_gap_size, representative_seq$assm_start)
representative_seq$adjusted_assm_end = ifelse(representative_seq$ref_gap_size <0,representative_seq$assm_end-representative_seq$ref_gap_size, representative_seq$assm_end)
representative_seq$ref_start = representative_seq$ref_start - 1

# adjusted_assm_start has to be greater than 0
representative_seq = representative_seq[representative_seq$adjusted_assm_start>0,]
representative_seq.gr = makeGRangesFromDataFrame(representative_seq, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")

# Also compare with EEE's datasets and add an additional column

# Break the INFO column
vcf$type = substr(str_split_fixed(vcf$INFO, ";", 40)[,1],8,10)
vcf$insert_size = as.numeric(substr(str_split_fixed(vcf$INFO, ";", 40)[,2],7, nchar(str_split_fixed(vcf$INFO, ";", 40)[,2])))
vcf$END = substr(str_split_fixed(vcf$INFO, ";", 40)[,3],5, nchar(str_split_fixed(vcf$INFO, ";", 40)[,3]))
vcf$repeat_type = substr(str_split_fixed(vcf$INFO, ";", 40)[,15],13, nchar(str_split_fixed(vcf$INFO, ";", 40)[,15]))
vcf$n_count = lapply(str_split(str_split_fixed(vcf$INFO, ";", 40)[,5],","), function(x) length(x))
vcf$contig_name = substr(str_split_fixed(vcf$INFO, ";", 40)[,12],8, nchar(str_split_fixed(vcf$INFO, ";", 40)[,12]))

# Keep insertions only
vcf_INS = vcf[vcf$type=="INS",]
vcf_INS$ref_start = as.numeric(vcf_INS$POS)
vcf_INS$ref_end = as.numeric(vcf_INS$END)+vcf_INS$insert_size
vcf_INS$n_count = unlist(vcf_INS$n_count)

# Make a GRange object
vcf_INS.gr = makeGRangesFromDataFrame(vcf_INS, seqnames.field = "CHROM", start.field = "ref_start", end.field = "ref_end")

# Overlap PB and our datasets
PB_rep_ov = findOverlaps(representative_seq.gr, vcf_INS.gr, type = "any")

# Add a new column describing overlaps with PB INS
# Need to loop through the results because some queries overlap more than 1 PB entry
representative_seq$PB = 0
representative_seq$PB_freq = 0
representative_seq$PB_qual = NA
representative_seq$PB_id = NA
representative_seq$PB_contig_name = NA

for (i in unique(PB_rep_ov@from)){
    df = PB_rep_ov[PB_rep_ov@from==i,]
    if (length(df)==1){
        representative_seq$PB[i] = vcf_INS$insert_size[df@to]
        representative_seq$PB_id[i] = vcf_INS$ID[df@to]
        representative_seq$PB_contig_name[i] = vcf_INS$contig_name[df@to]
        representative_seq$PB_freq[i] = vcf_INS$n_count[df@to]
        representative_seq$PB_qual[i] = vcf_INS$QUAL[df@to]
        
        next
    } else {
        expected_INS_size = representative_seq$insert_size[i]
        min_diff = which.min(abs(vcf_INS$insert_size[df@to]-expected_INS_size))
        representative_seq$PB[i] = vcf_INS$insert_size[df@to[min_diff]]
        representative_seq$PB_id[i] = paste(vcf_INS$ID[df@to], collapse = ";")
        representative_seq$PB_contig_name[i] = paste(vcf_INS$contig_name[df@to], collapse = ";")
        representative_seq$PB_freq[i] = paste(vcf_INS$n_count[df@to], collapse = ";")
        representative_seq$PB_qual[i] = paste(vcf_INS$QUAL[df@to], collapse = ";")
        
    }
}

representative_seq$PB = as.numeric(representative_seq$PB)
# PB validation is +/- 50bp
representative_seq$PB_validated = ifelse(abs(representative_seq$PB-representative_seq$insert_size)<=50, 1,0)

write.table(representative_seq, paste0(dir, "/discovery/assemblytics_representative_seq.txt"), col.names = T, row.names = F, quote = F, sep = '\t')












# Define high confident set
#representative_seq_conf = representative_seq[representative_seq$ngap_boundaries=="no" & representative_seq$sample_count>1 ,] 
#singleton_conf_no_ngap = representative_seq[representative_seq$ngap_boundaries=="no" & representative_seq$sample_count==1 & representative_seq$ngap==0 & (representative_seq$validated==1|representative_seq$PB_validated==1),]
#singleton_conf_ngap = representative_seq[representative_seq$ngap_boundaries=="no" & representative_seq$sample_count==1 & representative_seq$ngap!=0 & ((representative_seq$BN_size!=0 & representative_seq$insert_size>=1000) | representative_seq$PB!=0),]
#representative_seq_conf = rbind.data.frame(representative_seq_conf, singleton_conf_no_ngap, singleton_conf_ngap, stringsAsFactors = F)

# Define validated sequences
#representative_seq_val = representative_seq_conf[representative_seq_conf$validated==1,]

#write.table(representative_seq_conf, "../discovery/assemblytics_representative_seq_conf.txt", col.names = T, row.names = F, quote = F, sep = '\t')
#write.table(representative_seq_val, "../discovery/assemblytics_representative_seq_validated.txt", col.names = T, row.names = F, quote = F, sep = '\t')

# Get a list of singletons
#singletons = representative_seq[representative_seq$sample_count==1 & representative_seq$validated!=1,]
#write.table(representative_seq, "../discovery/assemblytics_singletons.txt", col.names = T, row.names = F, quote = F, sep = '\t')

#sum(representative_seq$insert_size)
#nrow(representative_seq)

#sum(representative_seq_conf$insert_size)
#nrow(representative_seq_conf)

#rep_seq = fread("../discovery/assemblytics_representative_seq.txt", stringsAsFactors = F, header = T)


