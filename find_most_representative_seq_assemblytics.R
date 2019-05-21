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
              help="EEE's vcf file"),
  make_option(c("-t", "--threads"), action="store", default=NA, type='numeric',
              help="number of threads")
  
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
threads = opt$threads

#source("insertion_filtering_functions.R")
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

assemblytics_combined = fread(paste0(dir, "/discovery/assemblytics_combined_results.txt"), stringsAsFactors = F, header = T)

# Read EEE's vcf file
vcf = read.table(vcf, stringsAsFactors = F)
colnames(vcf) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "CHM1", "CHM13", "HG00514", "HG00733", "NA12878", "HG02818", "NA19434", "HG01352", "HG02059", "NA12878", "HG04217", "HG02106", "HG00268", "HX1")
  
#assemblytics_combined_per_chr = assemblytics_combined[assemblytics_combined$ref_chr==CHR,] # just to make testing faster
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source")

# Convert the validated column to a numeric factor
assemblytics_combined$validated[is.na(assemblytics_combined$validated)] = -1
assemblytics_combined$validated[assemblytics_combined$validated=="yes"] = 1
assemblytics_combined$validated[assemblytics_combined$validated=="no"] = 0

# Find total non-redundant sequences 
assemblytics.gr = makeGRangesFromDataFrame(assemblytics_combined, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
assemblytics_assemblytics_ov_self = findOverlaps(assemblytics.gr, assemblytics.gr, type = "any", ignore.strand=T)

# Find a representative sequence in each group
## Representative is defined by the sequence with the min difference between insert_size and BN_size, and possibly no ngap
## Here we're only using alignments that are also identified by BN and at least 500+/- bp

edges = as.data.frame(assemblytics_assemblytics_ov_self)
net = graph_from_data_frame(d=edges, directed = F)
membership = components(net)$membership
components_count = count_components(net)
assemblytics_combined$component = membership # Add the component number to each INS

options(scipen=999)
write.table(assemblytics_combined, paste0(dir, "/discovery/assemblytics_combined_per_chr_results_with_component_group.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

CHR_LIST = paste0("chr",c(1:22,"X","Y"))
#CHR_LIST="chr1"
representative_seq_all = NULL
for (CHR in CHR_LIST){
    
  # Subset the file
  assemblytics_combined_per_chr <<- assemblytics_combined[assemblytics_combined$ref_chr==CHR,]
  # Search for the representative sequence for this non-redundant dataset
  # Rate-limiting step (make this run in parallel)
  representative_seq <<- NULL
  comp_list <<- unique(assemblytics_combined_per_chr$component)

  representative_seq = mclapply(1:length(comp_list), findRepresentativeSeq, mc.cores = threads)
  
  #head(representative_seq)
  # Store result in a dataframe
  representative_seq = data.frame(matrix(unlist(representative_seq), nrow=length(representative_seq), byrow=T), stringsAsFactors = F)
  header = c(colnames(assemblytics_combined_per_chr), "population", "sample_record", "sample_count", "sample_AFR", "sample_AMR", "sample_EAS", "sample_EUR", "sample_SAS", "sample_nonAFR", "percent_AFR", "percent_AMR", "percent_EAS", "percent_EUR", "percent_SAS", "percent_nonAFR", "edge_start", "edge_end")
  colnames(representative_seq) = header
  
  representative_seq_all = rbind.data.frame(representative_seq_all, representative_seq, stringsAsFactors = F)
  print(paste0(CHR,"done!"))
}  

write.table(representative_seq_all, paste0(dir, "/discovery/assemblytics_representative_seq_raw.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

# Re-read this file so column classes can be assigned correctly
representative_seq = fread(paste0(dir, "/discovery/assemblytics_representative_seq_raw.txt"), stringsAsFactors = F, header = T)

# Need to adjust the start coordinates
representative_seq$ref_start = representative_seq$ref_start - 1

# Add contig lenth to avoid "truncated seq" problem in the next step 
# Read file containing contig length info
fai = fread(paste0(dir,"/discovery/tmp_idx.txt"), stringsAsFactors = F)
colnames(fai) = c("assm_id", "contig_length", "sample", "haplo")

representative_seq = merge(representative_seq, fai, all.x = T, by.x = c("assm_id", "sample", "haplo"), by.y = c("assm_id", "sample", "haplo"))
# Remove if adjusted_assm_start or adjusted_assm_end is outside the contig length
representative_seq = representative_seq[representative_seq$adjusted_assm_start>0 & (representative_seq$adjusted_assm_end < representative_seq$contig_length),]
# Remove contig_length
representative_seq = representative_seq[,c(4:18,1,19:31,2:3,32:48)]

##------------------------------------------------PacBio comparison----------------------------------------------
# Make another GRange object for PB comparison
representative_seq.gr = makeGRangesFromDataFrame(representative_seq, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")

# Compare with EEE's SV datasets and add an additional column

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
representative_seq$PB_validated = ifelse(abs(representative_seq$PB-representative_seq$insert_size) <= 50, 1,0)

representative_seq$adjusted_assm_start = as.integer(representative_seq$adjusted_assm_start)
# Remove if adjusted_assm_start is still greater than adjusted_assm_end
discard = which(representative_seq$adjusted_assm_start>representative_seq$adjusted_assm_end)
representative_seq = representative_seq[-discard,]

# Turn off scientific notation
write.table(representative_seq, paste0(dir, "/discovery/assemblytics_representative_seq.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

#representative_seq = fread("../discovery/assemblytics_representative_seq.txt", stringsAsFactors = F, header = T)














#representative_seq = fread("../discovery/assemblytics_representative_seq.txt", stringsAsFactors = F, header = T)
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


