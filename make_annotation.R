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
## Notes: This script first annotates insertions then filters a set of high confident insertions
##
## ---------------------------------------------------------------------------------
##

# remove all previous variables
rm(list=ls())

# load up packages  
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(parallel))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))

# command line options
## ---------------------------------------------------------------------------------

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory"),
  make_option(c("-t", "--threads"), action="store", default=NA, type='numeric',
              help="number of threads"),
  make_option(c("-r", "--refFlat"), action="store", default=NA, type='character',
              help="path to refFlat file containing gene info"),
  make_option(c("-g", "--gwas"), action="store", default=NA, type='character',
              help="path to gwas variant file") 
  
)

opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
dir = opt$dir
threads = opt$threads
refFlat = opt$refFlat
gwas = opt$gwas

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

#refFlat = "../../../../KwokRaid02/karen/database/genome_annotation/03312019/refFlat.txt"
#dir="../../new_ref2/"
#gwas="../../../../KwokRaid02/karen/database/gwasCatalog.txt"

# ------------------------------------------------------------ Part 1: Annotation --------------------------------------------------------------
# read input sequences and databases
rep_seq = fread(paste0(dir, "/discovery/assemblytics_representative_seq.txt"), stringsAsFactors = F)

# read all the gene info and get exon coordinates
gene = fread(refFlat, stringsAsFactors = F)
gene = gene[,c(1:3,5,6,9,10,11)]
colnames(gene) = c("gene", "gene_id", "chr", "start", "end", "exon_count", "exon_start", "exon_end")

# read gwas variants
gwas = fread(gwas, stringsAsFactors = F, sep = '\t')
gwas = gwas[,c(2:4,11:12,15,18)]
colnames(gwas) = c("chr", "start", "end", "pheno", "populations", "gene", "Pval")

# read sample metadata
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source", "project")

# read trf output
trf = fread(paste0(dir ,"/discovery/final_fasta/trf_rep_seq_count.txt"), stringsAsFactors = F)
colnames(trf) = c("query", "len", "TRF_N_count", "TRF_N_perct")

# read gc content 
gc = fread(paste0(dir, "/discovery/final_fasta/gc_count.txt"), stringsAsFactors = F)
colnames(gc) = c("query", "gc_perct")

# read repeatmasker output
repeats = fread(paste0(dir, "/discovery/comp_repeat.txt"), stringsAsFactors = F)

# create exon file from the gene dataframe
exon_chr = rep(gene$chr, gene$exon_count)
exon_start = as.numeric(unlist(strsplit(gene$exon_start, ",")))
exon_end = as.numeric(unlist(strsplit(gene$exon_end, ",")))
exon_df = cbind.data.frame(exon_chr, exon_start, exon_end, stringsAsFactors = F)

# create grange objects for variables
rep_seq.gr = makeGRangesFromDataFrame(rep_seq, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end", ignore.strand = T)
gene.gr = makeGRangesFromDataFrame(gene)
exon.gr = makeGRangesFromDataFrame(exon_df, seqnames.field = "exon_chr", start.field = "exon_start", end.field = "exon_end")
gwas.gr = makeGRangesFromDataFrame(gwas)

# Find overlap
seq_gene_ov = findOverlaps(rep_seq.gr, gene.gr, type = "any")
seq_exon_ov = findOverlaps(rep_seq.gr, exon.gr, type = "any")
seq_gwas_ov = findOverlaps(rep_seq.gr, gwas.gr, type = "any")

rep_seq$gene = NA
rep_seq$gene[seq_gene_ov@from] = gene$gene[seq_gene_ov@to]

rep_seq$exon = 0
rep_seq$exon[seq_exon_ov@from] = 1

rep_seq$gwas = NA
rep_seq$gwas_populations = NA
rep_seq$gwas_Pval = NA
rep_seq$gwas[seq_gwas_ov@from] = gwas$pheno[seq_gwas_ov@to]
rep_seq$gwas_populations[seq_gwas_ov@from] = gwas$populations[seq_gwas_ov@to]
rep_seq$gwas_Pval[seq_gwas_ov@from] = gwas$Pval[seq_gwas_ov@to]

# Calculate distance to the nearest gene
min_dist_to_exon = NULL
min_dist_to_exon = mclapply(1:nrow(rep_seq), calc_min_dist_to_exon, mc.cores = threads)

# check if there're any error messages
job_err = min_dist_to_exon[1:nrow(rep_seq)]
stopifnot(length(which(grepl("Error", job_err)))==0)

rep_seq$min_dist_to_exon = unlist(min_dist_to_exon)

# Do a chi-sq analysis to see if the population distributions are uneven (for each insertion)
metadata_ethnicity_tbl = c(table(metadata$population)[c("AFR", "AMR", "EAS", "EUR", "SAS")])
metadata_ethnicity_tbl = ifelse(is.na(metadata_ethnicity_tbl), 0, metadata_ethnicity_tbl)
rep_seq$pop_chisq_P = apply(rep_seq[,c("sample_AFR", "sample_AMR", "sample_EAS", "sample_EUR", "sample_SAS")],1, function(x) chisq.test(rbind(x,metadata_ethnicity_tbl-x))$p.value)

# adjust p value for multiple testings
rep_seq$pop_chisq_P_adj = p.adjust(rep_seq$pop_chisq_P, method = "BH")

# add repeats annotation
rep_seq$repeat_class = repeats$repeat_class

#add trf annotation
rep_seq$key = paste0(rep_seq$assm_id,":", rep_seq$adjusted_assm_start,"-",rep_seq$adjusted_assm_end)
rep_seq = merge(rep_seq, trf, by.x = "key", by.y = "query")
rep_seq$TRF_N_count = ifelse(rep_seq$ngap==0, rep_seq$TRF_N_count, NA)
rep_seq$TRF_N_perct = ifelse(rep_seq$ngap==0, rep_seq$TRF_N_perct, NA)

# add gc%
rep_seq = merge(rep_seq, gc, by.x = "key", by.y = "query")

# remove key column
rep_seq = rep_seq[,2:ncol(rep_seq)]

# Write the annotated output
write.table(rep_seq, paste0(dir, "/discovery/assemblytics_representative_seq_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

# ------------------------------------------------ Part 2: Choose confident insertions -----------------------------------------------------

conf = rep_seq[rep_seq$ngap_boundaries!="yes" & rep_seq$edge_start<1 & rep_seq$edge_end<1 & rep_seq$ngap<=10,]
TR_discard = which(conf$TRF_N_perct>0.9 & conf$sample_perct<0.9)
small_variant_discard = which(conf$insert_size<50 & conf$sample_perct<0.5)
disc = c(unique(TR_discard, small_variant_discard))

if (length(disc)>=1){
  conf = conf[-disc,]
}

write.table(conf, paste0(dir, "/discovery/assemblytics_representative_seq_conf_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')









