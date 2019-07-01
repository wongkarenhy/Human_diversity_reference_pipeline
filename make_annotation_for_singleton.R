#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: make_annotation_for_singleton.R
##
## Purpose of script: Annotate singleton insertions
##
## Author: Karen Wong
##
## Date Created: 06/30/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------
##
## Notes: This script annotates all singletons 
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
              help="path to gwas variant file"),
  make_option(c("-c", "--chr"), action="store", default=NA, type='character',
              help="list of chromosomes")  
  
)

opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
dir = opt$dir
threads = opt$threads
refFlat = opt$refFlat
gwas = opt$gwas
chr_list = opt$chr

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

#refFlat = "../../../../KwokRaid02/karen/database/genome_annotation/03312019/refFlat.txt"
#dir="../"
#gwas="../../../../KwokRaid02/karen/database/gwasCatalog.txt"

# ------------------------------------------------------------ Annotation --------------------------------------------------------------
# read input sequences and databases
rep_seq = NULL
seqtk = NULL
for (i in 1:length(chr_list)){
  
  chr = chr_list[i]
  # Read individual processed assemblytics file
  rep_seq[[i]] = fread(paste0(dir, "/singleton/singleton_chr", chr,".txt"), stringsAsFactors = F)
  seqtk[[i]] = fread(paste0(dir, "/singleton/seqtk_comp_chr", chr, '.txt'), stringsAsFactors = F)
  
}

rep_seq = ldply(rep_seq, data.frame)
seqtk_df = ldply(seqtk, data.frame)

colnames(seqtk_df) = c("name", "length", "A", "C" ,"G", "T", "two", "three", "N", "CpG", "tv", "ts", "CpG_ts")
seqtk_df$total_unamb_unmask = rowSums(seqtk_df[,c("A", "C", "G", "T")])

# read all the gene info and get exon coordinates
gene = fread(refFlat, stringsAsFactors = F)
gene = gene[,c(1:3,5,6,9,10,11)]
colnames(gene) = c("gene", "gene_id", "chr", "start", "end", "exon_count", "exon_start", "exon_end")

# read gwas variants
gwas = fread(gwas, stringsAsFactors = F, sep = '\t')
gwas = gwas[,c(2:4,11:12,15,18)]
colnames(gwas) = c("chr", "start", "end", "pheno", "populations", "gene", "Pval")

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

# split seqtk_df$name into separate columns for merging in the next step
seqtk_df$assm_id = str_split_fixed(seqtk_df$name,":|-", 3)[,1]
seqtk_df$adjusted_assm_start = str_split_fixed(seqtk_df$name,":|-", 3)[,2]
seqtk_df$adjusted_assm_end = str_split_fixed(seqtk_df$name,":|-", 3)[,3]
singleton_df = merge(rep_seq, seqtk_df, by = c("assm_id", "adjusted_assm_start", "adjusted_assm_end"))

# Write the annotated output
write.table(singleton_df, paste0(dir, "/discovery/singleton_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')









