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
## Notes: This script annotates insertions 
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
suppressMessages(library(plyr))


# command line options
## ---------------------------------------------------------------------------------

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory"),
  make_option(c("-t", "--threads"), action="store", default=NA, type='numeric',
              help="number of threads"),
  make_option(c("-r", "--refFlat"), action="store", default=NA, type='character',
              help="path to refFlat file containing gene info"),
  make_option(c("-e", "--exon"), action="store", default=NA, type='character',
              help="path to coding exon, 5' and 3' UTR merged bed file"),
  make_option(c("-f", "--func"), action="store", default=NA, type='character',
              help="path to functional element bed file")
  
)

opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
dir = opt$dir
threads = opt$threads
refFlat = opt$refFlat
coding_utr = opt$exon
func = opt$func

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

#refFlat = "../../../../KwokRaid02/karen/database/genome_annotation/08022019/ncbiRefSeq_ucsc_track_all_fields.bed"
#coding_utr = "../../../../KwokRaid02/karen/database/genome_annotation/08022019/ncbiRefSeq_ucsc_track_coding_exon_3_5_UTR_merged.bed"
#dir="../"
#func = "../../../../KwokRaid02/karen/database/genome_annotation/08022019/refSeqFuncElems_all_fields.bed"

# ------------------------------------------------------------ Annotation --------------------------------------------------------------
#chr_list = c(1:22, "X", "Y")

# read gene info
refFlat = fread(refFlat, stringsAsFactors = F, select = c("chrom", "txStart", "txEnd", "name2"))
#colnames(refFlat) = c("gene_symbol", "chr", "start", "end")

# read coding exon bed file
coding_utr = fread(coding_utr, stringsAsFactors = F, select = c(1:3,7))
colnames(coding_utr) = c("chr", "start", "end", "type")

# read functional element annotation file
func = fread(func, stringsAsFactors = F, select = c("#chrom", "chromStart", "chromEnd", "name"))
colnames(func)[1] = "chr"

# read sample metadata
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source", "project")

# # read trf output
# trf = fread(paste0(dir ,"/discovery/final_fasta/trf_rep_seq_count.txt"), stringsAsFactors = F)
# colnames(trf) = c("query", "len", "TRF_N_count", "TRF_N_perct")
# 
# # read gc content 
# gc = fread(paste0(dir, "/discovery/final_fasta/gc_count.txt"), stringsAsFactors = F)
# colnames(gc) = c("query", "gc_perct")
# 
# # read repeatmasker output
# repeats = fread(paste0(dir, "/discovery/comp_repeat.txt"), stringsAsFactors = F)

# rep_seq_all = NULL
# for (i in 1:length(chr_list)){
# 
#     chr = chr_list[i]
#     print(chr)
    # read input sequences and databases
    rep_seq = fread(paste0(dir, "/discovery/all_raw.txt"), stringsAsFactors = F)
    colnames(rep_seq) = c("ref_chr","ref_start","ref_end","insert_size","strand","ref_gap_size","q_gap_size","assm_coords","assm_id","assm_start","assm_end","adjusted_coords","adjusted_assm_start","adjusted_assm_end","adjusted_insert_size","gap_ratio","BN_size","label_dist","BN_start","BN_end","BN_enzyme","BN_validated","ngap","ngap_boundaries","ngap_boundaries_size_left","ngap_boundaries_size_right","ngap_perct","sample","haplo","method","anchors_r","INS_id","component","type","idx","PB_only","edge_start","edge_end","avg_edge","edge_start2","edge_end2","avg_edge2","cluster_id","first_cluster_hap_size","first_cluster_hap_perct","first_cluster_sample_size","sample_count","sample_perct","sample_record","sample_AFR","sample_AMR","sample_EAS","sample_EUR","sample_SAS","sample_nonAFR","percent_AFR","percent_AMR","percent_EAS","percent_EUR","percent_SAS","percent_nonAFR")

    # if type==BN and BN_validated !=1 --> 2
    rep_seq$BN_validated = ifelse(rep_seq$BN_validated!=1 & rep_seq$type=="BN", 2, rep_seq$BN_validated)
    
    # create grange objects for variables
    rep_seq.gr = makeGRangesFromDataFrame(rep_seq, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end", ignore.strand = T)
    refFlat.gr = makeGRangesFromDataFrame(refFlat)
    coding_utr.gr = makeGRangesFromDataFrame(coding_utr)
    func.gr = makeGRangesFromDataFrame(func, seqnames.field = "chr", start.field = "chromStart", end.field = "chromEnd", ignore.strand = T)
    
    # Find overlap
    seq_refFlat_ov = findOverlaps(rep_seq.gr, refFlat.gr, type = "any")
    seq_coding_utr_ov = findOverlaps(rep_seq.gr, coding_utr.gr, type = "any")
    seq_func_ov = findOverlaps(rep_seq.gr, func.gr, type = "any")
    
    rep_seq$gene = NA
    rep_seq$gene[seq_refFlat_ov@from] = refFlat$name2[seq_refFlat_ov@to]
    
    rep_seq$genom_anno = "intron"
    rep_seq$genom_anno[seq_coding_utr_ov@from] = coding_utr$type[seq_coding_utr_ov@to]
    
    rep_seq$genom_anno = ifelse(is.na(rep_seq$gene) & rep_seq$genom_anno=="intron", "intergenic", rep_seq$genom_anno)
    
    rep_seq$func = NA
    rep_seq$func[seq_func_ov@from] = func$name[seq_func_ov@to]

    # rep_seq_all[[i]] = rep_seq
# } 
# 
# rep_seq_all_combined = ldply(rep_seq_all, data.frame)

# rep_seq$gwas = NA
# rep_seq$gwas_populations = NA
# rep_seq$gwas_Pval = NA
# rep_seq$gwas[seq_gwas_ov@from] = gwas$pheno[seq_gwas_ov@to]
# rep_seq$gwas_populations[seq_gwas_ov@from] = gwas$populations[seq_gwas_ov@to]
# rep_seq$gwas_Pval[seq_gwas_ov@from] = gwas$Pval[seq_gwas_ov@to]

# Calculate distance to the nearest gene
# min_dist_to_exon = NULL
# min_dist_to_exon = mclapply(1:nrow(rep_seq), calc_min_dist_to_exon, mc.cores = threads)
# 
# # check if there're any error messages
# job_err = min_dist_to_exon[1:nrow(rep_seq)]
# stopifnot(length(which(grepl("Error", job_err)))==0)
# 
# rep_seq$min_dist_to_exon = unlist(min_dist_to_exon)
# 
# # Do a chi-sq analysis to see if the population distributions are uneven (for each insertion)
# metadata_ethnicity_tbl = c(table(metadata$population)[c("AFR", "AMR", "EAS", "EUR", "SAS")])
# metadata_ethnicity_tbl = ifelse(is.na(metadata_ethnicity_tbl), 0, metadata_ethnicity_tbl)
# rep_seq$pop_chisq_P = apply(rep_seq[,c("sample_AFR", "sample_AMR", "sample_EAS", "sample_EUR", "sample_SAS")],1, function(x) chisq.test(rbind(x,metadata_ethnicity_tbl-x))$p.value)
# 
# # adjust p value for multiple testings
# rep_seq$pop_chisq_P_adj = p.adjust(rep_seq$pop_chisq_P, method = "BH")

# add repeats annotation
rep_seq[,c("repeatClass","SINE", "LINE", "LTR", "DNA", "small_RNA", "satellites","simpleRepeat", "lowComplex", "repeatEnriched_RM","trf", "gc")] = NA
 
#rep_seq = rep_seq[1:100,]

rep_seq_repeats = NULL
repeatAnnotate = function(i){   
#for (i in 1:nrow(rep_seq)){
  
    key = rep_seq$INS_id[i]
    try = tryCatch(fread(paste0(dir, "/discovery/repeatmasker/",key,"/",key,".fa.tbl"), skip = 9, nrows = 21, blank.lines.skip = T, fill = T, sep = '\t'), error=function(e) NULL)
    trf = tryCatch(fread(paste0(dir, "/discovery/repeatmasker/",key,"/trf_rep_seq_count.txt"), drop = 1:3), error=function(e) NULL)
    if (!is.null(trf)){
      rep_seq$trf[i] = trf$V4
    }
    gc = tryCatch(fread(paste0(dir, "/discovery/repeatmasker/",key,"/gc_count.txt"), drop = 1), error=function(e) NULL)
    if (!is.null(gc)){
      rep_seq$gc[i] = gc$V2
    }
    if (is.null(try)){
      rep_seq$repeatClass[i] = "none"
      rep_seq$repeatEnriched_RM[i] = 0
    }else{
      print(i)
      try = as.list(try[c(1,4,8,13,18,19,20,21)])
      names(try) = "classes"
      
      for (k in 1:8){
        
          col = k+65
          rep_seq[[i,col]] = rev(strsplit(try$classes[k], "\\s+")[[1]])[2]
      }
      rep_seq$repeatClass[i] = names(which.max(rep_seq[i,66:73]))
      rep_seq$repeatEnriched_RM[i] = ifelse(sum(as.numeric(rep_seq[i,66:73]))>50,1,0)
      
    }
    
    return(rep_seq[i])
}

rep_seq_repeats = mclapply(1:nrow(rep_seq), repeatAnnotate, mc.cores = threads)
rep_seq_repeats = ldply(rep_seq_repeats, data.frame)

job_err = res[1:length(rep_seq)]
stopifnot(length(which(grepl("Error", job_err)))==0)


        
    
# rep_seq$repeat_class = repeats$repeat_class
# 
# #add trf annotation
# rep_seq$key = paste0(rep_seq$assm_id,":", rep_seq$adjusted_assm_start,"-",rep_seq$adjusted_assm_end)
# rep_seq = merge(rep_seq, trf, by.x = "key", by.y = "query")
# rep_seq$TRF_N_count = ifelse(rep_seq$ngap==0, rep_seq$TRF_N_count, NA)
# rep_seq$TRF_N_perct = ifelse(rep_seq$ngap==0, rep_seq$TRF_N_perct, NA)
# 
# # add gc%
# rep_seq = merge(rep_seq, gc, by.x = "key", by.y = "query")
# 
# # remove key column
# rep_seq = rep_seq[,2:ncol(rep_seq)]

# Write the annotated output
write.table(rep_seq_repeats, paste0(dir, "/discovery/assemblytics_representative_seq_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')








