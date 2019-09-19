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
  make_option(c("-s", "--transcript"), action="store", default=NA, type='character',
              help="path to refseq transcript file"),
  make_option(c("-e", "--exon"), action="store", default=NA, type='character',
              help="path to coding exon, 5' and 3' UTR merged bed file"),
  make_option(c("-f", "--func"), action="store", default=NA, type='character',
              help="path to functional element bed file"),
  make_option(c("-b", "--build"), action="store", default=NA, type='character',
              help="path to ensemble regulatory build")
  
  
)

opt = parse_args(OptionParser(option_list=option_list))
## ---------------------------------------------------------------------------------

# Parse user input
dir = opt$dir
threads = opt$threads
refFlat = opt$refFlat
coding_utr = opt$exon
func = opt$func
build = opt$build

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

#transcript = "../../../../KwokRaid02/karen/database/genome_annotation/08022019/ncbiRefSeq_ucsc_track_all_transcripts.bed"
#coding_utr = "../../../../KwokRaid02/karen/database/genome_annotation/08022019/ncbiRefSeq_ucsc_track_coding_exon_3_5_UTR_merged_mod_gene_symbol.bed"
#dir="../"
#func = "../../../../KwokRaid02/karen/database/genome_annotation/08022019/refSeqFuncElems_all_fields.bed"
#build = "../../../../KwokRaid02/karen/database/genome_annotation/08022019/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff"

# ------------------------------------------------------------ Annotation --------------------------------------------------------------
chr_list = c(1:22, "X", "Y")

# read coding exon bed file
coding_utr = fread(coding_utr, stringsAsFactors = F)
coding = coding_utr[coding_utr$type=="coding_exon",]
utr5 = coding_utr[coding_utr$type=="5UTR",]
utr3 = coding_utr[coding_utr$type=="3UTR",]

# read transcript bed file (this is mainly to define intronic regions)
transcript = fread(transcript, stringsAsFactors = F)

# read functional element annotation file
func = fread(func, stringsAsFactors = F, select = c("#chrom", "chromStart", "chromEnd", "name"))
colnames(func)[1] = "chr"

# read ensemble regulatory build file
reg = fread(build, stringsAsFactors = F, select = c(1,3:5), sep = "\t")
colnames(reg) = c("chr", "reg_type", "start", "end")
reg = reg[reg$chr %in% chr_list,]
reg$chr = paste0("chr", reg$chr)

# Split reg into different reg_type
reg_CTCF = reg[reg$reg_type=="CTCF_binding_site",]
reg_enhancer = reg[reg$reg_type=="enhancer",]
reg_open_chr = reg[reg$reg_type=="open_chromatin_region",]
reg_promoter = reg[reg$reg_type=="promoter",]
reg_promoter_f = reg[reg$reg_type=="promoter_flanking_region",]
reg_TF = reg[reg$reg_type=="TF_binding_site",]

# read sample metadata
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population", "source", "project")

rep_seq = fread(paste0(dir, "/discovery/all_raw.txt"), stringsAsFactors = F)
colnames(rep_seq) = c("ref_chr","ref_start","ref_end","insert_size","strand","ref_gap_size","q_gap_size","assm_coords","assm_id","assm_start","assm_end","adjusted_coords","adjusted_assm_start","adjusted_assm_end","adjusted_insert_size","gap_ratio","BN_size","label_dist","BN_start","BN_end","BN_enzyme","BN_validated","ngap","ngap_boundaries","ngap_boundaries_size_left","ngap_boundaries_size_right","ngap_perct","sample","haplo","method","anchors_r","INS_id","component","type","idx","PB_only","edge_start","edge_end","avg_edge","edge_start2","edge_end2","avg_edge2","cluster_id","first_cluster_hap_size","first_cluster_hap_perct","first_cluster_sample_size","sample_count","sample_perct","sample_record","sample_AFR","sample_AMR","sample_EAS","sample_EUR","sample_SAS","sample_nonAFR","percent_AFR","percent_AMR","percent_EAS","percent_EUR","percent_SAS","percent_nonAFR")

# if type==BN and BN_validated !=1 --> 2
rep_seq$BN_validated = ifelse(rep_seq$BN_validated!=1 & rep_seq$type=="BN", 2, rep_seq$BN_validated)

# Need to split the dataframe into 2 (one df with pos ref gap size)
rep_seq_pos = rep_seq[rep_seq$ref_gap_size>5,]
rep_seq_pos$left_start = rep_seq_pos$ref_start-1
rep_seq_pos$left_end = rep_seq_pos$ref_start +1
rep_seq_pos$right_start =  rep_seq_pos$ref_end-1
rep_seq_pos$right_end = rep_seq_pos$ref_end +1
rep_seq_other = rep_seq[rep_seq$ref_gap_size<=5,]

# create grange objects for variables
rep_seq_pos_l.gr = makeGRangesFromDataFrame(rep_seq_pos, seqnames.field = "ref_chr", start.field = "left_start", end.field = "left_end", ignore.strand = T)
rep_seq_pos_r.gr = makeGRangesFromDataFrame(rep_seq_pos, seqnames.field = "ref_chr", start.field = "right_start", end.field = "right_end", ignore.strand = T)
rep_seq_other.gr = makeGRangesFromDataFrame(rep_seq_other, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end", ignore.strand = T)
func.gr = makeGRangesFromDataFrame(func, seqnames.field = "chr", start.field = "chromStart", end.field = "chromEnd", ignore.strand = T)

for (i in c("transcript", "coding", "utr5", "utr3", "reg_CTCF", "reg_enhancer", "reg_open_chr", "reg_promoter", "reg_promoter_f", "reg_TF")){
  assign(paste0(i,".gr"), makeGRangesFromDataFrame(eval(parse(text=i))))
}

# Find overlap
for (i in c("rep_seq_pos_l", "rep_seq_pos_r", "rep_seq_other")){
  
  for (f in c("transcript","coding", "utr5","utr3","func", "reg_CTCF", "reg_enhancer", "reg_open_chr", "reg_promoter", "reg_promoter_f","reg_TF")){
    
    assign(paste0(i,"_",f,"_ov"), findOverlaps(eval(parse(text=paste0(i,".gr"))), eval(parse(text=paste0(f,".gr"))), type = "any"))
    
  }

}

### Annotate the ref_seq_other
rep_seq_other[,c("gene","genom_anno", "func", "reg_CTCF", "reg_enhancer", "reg_open_chr", "reg_promoter", "reg_promoter_f", "reg_TF", "repeatClass","SINE", "LINE", "LTR", "DNA", "small_RNA", "satellites","simpleRepeat", "lowComplex", "repeatEnriched_RM","trf", "gc")] = 0
rep_seq_pos[,c("gene","genom_anno", "func", "reg_CTCF", "reg_enhancer", "reg_open_chr", "reg_promoter", "reg_promoter_f", "reg_TF", "repeatClass","SINE", "LINE", "LTR", "DNA", "small_RNA", "satellites","simpleRepeat", "lowComplex", "repeatEnriched_RM","trf", "gc")] = 0

# Genic annotations
rep_seq_other$gene[rep_seq_other_transcript_ov@from] = transcript$name2[rep_seq_other_transcript_ov@to]
rep_seq_other$genom_anno[rep_seq_other$gene==0] = "intergenic"
rep_seq_other$gene[rep_seq_other$gene==0] = NA
rep_seq_other$gene[rep_seq_other_utr5_ov@from] = utr5$gene_symbol[rep_seq_other_utr5_ov@to]
rep_seq_other$genom_anno[rep_seq_other_utr5_ov@from] = "5UTR"
rep_seq_other$gene[rep_seq_other_utr3_ov@from] = utr3$gene_symbol[rep_seq_other_utr3_ov@to]
rep_seq_other$genom_anno[rep_seq_other_utr3_ov@from] = "3UTR"
rep_seq_other$gene[rep_seq_other_coding_ov@from] = coding$gene_symbol[rep_seq_other_coding_ov@to]
rep_seq_other$genom_anno[rep_seq_other_coding_ov@from] = "coding"
# all other should be intronic
rep_seq_other$genom_anno[rep_seq_other$genom_anno==0 & !is.na(rep_seq_other$gene)] = "intron"
# functional annotations
rep_seq_other$func[rep_seq_other_func_ov@from] = func$name[rep_seq_other_func_ov@to]
rep_seq_other$reg_CTCF[rep_seq_other_reg_CTCF_ov@from] = 1
rep_seq_other$reg_enhancer[rep_seq_other_reg_enhancer_ov@from] = 1
rep_seq_other$reg_open_chr[rep_seq_other_reg_open_chr_ov@from] = 1
rep_seq_other$reg_promoter[rep_seq_other_reg_promoter_ov@from] = 1
rep_seq_other$reg_promoter_f[rep_seq_other_reg_promoter_f_ov@from] = 1
rep_seq_other$reg_TF[rep_seq_other_reg_TF_ov@from] = 1

#------------------------------------------------------------------------------------------------------------------------------------
### Now we need to deal with rep_seq_pos

# Genic annotations
rep_seq_pos$gene[rep_seq_pos_l_transcript_ov@from] = transcript$name2[rep_seq_pos_l_transcript_ov@to]
rep_seq_pos$gene[rep_seq_pos_r_transcript_ov@from] = transcript$name2[rep_seq_pos_r_transcript_ov@to]
rep_seq_pos$genom_anno[rep_seq_pos$gene==0] = "intergenic"
rep_seq_pos$gene[rep_seq_pos$gene==0] = NA

for (i in c("l", "r")){
  
  rep_seq_pos$gene[eval(parse(text=paste0("rep_seq_pos_",i,"_utr5_ov")))@from] = utr5$gene_symbol[eval(parse(text=paste0("rep_seq_pos_",i,"_utr5_ov")))@to]
  rep_seq_pos$genom_anno[eval(parse(text=paste0("rep_seq_pos_",i,"_utr5_ov")))@from] = "5UTR"
  rep_seq_pos$gene[eval(parse(text=paste0("rep_seq_pos_",i,"_utr3_ov")))@from] = utr3$gene_symbol[eval(parse(text=paste0("rep_seq_pos_",i,"_utr3_ov")))@to]
  rep_seq_pos$genom_anno[eval(parse(text=paste0("rep_seq_pos_",i,"_utr3_ov")))@from] = "3UTR"
  rep_seq_pos$gene[eval(parse(text=paste0("rep_seq_pos_",i,"_coding_ov")))@from] = coding$gene_symbol[eval(parse(text=paste0("rep_seq_pos_",i,"_coding_ov")))@to]
  rep_seq_pos$genom_anno[eval(parse(text=paste0("rep_seq_pos_",i,"_coding_ov")))@from] = "coding"
  # functional annotations
  rep_seq_pos$func[eval(parse(text=paste0("rep_seq_pos_",i,"_func_ov")))@from] = func$name[eval(parse(text=paste0("rep_seq_pos_",i,"_func_ov")))@to]
  rep_seq_pos$reg_CTCF[eval(parse(text=paste0("rep_seq_pos_",i,"_reg_CTCF_ov")))@from] = 1
  rep_seq_pos$reg_enhancer[eval(parse(text=paste0("rep_seq_pos_",i,"_reg_enhancer_ov")))@from] = 1
  rep_seq_pos$reg_open_chr[eval(parse(text=paste0("rep_seq_pos_",i,"_reg_open_chr_ov")))@from] = 1
  rep_seq_pos$reg_promoter[eval(parse(text=paste0("rep_seq_pos_",i,"_reg_promoter_ov")))@from] = 1
  rep_seq_pos$reg_promoter_f[eval(parse(text=paste0("rep_seq_pos_",i,"_reg_promoter_f_ov")))@from] = 1
  rep_seq_pos$reg_TF[eval(parse(text=paste0("rep_seq_pos_",i,"_reg_TF_ov")))@from] = 1
  
}
# all other should be intronic
rep_seq_pos$genom_anno[rep_seq_pos$genom_anno==0 & !is.na(rep_seq_pos$gene)] = "intron"

# Remove the 4 added columns 
keep = which(!colnames(rep_seq_pos) %in% c("left_start", "left_end", "right_start", "right_end"))
rep_seq_pos = rep_seq_pos[,..keep]

### Merge the two df back together
rep_seq = rbind.data.frame(rep_seq_pos, rep_seq_other, stringsAsFactors = F)

### Repeat annotations and GC analysis
rep_seq_repeats = NULL
repeatAnnotate = function(i){   
#for (i in 1:nrow(rep_seq)){
  
    key = rep_seq$INS_id[i]
    try = tryCatch(fread(paste0(dir, "/discovery/repeatmasker/",rep_seq$sample[i],"/",key,"/",key,".fa.tbl"), skip = 9, nrows = 21, blank.lines.skip = T, fill = T, sep = '\t'), error=function(e) NULL)
    gc = tryCatch(fread(paste0(dir, "/discovery/repeatmasker/",rep_seq$sample[i],"/",key,"/gc_count.txt")), error=function(e) NULL)
    if (!is.null(gc)){
      colnames(gc) = c("id", "len", "A", "C", "G", "T", "ambiguous", "ignore", "N", "CpG", "tv", "ts", "CpGpair")
      seq_N = as.numeric(gc$N)
      rep_seq$gc[i] = (gc$C + gc$G)/(gc$len-gc$N)
    }
    
    trf = tryCatch(fread(paste0(dir, "/discovery/repeatmasker/",rep_seq$sample[i],"/",key,"/trf_rep_seq_count.txt")), error=function(e) NULL)
    if (!is.null(trf)){
      colnames(trf) = c("id", "len", "A", "C", "G", "T", "ambiguous", "ignore", "N", "CpG", "tv", "ts", "CpGpair")
      rep_seq$trf[i] = (trf$N - seq_N)/(trf$len - seq_N)
    }
    
    if (is.null(try)){
      rep_seq$repeatClass[i] = "none"
      rep_seq$repeatEnriched_RM[i] = 0
    }else{
      print(i)
      try = as.list(try[c(1,4,8,13,18,19,20,21)])
      names(try) = "classes"
      
      for (k in 1:8){
        
          col = k+71
          rep_seq[[i,col]] = rev(strsplit(try$classes[k], "\\s+")[[1]])[2]
      }
      rep_seq$repeatClass[i] = names(which.max(rep_seq[i,72:79]))
      rep_seq$repeatEnriched_RM[i] = ifelse(sum(as.numeric(rep_seq[i,72:79]))>50,1,0)
      
    }
    
    return(rep_seq[i])
}

rep_seq_repeats = mclapply(1:nrow(rep_seq), repeatAnnotate, mc.cores = threads)
job_err = rep_seq_repeats[1:length(rep_seq)]
stopifnot(length(which(grepl("Error", job_err)))==0)

rep_seq_repeats = ldply(rep_seq_repeats, data.frame)

### Write the annotated output
write.table(rep_seq_repeats, paste0(dir, "/discovery/assemblytics_representative_seq_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')



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

#rep_seq = rep_seq[1:100,]



#df = fread("../discovery/assemblytics_representative_seq_annotated.txt", stringsAsFactors = F)




