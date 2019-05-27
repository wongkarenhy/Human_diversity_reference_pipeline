#!/usr/lib64/R/bin/Rscript
# This script annotates the representative insertion sequences
# This script has 2 parts
# Part1: Annotation
# --------------------------------------------------------------------------------------------------------------------------------------------
# ***************** 1. If the insertion is a singleton and without ngap, it has to have orthogonal validation by PB (+/-50bp) or BN (+/-700bp)
# ***************** 2. If a singleton has ngap, it is considered validated as long as the breakpoints overlap PB or BN
# ***************** 3. Not counting anything if the Ngap is at boundaries (hope to fill in the Ngap and redefine breakpoints with PB long reads)
# ***************** Sequences that we cannot use for sure: Singleton with no overlapping PB or BN entries***************************************
# --------------------------------------------------------------------------------------------------------------------------------------------
# Remove all previous variables
rm(list=ls())
######################################################################################################################
#command line options
library(optparse)

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory")
  
)
opt = parse_args(OptionParser(option_list=option_list))
################################################################################################################
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))

# Parse user input
dir = opt$dir

# ------------------------------------------------------------ Part 1: Annotation --------------------------------------------------------------
# Read input sequences and databases
seq = fread(paste0(dir, "/discovery/assemblytics_representative_seq.txt"), header = T, stringsAsFactors = F)
gene = fread("/media/KwokRaid02/karen/database/genome_annotation/03312019/refFlat.txt", stringsAsFactors = F)
exon = fread("/media/KwokRaid02/karen/database/genome_annotation/refgene_exon.bed", stringsAsFactors = F)
non_coding = fread("/media/KwokRaid02/karen/database/genome_annotation/01062019/lincRNAsTranscripts.txt", stringsAsFactors = F)
gencodeV29_nrRNA = fread("/media/KwokRaid02/karen/database/genome_annotation/09282018/gencode.v29.long_noncoding_RNAs_simplified.gtf", stringsAsFactors = F)
gencodeV29_pseudo = fread("/media/KwokRaid02/karen/database/genome_annotation/09282018/gencode.v29.2wayconspseudos_simplified.gtf", stringsAsFactors = F)
gwas = fread("/media/KwokRaid02/karen/database/gwasCatalog.txt", stringsAsFactors = F, sep = '\t')
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
repeatMasker = read.table(paste0(dir, "/discovery/final_fasta/repeats/assemblytics_representative_seq.fa.out"), stringsAsFactors = F, skip = 3, comment.char = "*")
trf = fread(paste0(dir ,"/discovery/final_fasta/trf_rep_seq_count.txt"), stringsAsFactors = F)

# Subset the proper columns
gene = gene[,c(1,3,5,6,10,11)]
colnames(gene) = c("gene", "chr", "start", "end", "exon_start", "exon_end")
exon = exon[,c(1:3)]
colnames(exon) = c("chr", "start", "end")
non_coding = non_coding[,c(2,3,5,6)]
colnames(non_coding) = c("name","chr", "start", "end")
colnames(gencodeV29_nrRNA) = c("chr", "start", "end")
colnames(gencodeV29_pseudo) = c("chr", "start", "end")
gwas = gwas[,c(2:4,11:12,15,18)]
colnames(gwas) = c("chr", "start", "end", "pheno", "populations", "gene", "Pval")
colnames(metadata) = c("sample", "sex", "fastq", "bam", "assembly", "BN", "enzyme", "supernova_ver", "alt_name", "nucmer", "population")
repeatMasker = repeatMasker[,c(5:7,10:11)]
colnames(repeatMasker) = c("name", "start", "end", "type", "group")
colnames(trf) = c("query", "len", "N_count", "N_perct")

# Overlap seq with refFlat
seq.gr = makeGRangesFromDataFrame(seq, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end", ignore.strand = T)
gene.gr = makeGRangesFromDataFrame(gene)
exon.gr = makeGRangesFromDataFrame(exon)
non_coding.gr = makeGRangesFromDataFrame(non_coding)
gencodeV29_nrRNA.gr = makeGRangesFromDataFrame(gencodeV29_nrRNA)
gencodeV29_pseudo.gr = makeGRangesFromDataFrame(gencodeV29_pseudo)
gwas.gr = makeGRangesFromDataFrame(gwas)

# Find overlap
seq_gene_ov = findOverlaps(seq.gr, gene.gr, type = "any")
seq_exon_ov = findOverlaps(seq.gr, exon.gr, type = "any")
seq_non_coding_ov = findOverlaps(seq.gr, non_coding.gr, type = "any")
seq_gencodeV29_nrRNA_ov = findOverlaps(seq.gr, gencodeV29_nrRNA.gr, type = "any")
seq_gencodeV29_pseudo_ov = findOverlaps(seq.gr, gencodeV29_pseudo.gr, type = "any")
seq_gwas_ov = findOverlaps(seq.gr, gwas.gr, type = "any")

seq$gene = NA
seq$gene[seq_gene_ov@from] = gene$gene[seq_gene_ov@to]

seq$exon = 0
seq$exon[seq_exon_ov@from] = 1

seq$non_coding = 0
seq$non_coding[seq_non_coding_ov@from] = 1

seq$gencodeV29_nrRNA = 0
seq$gencodeV29_nrRNA[seq_gencodeV29_nrRNA_ov@from] = 1

seq$gencodeV29_pseudo = 0
seq$gencodeV29_pseudo[seq_gencodeV29_pseudo_ov@from] = 1

seq$gwas = NA
seq$gwas_populations = NA
seq$gwas_Pas = NA
seq$gwas[seq_gwas_ov@from] = gwas$pheno[seq_gwas_ov@to]
seq$gwas_populations[seq_gwas_ov@from] = gwas$populations[seq_gwas_ov@to]
seq$gwas_Pas[seq_gwas_ov@from] = gwas$Pval[seq_gwas_ov@to]

# Calculate distance to the nearest gene
seq$min_dist_to_exon = NA
for (i in 1:nrow(seq)){
  
  breakpoint1 = seq$ref_start[i]
  breakpoint2 = seq$ref_end[i]
  specific_chr = seq$ref_chr[i]
  # Subset the exon file
  exon_subset = exon[which(exon$chr==specific_chr),]
  
  # Compute the minimum distance from either breakpoint to the nearest exon
  seq$min_dist_to_exon[i] = min(abs(exon_subset$start - breakpoint1), abs(exon_subset$end - breakpoint1), abs(exon_subset$start - breakpoint2), abs(exon_subset$end - breakpoint2))
  
}

# Do a chi-sq analysis to see if the population distributions are uneven (for each insertion)
pop_total_count = as.data.frame(table(metadata$population))
seq$pop_chisq_P = apply(seq[,c("sample_AFR", "sample_AMR", "sample_EAS", "sample_EUR", "sample_SAS")],1, function(x) chisq.test(rbind(x,pop_total_count$Freq))$p.value)

# Adjust p value for multiple testings
seq$pop_chisq_P_adj = p.adjust(seq$pop_chisq_P, method = "BH")

# Annotate repeats
repeatMasker$key_id = str_split_fixed(repeatMasker$name, ":|-", 3)[,1]
repeatMasker$key_start = str_split_fixed(repeatMasker$name, ":|-", 3)[,2]

repeatMasker = merge(repeatMasker, seq[,c("assm_id", "adjusted_assm_start", "insert_size", "ref_gap_size", "q_gap_size")], all.x = T, by.x = c("key_id", "key_start"), by.y = c("assm_id", "adjusted_assm_start"))
repeatMasker$total_size = ifelse(repeatMasker$ref_gap_size<0, repeatMasker$insert_size-(repeatMasker$ref_gap_size*2), repeatMasker$insert_size)
repeatMasker$repeat_frac = (repeatMasker$end - repeatMasker$start)/repeatMasker$total_size

repeats = NULL
for (i in unique(repeatMasker$name)){
  df = repeatMasker[repeatMasker$name==i,]
  if (nrow(df)==1){
      df$repeat_group = df$group
      df$repeat_type = df$type
      df$repeat_len = df$end - df$start
  } else{
      df$repeat_group = df$group[which.max(df$repeat_frac)]
      df$repeat_type = df$type[which.max(df$repeat_frac)]
      df$repeat_len = df$end[which.max(df$repeat_frac)] - df$start[which.max(df$repeat_frac)]
      df$repeat_frac = sum(df$repeat_frac)
      df$repeat_frac = ifelse(df$repeat_frac>1, 1, df$repeat_frac)
  }
  df = df[1,c("key_id", "key_start", "name", "repeat_frac", "repeat_type","repeat_group", "repeat_len")]
  repeats = rbind.data.frame(repeats, df, stringsAsFactors = F)
}

# Make sure the keys are of the same class
repeats$key_id = as.character(repeats$key_id)
repeats$key_start = as.numeric(repeats$key_start)
seq$assm_id = as.character(seq$assm_id)
seq$adjusted_assm_start = as.numeric(seq$adjusted_assm_start)
# Merge repeats dataframe back to seq
seq_merged = merge(seq, repeats, all.x = T, by.x = c("assm_id", "adjusted_assm_start"), by.y = c("key_id", "key_start"))

# Merge TRF percent for each entry
trf$key_id = as.character(str_split_fixed(trf$query, ":|-", 3)[,1])
trf$key_start = as.numeric(str_split_fixed(trf$query, ":|-", 3)[,2])
seq_merged = merge(seq_merged, trf, all.x = T, by.x = c("assm_id", "adjusted_assm_start"), by.y = c("key_id", "key_start"))

# Remove unucessary columns
seq_merged = seq_merged[,c(3:17,1,18:23,2,24:74)]

# Write the annotated output
write.table(seq_merged, paste0(dir, "/discovery/assemblytics_representative_seq_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')

# Define a conf list for constructing the reference
conf = seq_merged[seq_merged$ngap_boundaries=="no" & seq_merged$N_perct<0.9,]
conf = conf[(conf$sample_count>1 | (conf$sample_count==1 & conf$ngap==0 & conf$PB_validated==1 & conf$haplo!="unphased")), ]
write.table(conf, paste0(dir, "/discovery/assemblytics_representative_seq_conf_annotated.txt"), col.names = T, row.names = F, quote = F, sep = '\t')



# ------------------------------------------------------------ Part 2: Ngap filling -----------------------------------------------------------

# Get a list of potentially resolvable ngaps
# Resolvable 
#ngap = seq_merged[seq_merged$ngap!=0 & !is.na(seq_merged$PB_id),]
#write.table(ngap, "../discovery/assemblytics_representative_seq_resolvable_ngap.txt", col.names = T, row.names = F, quote = F, sep = '\t')

# Get a list of all sequences with Ngaps
#seq_merged = fread("../discovery/assemblytics_representative_seq_annotated.txt", stringsAsFactors = F, sep = '\t')
#all_ngap = seq_merged[seq_merged$ngap!=0,]
#write.table(all_ngap, "../discovery/assemblytics_representative_seq_all_ngap.txt", col.names = T, row.names = F, quote = F, sep = '\t')
#send_eleanor = all_ngap[,c(1:8,22,23)]
#send_eleanor$haplo = ifelse(send_eleanor$haplo=="2.1;2.2", "2.1", send_eleanor$haplo)
#write.table(send_eleanor, "../discovery/assemblytics_representative_seq_send_eleanor.txt", col.names = T, row.names = F, quote = F, sep = '\t')


# ------------------------------------------------------------ Part 3: Manual plotting -----------------------------------------------------------
# Subset PB/BN varified sequences and find AluY and AluS contribution
#verified = seq_merged[(seq_merged$PB_validated==1 | seq_merged$validated==1) & seq_merged$ngap==0,]
#verified_alu = verified[grep("Alu", verified$repeat_type),]
#verified_alu$alu = substr(verified_alu$repeat_type, 1,4)
#verified_alu$sample_count = as.factor(verified_alu$sample_count)
#verified_aluY = verified_alu[verified_alu$alu=="AluY",]
#verified_aluS = verified_alu[verified_alu$alu=="AluS",]

#ggplot(verified_alu_m, aes(fill=alu, value)) + geom_density(alpha=0.3) + scale_y_sqrt()
#mean_aluY = mean(verified_aluY$insert_size)
#mean_aluS = mean(verified_aluS$insert_size)

#ggplot(verified_aluY, aes(insert_size)) + geom_histogram(bins = 100) + scale_x_continuous(breaks=seq(0,7000,500), limits = c(0,7000)) +
#  geom_vline(xintercept = mean_aluY, color = "red", linetype="dashed") + xlab("Insertion Size (bp)") + ylab("Count") +
#  ggtitle("AluY Associated Insertions") + geom_text(x=mean_aluY+400, y=900, label=round(mean_aluY,2), color = "red")

#ggplot(verified_aluS, aes(insert_size)) + geom_histogram(bins = 100)+ scale_x_continuous(breaks=seq(0,7000,500), limits = c(0,7000)) +
#  geom_vline(xintercept = mean_aluS, color = "red", linetype="dashed") + xlab("Insertion Size (bp)") + ylab("Count") +
#  ggtitle("AluS Associated Insertions") + geom_text(x=mean_aluS+400, y=60, label=round(mean_aluS,2), color = "red")










# plot ref_gap_size
#hist(seq$ref_gap_size, prob = T, nclass = 5000, xlim = c(-500,500), 
#     border = "white", col = "black", xlab = "Reference gap size (bp)", main = NA, cex.lab = 1.3)
#abline(v=mean(seq$ref_gap_size), col = "red", lwd = 2, lty = 3)

# Plot min_dist_to_exon
#hist((seq$min_dist_to_exon/1000), nclass = 5000, probability = T, col = "black", border = "white", 
#     xlab = "Min distance to the nearest exon (kb)", main = NA, cex.lab = 1.3, xlim = c(0,100))
#lines(density(seq$min_dist_to_exon/1000), col = "red", lwd = 2)

# Plot population frequency




#anno = fread("../discovery/assemblytics_representative_seq_annotated.txt", header = T, stringsAsFactors = F)



