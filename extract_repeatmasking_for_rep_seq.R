#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: extract_repeatmasking_for_rep_seq.R
##
## Purpose of script: Pick the dominant repeat class per insertion
##
## Author: Karen Wong
##
## Date Created: 06/10/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------

# remove all previous variables
rm(list=ls())

# load up packages  
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(parallel))

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

# read representative seq dataframe
rep_seq = fread(paste0(dir, "/discovery/assemblytics_representative_seq.txt"), stringsAsFactors = F)

repeat_class = NULL
processRepeatMasker = function(i) {
#for (i in 1:nrow(rep_seq)){
    
    print(i)
  
    comp = rep_seq$component[i]
    id = rep_seq$INS_id[i]
    
    repeatRes <- try(read.table(paste0(dir, "/tmp/",comp,"/",comp,".fa.out"), stringsAsFactors = F, skip = 3, comment.char = "*"), silent = T)

    if (inherits(repeatRes, 'try-error')){
        repeat_class = NA
    
    #repeatRes = read.table(paste0(dir, "/tmp/",comp,"/",comp,".fa.out"), stringsAsFactors = F, skip = 3, comment.char = "*")
    
    #repeatRes = fread(paste0(dir, "/tmp/",comp,"/",comp,".fa.out"), stringsAsFactors = F, fill=TRUE, sep=' ')
    
    #if (nrow(repeatRes)==0){
    #     repeat_class = NA
    } else{
      
      # add a header
      colnames(repeatRes) = c("bit", "perc_div", "perc_del", "perc_ins", "query", "begin", "end", "q_left", "strand", "matching_repeat", "repeat_class", "begin_rep", "end_rep", "repet_left", "ID")
      
      target = repeatRes[grepl(id, repeatRes$query),]
      
      if (nrow(target)==0){
        
          repeat_class = NA
          
      } else{
          
          target$repeat_len = as.numeric(target$end) - as.numeric(target$begin)
          repeat_class = names(which.max(tapply(target$repeat_len, target$repeat_class, sum)))
        
      }
      
      
    }
    
    return(repeat_class)
}
repeat_class = mclapply(1:nrow(rep_seq), processRepeatMasker, mc.cores = threads)
#repeat_class = mclapply(1:100, processRepeatMasker, mc.cores = threads)

job_err = repeat_class[1:nrow(rep_seq)]
stopifnot(length(which(grepl("Error", job_err)))==0)

comp_repeat = cbind.data.frame(rep_seq$component, rep_seq$INS_id, rep_seq$cluster_id, unlist(repeat_class), stringsAsFactors = F)
colnames(comp_repeat) = c("component", "INS_id", "cluster_id", "repeat_class")

write.table(comp_repeat, paste0(dir, "/discovery/comp_repeat.txt"), col.names = T, row.names = F, quote = F, sep = '\t')


