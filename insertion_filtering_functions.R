readBN = function(sample, BN_path){
  
  # Validate the BN SV if larger than 1kb
  BN = read.table(paste0(BN_path , "/", BN_list[grepl(sample, BN_list)]), stringsAsFactors = F)
  BN = BN[,c(3,7,8,10,18,25,5,6,26)]
  colnames(BN) = c("ref_chr", "ref_start", "ref_end", "type", "zygosity", "size", "map_start", "map_end", "enzyme")
  BN = BN[BN$type=="insertion",] # Only subset insertions
  BN$label_dist = BN$map_end - BN$map_start
  
  # Change the chr name from 1-24 to chr1-chrY
  BN$ref_chr = paste0("chr",BN$ref_chr)
  BN$ref_chr = gsub("chr23", "chrX", BN$ref_chr)
  BN$ref_chr = gsub("chr24", "chrY", BN$ref_chr)
  
  # BN smap has imprecise breakpoints when the SVs are large, generally when the SVs involve more than 1 label
  # Thus, for all insertions, add 5k to each end to make sure they overlap with the insertion data
  BN$ref_start = BN$ref_start - 5000
  BN$ref_end = BN$ref_end + 5000
  
  return(BN)
  
}


removeSegdupSVblacklist = function(assemblytics, assemblytics.gr, segdup.gr, sv_bl.gr){
  
  #Check segdup and remove if overlapping
  segdup_ov = findOverlaps(assemblytics.gr, segdup.gr, type = "any")
  sv_bl_ov = findOverlaps(assemblytics.gr, sv_bl.gr, type = "any")
  # Discard entries if overlapping segdup
  discard = unique(c(segdup_ov@from, sv_bl_ov@from))
  
  if (length(discard)!=0){
    assemblytics = assemblytics[-discard,]
  }
  
  return(assemblytics)
}


overlapBN = function(BN_sample, assemblytics, assemblytics.gr){
  assemblytics$BN_size = 0
  assemblytics$label_dist = 0
  assemblytics$BN_start = NA # This start corresponds to the reference position
  assemblytics$BN_end = NA # This end corresponds to the reference position
  assemblytics$BN_enzyme = NA
  
  if (!is.null(BN_sample)){
    
    # Find overlap between the assemblytics and BN
    # Need to recreate assemblytics.gr
    BN_sample.gr = makeGRangesFromDataFrame(BN_sample, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
    BN_ov = findOverlaps(assemblytics.gr, BN_sample.gr, type="any")
    
    diff_DF = NULL
    for (k in unique(BN_ov@from)){
      diff = NULL
      df = BN_ov[BN_ov@from==k]
      assemblytics_size = assemblytics$insert_size[k]
      diff$value = abs(as.numeric(BN_sample$size[df@to])-as.numeric(assemblytics_size))
      diff$enzyme = BN_sample$enzyme[df@to]
      diff$dist = BN_sample$label_dist[df@to]
      
      diff_df = as.data.frame(diff)
      diff_DF = rbind.data.frame(diff_DF, diff_df)
      
      # If same enzyme, pick the one with min value
      if (any((table(diff$enzyme))>1)){
        doublet = diff$enzyme[(table(diff$enzyme))>1]
        discard = which.min(diff$value[diff$enzyme==doublet])
        discard_pos = which(diff$enzyme==doublet)[discard]
        diff$dist[discard_pos] = NA
      }
      
      # If diff enzyme, pick the one with min label_dist
      BN_keep = which.min(diff$dist)
      assemblytics$BN_size[k] = BN_sample$size[df@to[BN_keep]]
      assemblytics$label_dist[k] = BN_sample$label_dist[df@to[BN_keep]]
      assemblytics$BN_start[k] = BN_sample$ref_start[df@to[BN_keep]]
      assemblytics$BN_end[k] = BN_sample$ref_end[df@to[BN_keep]]
      assemblytics$BN_end[k] = BN_sample$ref_end[df@to[BN_keep]]
      assemblytics$BN_enzyme[k] = diff$enzyme[BN_keep]
    }
    
  }
  return(assemblytics)
}


overlapNgap = function(assemblytics, ngapFile) {

  # prepare assemblytics columns for overlapping
  assemblytics$adjusted_assm_start = as.numeric(assemblytics$adjusted_assm_start)
  assemblytics$adjusted_assm_end = as.numeric(assemblytics$adjusted_assm_end)
  assemblytics$boundary_left_start = as.numeric(assemblytics$adjusted_assm_start) - 10
  assemblytics$boundary_left_end = as.numeric(assemblytics$adjusted_assm_start) + 10
  assemblytics$boundary_right_start = as.numeric(assemblytics$adjusted_assm_end) - 10
  assemblytics$boundary_right_end = as.numeric(assemblytics$adjusted_assm_end) + 10
  
  if (is.null(ngapFile)){
    
    assemblytics$ngap = 0
    assemblytics$ngap_boundaries = 0
    assemblytics$ngap_boundaries_size_left = 0
    assemblytics$ngap_boundaries_size_right = 0
    
  } else {
    
    colnames(ngapFile) = c("ref_assm", "start", "end", "ngap_size")
    ngap.gr = makeGRangesFromDataFrame(ngapFile, seqnames.field = "ref_assm", start.field = "start", end.field = "end")

    assemblytics_ngap.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "adjusted_assm_start", end.field = "adjusted_assm_end")
    assemblytics_boundary_left.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "boundary_left_start", end.field = "boundary_left_end")
    assemblytics_boundary_right.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "boundary_right_start", end.field = "boundary_right_end")
    
    # Actual overlap of Ngap
    assemblytics_ngap_ngap_ov = findOverlaps(assemblytics_ngap.gr, ngap.gr, type="any")
    assemblytics_boundary_left_ngap_ov = findOverlaps(assemblytics_boundary_left.gr, ngap.gr, type="any")
    assemblytics_boundary_right_ngap_ov = findOverlaps(assemblytics_boundary_right.gr, ngap.gr, type="any")
    
    assemblytics$ngap = 0
    # Report the combined Ngap size in that overlapping range
    for (n in unique(assemblytics_ngap_ngap_ov@from)){
      
      df = assemblytics_ngap_ngap_ov[assemblytics_ngap_ngap_ov@from==n]
      #assemblytics_size = assemblytics$insert_size[k]
      ngap_combined_size = sum(as.numeric(ngapFile$ngap_size[df@to]))
      assemblytics$ngap[n] = ngap_combined_size
      
    }
    
    assemblytics$ngap_boundaries = "no"
    assemblytics$ngap_boundaries[assemblytics_boundary_left_ngap_ov@from] = "yes"
    assemblytics$ngap_boundaries[assemblytics_boundary_right_ngap_ov@from] = "yes"
    
    # Need to get the size of Ngap at boundaries
    # This is different from the total Ngap size
    assemblytics$ngap_boundaries_size_left = 0
    assemblytics$ngap_boundaries_size_right = 0
    assemblytics$ngap_boundaries_size_left[assemblytics_boundary_left_ngap_ov@from] = ngapFile$ngap_size[assemblytics_boundary_left_ngap_ov@to]
    assemblytics$ngap_boundaries_size_right[assemblytics_boundary_right_ngap_ov@from] = ngapFile$ngap_size[assemblytics_boundary_right_ngap_ov@to]
    
    # Only keep if insert_size is larger than ngap
    assemblytics = assemblytics[assemblytics$insert_size > assemblytics$ngap, ]
    
  }
  
  assemblytics$ngap_perct = assemblytics$ngap/assemblytics$adjusted_insert_size
  
  return(assemblytics)
}


# Remove redundant ref coordinates
removeRedundantInsertion = function(assemblytics){
  
  dup_start = names(which(table(assemblytics$ref_start) > 1))
  dup_end = names(which(table(assemblytics$ref_end) > 1))
  
  assemblytics = assemblytics[(!assemblytics$ref_start %in% dup_start) & (!assemblytics$ref_end %in% dup_end),]
  #df = assemblytics[(assemblytics$ref_start %in% dup_start) | (assemblytics$ref_end %in% dup_end),]
  #keep = which(!duplicated(assemblytics[,c("ref_chr", "ref_start", "ref_end")]))
  #assemblytics = assemblytics[keep,]
  
  return(assemblytics)
}

processAssmCoords = function(assemblytics){
    
  # Split the assm_coords into 3 columns
  assemblytics$assm_id = str_split_fixed(assemblytics$assm_coords, ":|-", 4)[,1]
  assemblytics$assm_start = as.numeric(str_split_fixed(assemblytics$assm_coords, ":|-", 4)[,2]) 
  assemblytics$assm_end = as.numeric(str_split_fixed(assemblytics$assm_coords, ":|-", 4)[,3]) 
  
  assemblytics$tmp = assemblytics$assm_start
  
  assemblytics$assm_start = ifelse(assemblytics$tmp > assemblytics$assm_end, assemblytics$assm_end, assemblytics$assm_start)
  assemblytics$assm_end = ifelse(assemblytics$tmp > assemblytics$assm_end, assemblytics$tmp, assemblytics$assm_end)
  
  assemblytics = assemblytics[, -which(names(assemblytics) %in% "tmp")]
  
}

processAdjustedCoords = function(assemblytics){
  
    # Parse the adjusted_coordinates column
    # **1. Remove INS if adjusted_coordinates can't be computed given ref_gap_size is negative
    discard = grep("\\[", assemblytics$adjusted_coords)
    if (length(discard)>0){
      assemblytics = assemblytics[-discard,]
    }
    assemblytics$adjusted_assm_start = NA
    assemblytics$adjusted_assm_end = NA
    
    # **2. Split assemblytics into 2 dataframes (negative ref_gap_size and positive ref_gap_size)
    # we will only use the adjusted_coordinates column for neg_ref_gap
    neg_ref_gap = assemblytics[assemblytics$ref_gap_size<0,]
    pos_ref_gap = assemblytics[assemblytics$ref_gap_size>=0,]
    
    # **2. Make sure only two sets of coordinates are reported
    # in neg_ref_gap, we need to make sure assm_start and assm_end are found in the adjusted_coordinates column
    neg_ref_gap_list = sapply(strsplit(neg_ref_gap$adjusted_coords, ";"), unique)
    discard = NULL
    for (s in 1:length(neg_ref_gap_list)){
      
      if (length(neg_ref_gap_list[[s]]!=2)){
          discard = c(discard, s)
          next
      }
      
      startPattern = neg_ref_gap$assm_start[s]
      endPattern = neg_ref_gap$assm_end[s]
      
      adjusted_assm_start = str_split_fixed(neg_ref_gap_list[[s]][grep(startPattern,neg_ref_gap_list[[s]])], "-", 2)
      if (nrow(adjusted_assm_start)==1) {
        neg_ref_gap$adjusted_assm_start[s] = adjusted_assm_start[which(adjusted_assm_start != startPattern)]
      } else if (nrow(adjusted_assm_start)==2) {
        neg_ref_gap$adjusted_assm_start[s] = min(adjusted_assm_start[which(!adjusted_assm_start %in% startPattern)])
        neg_ref_gap$adjusted_assm_end[s] = max(adjusted_assm_start[which(!adjusted_assm_start %in% startPattern)])
        next
      }
      else discard = c(discard, s)
      
      adjusted_assm_end = str_split_fixed(neg_ref_gap_list[[s]][grep(endPattern,neg_ref_gap_list[[s]])], "-", 2)
      if (nrow(adjusted_assm_end)==1) neg_ref_gap$adjusted_assm_end[s] = adjusted_assm_end[which(adjusted_assm_end != endPattern)]
      else discard = c(discard, s)
    }
    neg_ref_gap = neg_ref_gap[-discard,]
    
    neg_ref_gap$tmp_s = as.numeric(neg_ref_gap$adjusted_assm_start)
    neg_ref_gap$tmp_e = as.numeric(neg_ref_gap$adjusted_assm_end)
    
    neg_ref_gap$adjusted_assm_end[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)] = neg_ref_gap$tmp_s[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)]
    neg_ref_gap$adjusted_assm_start[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)] = neg_ref_gap$tmp_e[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)]
    neg_ref_gap = subset(neg_ref_gap, select = -c(tmp_s, tmp_e))
    
    # **3. If pos_ref_gap, adjusted_assm_start and adjusted_assm_end are identical to assm_start and assm_end
    pos_ref_gap$adjusted_assm_start = as.numeric(pos_ref_gap$assm_start)
    pos_ref_gap$adjusted_assm_end = as.numeric(pos_ref_gap$assm_end)
    
    # **4. Put the two lists together
    assemblytics = rbind.data.frame(pos_ref_gap, neg_ref_gap, stringsAsFactors = F)
    
    # **5. Add a column called adjusted insert size
    assemblytics$adjusted_insert_size = as.numeric(assemblytics$adjusted_assm_end) - as.numeric(assemblytics$adjusted_assm_start)
    
    return(assemblytics)
}

checkAssmCoordsWithinLen = function(assemblytics, faidx, sample, haplo){
  
  # subset the index file
  keep = which(faidx$sample==sample & faidx$haplo==haplo)
  fai = faidx[keep,]
  
  assemblytics = merge(assemblytics,fai[,c("assm_id", "scaffold_length")], by = "assm_id")
  assemblytics = assemblytics[(as.numeric(assemblytics$adjusted_assm_start)>0 & as.numeric(assemblytics$adjusted_assm_end)<assemblytics$scaffold_length),]
  
  return(assemblytics)
}

assignComponent=function(i)  {
  
  # break dataframe by chr
  assemblytics_combined_per_chr = assemblytics_combined[assemblytics_combined$ref_chr==i,]
  
  # find overlapping reference ranges
  assemblytics.gr = makeGRangesFromDataFrame(assemblytics_combined_per_chr, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
  assemblytics_assemblytics_ov_self = findOverlaps(assemblytics.gr, assemblytics.gr, type = "any", ignore.strand=T)
  
  edges = as.data.frame(assemblytics_assemblytics_ov_self)
  net = graph_from_data_frame(d=edges, directed = F)
  membership = components(net)$membership
  components_count = count_components(net)
  assemblytics_combined_per_chr$component = membership # add the component number to each INS
  
  return(assemblytics_combined_per_chr)
}


BN_validate = function(df){
  
  df$BN_validated = NA
  df$lower_bound = NA
  df$upper_bound = NA
  # BN validation <== this requires a set of rules
  # 1. if SVs are smaller than 1kb, size difference has to be <=200bp to be validated
  # 2. if SVs are bigger than 1kb, size difference has to be <=700bp to be considered validated
  # 3. If label_dist is greater than 40kb and if the SV sizes are discordant, it will not get penalized or validated
  
  if (all(is.na(df$BN_size))){
    
    df$BN_size = 0
    df$label_dist = 0
    
  } else {
    
    df$lower_bound[df$insert_size<1000] = df$insert_size[df$insert_size<1000] - 200
    df$upper_bound[df$insert_size<1000] = df$insert_size[df$insert_size<1000] + 200
    df$lower_bound[df$insert_size>=1000] = df$insert_size[df$insert_size>=1000] - 700
    df$upper_bound[df$insert_size>=1000] = df$insert_size[df$insert_size>=1000] + 700
    
    df$BN_validated = ifelse(df$BN_size>=df$lower_bound & df$BN_size<=df$upper_bound, 1, 0)
    df$BN_validated[df$BN_size==0] = -1
    
    # ADD ANOTHER SCENARIO WHEN MORE THAN 1 INSERTIONS ARE SPANNED BY 1 BN SV CALL
    ## BN SIZE IS THE OVERAL SIZE CHANGE
    df_BN.gr = makeGRangesFromDataFrame(df[!is.na(df$BN_start),], seqnames.field = "ref_chr", start.field = "BN_start", end.field = "BN_end")
    df_ref.gr = makeGRangesFromDataFrame(df, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
    df_BN_ref_ov = findOverlaps(df_ref.gr, df_BN.gr, type = "within")
    
    keep = which(duplicated(df_BN_ref_ov@to))
    for (k in keep){
      s_value = df_BN_ref_ov@to[k]
      rows = df_BN_ref_ov[which(df_BN_ref_ov@to==s_value)]
      sub_df = df[rows@from,]
      total_insert_size = sum(sub_df$insert_size)
      diff = abs(total_insert_size - sub_df$BN_size[1])
      
      if (all(sub_df$BN_validated=="no") & diff<=500){
        df[rows@from,"BN_validated"]="yes"
      }
    }
    # If insert_size is smaller than 500bp, BN validation is not used
    df$BN_validated[df$insert_size<=500] = -1
    
    # Finally, check the label_dist
    #df$BN_validated[df$label_dist>=40000 & df$BN_validated=="no"] = -1
    
  }
  return(df)
  
}


cal_edge_size = function(l, standard_inc=0.5){
  lenl = length(l)
  if (lenl<2) edge = NA
  else {
    l_sort = sort(l)
    edge = sum(log10(l_sort[2:lenl]-l_sort[1:(lenl-1)]+1))/((lenl-1)*log10(standard_inc+1))
    edge = 10^edge - 1
  }
  return(edge)
}

findRepresentativeSeq = function(assm_sub_df, res){
  
  rep_seq = NULL
  target_id = str_split_fixed(res$group[i], ";|,", 2)[,1]

  # retrieve line
  rep_seq[[1]] = assm_sub_df[assm_sub_df$INS_id==target_id,]
  rep_seq[[1]]$type = 0 # type 0 means it'll stay in the main analysis pipeline; 1 means it'll not be included in the main analysis (supplementary dataset)
  
  # identify all clusters
  all_clus = strsplit(res$group[i], ",")
  
  # extract the main cluster
  main_clus_id = unlist(strsplit(all_clus[[1]][1], ";"))
  main_clus_df = assm_sub_df[assm_sub_df$INS_id %in% main_clus_id,]
  
  # calculate edge score for main cluster if unique(sample)>1
  if (length(unique(main_clus_df$sample))>1){
    
      main_edge_df = main_clus_df[main_clus_df$ngap_boundaries_size_left==0 & main_clus_df$ngap_boundaries_size_right==0,]
      rep_seq[[1]]$edge_start = cal_edge_size(main_edge_df$ref_start)
      rep_seq[[1]]$edge_end = cal_edge_size(main_edge_df$ref_end)
      
  } else{
    
      rep_seq[[1]]$edge_start = NA
      rep_seq[[1]]$edge_end = NA
      
  }
  
  clus_start = min(main_clus_df$ref_start)
  clus_end = max(main_clus_df$ref_end)
  
  rep_seq[[1]]$cluster_id = paste0(rep_seq[[1]]$component, "_1" )
  
  rep_seq[[1]] = c(rep_seq[[1]], as.list(annotateRepSeq(main_clus_df))) 
  
  # If there are more than one clusters
  if (res$cluster[i]>1){
    
    for (l in 2:res$cluster[i]){ # for all other clusters...

        cluster_id = unlist(strsplit(all_clus[[1]][l], ";"))
        clus_df = assm_sub_df[assm_sub_df$INS_id %in% cluster_id,]
  
        rep_seq[[l]] = assm_sub_df[assm_sub_df$INS_id == cluster_id[1],]
        
        curr_start = min(clus_df$ref_start)
        curr_end = max(clus_df$ref_end)
        
        # check to see if other clusters overlap with existing cluster based on ref coordinates
        rep_seq[[l]]$type = ifelse(all(clus_end <  curr_start) | all(clus_start > curr_end), 0, 1)

        clus_start = c(clus_start, curr_start)
        clus_end = c(clus_end, curr_end)

        # for these clusters, calculate edge scores
        if (length(unique(clus_df$sample))>1){
          
          edge_df = clus_df[clus_df$ngap_boundaries_size_left==0 & clus_df$ngap_boundaries_size_right==0,]
          rep_seq[[l]]$edge_start = cal_edge_size(edge_df$ref_start)
          rep_seq[[l]]$edge_end = cal_edge_size(edge_df$ref_end)
          
        } else{
          
          rep_seq[[l]]$edge_start = NA
          rep_seq[[l]]$edge_end = NA
          
        }
        
        rep_seq[[l]]$cluster_id = paste0(rep_seq[[l]]$component, "_", l)
        
        rep_seq[[l]] = c(rep_seq[[l]], as.list(annotateRepSeq(clus_df))) 
        
      }

  }

  rep_seq = ldply(rep_seq, data.frame)
  return(rep_seq)
  
}

annotateRepSeq = function(clus_df){
  
    # add ethnicity info to rep_seq
    keep = which(metadata$sample %in% unique(clus_df$sample))
    sample_ethnicity = (metadata[keep,"population"])
    
    total_sample_count = nrow(metadata)
    sample_count = length(keep)
    sample_perct = sample_count/total_sample_count
    sample_record = unlist(clus_df$INS_id)
    sample_record = paste(unique(sample_record), collapse = ';')
    ethnicity_tbl = c(table(sample_ethnicity)[c("AFR", "AMR", "EAS", "EUR", "SAS")])
    ethnicity_tbl = ifelse(is.na(ethnicity_tbl), 0, ethnicity_tbl)
    
    metadata_ethnicity_tbl = c(table(metadata$population)[c("AFR", "AMR", "EAS", "EUR", "SAS")])
    metadata_ethnicity_tbl = ifelse(is.na(metadata_ethnicity_tbl), 0, metadata_ethnicity_tbl)
    
    ethnicity_prct = ethnicity_tbl/metadata_ethnicity_tbl
    
    sample_nonAFR = sum(ethnicity_tbl[c(2:5)])
    percent_nonAFR = sample_nonAFR/nrow(metadata[metadata$population!="AFR",])
    
    anno = c(sample_count, sample_perct, sample_record, as.vector(ethnicity_tbl), sample_nonAFR, as.vector(ethnicity_prct), percent_nonAFR)
    names(anno) = c("sample_count", "sample_perct", "sample_record", "sample_AFR", "sample_AMR", "sample_EAS", "sample_EUR", "sample_SAS", "sample_nonAFR", 'percent_AFR', "percent_AMR", 'percent_EAS', "percent_EUR", "percent_SAS", "percent_nonAFR")
    
    return(anno)

}













calc_min_dist_to_exon = function(i){
  
  print(i)
  breakpoint1 = rep_seq$ref_start[i]
  breakpoint2 = rep_seq$ref_end[i]
  target_chr = rep_seq$ref_chr[i]
  # Subset the exon file
  exon_subset = exon_df[which(exon_df$exon_chr==target_chr),]
  
  # Compute the minimum distance from either breakpoint to the nearest exon
  min_dist_to_exon = min(abs(exon_subset$exon_start - breakpoint1), abs(exon_subset$exon_end - breakpoint1), abs(exon_subset$exon_start - breakpoint2), abs(exon_subset$exon_end - breakpoint2))
  
  return(min_dist_to_exon)
}

checkResAssemblyticsComponentConcord = function(assemblytics, res){
  
  discard = which(! assemblytics$component %in% res$component)
  if (length(discard)!=0){
    
    assemblytics = assemblytics[-discard,]
    
  }
  
  # sort assemblytics by component
  assemblytics = assemblytics[order(assemblytics$component),]
  
  # if (any(unique(assemblytics$component) != unique(res$component))){
  #     stop("Component not identical between assemblytics and multiple alignment results")
  # }
  
  return(assemblytics)

}

# 
# findRepresentativeSeq = function(i){
#   
#     print(i)
#     #print("haha")
#     
#     comp = comp_list[i]
#     df = assemblytics_combined_per_chr[(assemblytics_combined_per_chr$component==comp),]
#     
#     # Duplicate the lines of haplo=="unphased"
#     dup = df[df$haplo=="unphased", ]
#     df = rbind.data.frame(df, dup, stringsAsFactors = F)
#     
#     # Initialize the edge variables
#     edge_start = NA
#     edge_end = NA
#     
#     # testing
#     df_clean = (df[df$ngap==0,c("ref_start", "ref_end")])
#     if (nrow(df_clean) > 2){
#       
#         edge_start = cal_edge_size(df_clean$ref_start)
#         edge_end = cal_edge_size(df_clean$ref_end)
#         #edge_df = rbind.data.frame(edge_df, c(component, edge_start, edge_end))
#       
#     }
#     
#     # Add ethnicity info to df and analyze populations
#     df = merge(df,metadata[,c("sample", "population")], by = "sample", all.x = T)
#     sample_df = unique(df[,c("sample", "population")])
#     sample_count = nrow(sample_df)
#     
#     for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS")){
#       assign(paste0("sample_", pop), nrow(sample_df[sample_df$population==pop,]))
#       assign(paste0("percent_", pop), round(eval(parse(text=paste0("sample_",pop)))/length(which(metadata$population==pop)),3))
#     }
#     
#     sample_record = unlist(paste0(df$INS_id,"_",df$population))
#     sample_record = paste(unique(sample_record), collapse = ';')
#     
#     sample_nonAFR = sum(sample_AMR, sample_EAS, sample_EUR, sample_SAS)
#     percent_nonAFR = sample_nonAFR/nrow(metadata[metadata$population!="AFR",])
#       
#     # Always use sequences without ngap
#     if (length(which(df$ngap==0)) >= 1){
#         df = df[df$ngap==0,]
#     } else { # If sequences have ngaps, sort based on shortest ngap
#         df = df[order(df$BN_validated, decreasing = T),]
#         df = df[order(df$ngap, df$ngap_boundaries),]
#     }
#     
#     # Choose the most represented breakpoints, then insert size
#     start_end_size_df = as.data.frame(table(df[,c("ref_start", "ref_end", "insert_size")]))
#     rep_bkpt_size = start_end_size_df[which.max(start_end_size_df$Freq),c("ref_start", "ref_end", "insert_size")]
#     
#     df = df[df$ref_start==rep_bkpt_size$ref_start & df$ref_end==rep_bkpt_size$ref_end & df$insert_size==rep_bkpt_size$insert_size, ]
#     
#     newline = df[1,c(2:23,1,24:32)]
#     newline = cbind.data.frame(newline, sample_record, sample_count, sample_AFR, sample_AMR, sample_EAS, sample_EUR, sample_SAS, sample_nonAFR, percent_AFR, percent_AMR, percent_EAS, percent_EUR, percent_SAS, percent_nonAFR, edge_start, edge_end)
#     #print(newline)
#     #representative_seq = rbind.data.frame(representative_seq, newline, stringsAsFactors = F)
#     representative_seq[[i]] <<- newline
#     #return(representative_seq)
#     #newline
# }
# 
