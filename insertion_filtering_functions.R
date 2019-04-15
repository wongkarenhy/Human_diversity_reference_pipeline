readBN = function(sample, BN_path){
  
  # Validate the BN SV if larger than 1kb
  BN = read.table(paste0(BN_path, sample,".smap"), stringsAsFactors = F)
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

removeSegdup = function(lastz_bkpt, lastz_bkpt.gr, segdup.gr, sv_bl.gr){
  
  #Check segdup and remove if overlapping
  segdup_ov = findOverlaps(lastz_bkpt.gr, segdup.gr, type = "any")
  sv_bl_ov = findOverlaps(lastz_bkpt.gr, sv_bl.gr, type = "any")
  # Discard entries if overlapping segdup
  discard = unique(c(segdup_ov@from, sv_bl_ov@from))
  
  if (length(discard)!=0){
    lastz_bkpt = lastz_bkpt[-discard,]
  }
  
  return(lastz_bkpt)
}

overlapBN = function(BN_sample, lastz_bkpt, lastz_bkpt.gr){
  lastz_bkpt$BN_size = 0
  lastz_bkpt$label_dist = 0
  lastz_bkpt$BN_start = NA # This start corresponds to the reference position
  lastz_bkpt$BN_end = NA # This end corresponds to the reference position
  lastz_bkpt$BN_enzyme = NA
  
  if (!is.null(BN_sample)){
    
    # Find overlap between the lastz_bkpt and BN
    # Need to recreate lastz_bkpt.gr
    BN_sample.gr = makeGRangesFromDataFrame(BN_sample, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
    BN_ov = findOverlaps(lastz_bkpt.gr, BN_sample.gr, type="any")
    
    diff_DF = NULL
    for (k in unique(BN_ov@from)){
      diff = NULL
      df = BN_ov[BN_ov@from==k]
      lastz_size = lastz_bkpt$insert_size[k]
      diff$value = abs(as.numeric(BN_sample$size[df@to])-as.numeric(lastz_size))
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
      lastz_bkpt$BN_size[k] = BN_sample$size[df@to[BN_keep]]
      lastz_bkpt$label_dist[k] = BN_sample$label_dist[df@to[BN_keep]]
      lastz_bkpt$BN_start[k] = BN_sample$ref_start[df@to[BN_keep]]
      lastz_bkpt$BN_end[k] = BN_sample$ref_end[df@to[BN_keep]]
      lastz_bkpt$BN_end[k] = BN_sample$ref_end[df@to[BN_keep]]
      lastz_bkpt$BN_enzyme[k] = diff$enzyme[BN_keep]
    }
    
  }
  return(lastz_bkpt)
}

overlapNgap = function(lastz_bkpt, ngap, lastz_ngap.gr, ngap.gr, lastz_boundary_left.gr, lastz_boundary_right.gr) {
  
  lastz_ngap_ngap_ov = findOverlaps(lastz_ngap.gr, ngap.gr, type="any")
  lastz_boundary_left_ngap_ov = findOverlaps(lastz_boundary_left.gr, ngap.gr, type="any")
  lastz_boundary_right_ngap_ov = findOverlaps(lastz_boundary_right.gr, ngap.gr, type="any")
  
  lastz_bkpt$ngap = 0
  # Report the combined Ngap size in that overlapping range
  for (n in unique(lastz_ngap_ngap_ov@from)){
    
    df = lastz_ngap_ngap_ov[lastz_ngap_ngap_ov@from==n]
    #lastz_size = lastz_bkpt$insert_size[k]
    ngap_combined_size = sum(as.numeric(ngap$ngap_size[df@to]))
    lastz_bkpt$ngap[n] = ngap_combined_size
    
  }
  
  lastz_bkpt$ngap_boundaries = "no"
  lastz_bkpt$ngap_boundaries[lastz_boundary_left_ngap_ov@from] = "yes"
  lastz_bkpt$ngap_boundaries[lastz_boundary_right_ngap_ov@from] = "yes"
  
  # Need to get the size of Ngap at boundaries
  # This is different from the total Ngap size
  lastz_bkpt$ngap_boundaries_size_left = 0
  lastz_bkpt$ngap_boundaries_size_right = 0
  lastz_bkpt$ngap_boundaries_size_left[lastz_boundary_left_ngap_ov@from] = ngap$ngap_size[lastz_boundary_left_ngap_ov@to]
  lastz_bkpt$ngap_boundaries_size_right[lastz_boundary_right_ngap_ov@from] = ngap$ngap_size[lastz_boundary_right_ngap_ov@to]
  
  # Only keep if insert_size is larger than ngap
  lastz_bkpt = lastz_bkpt[lastz_bkpt$insert_size > lastz_bkpt$ngap, ]
  
  return(lastz_bkpt)
}

# Remove redundant ref coordinates
removeRedundantInsertion = function(lastz_bkpt){
  keep = which(!duplicated(lastz_bkpt[,c("ref_chr", "ref_start", "ref_end")]))
  lastz_bkpt = lastz_bkpt[keep,]
  
  return(lastz_bkpt)
}

# Combine the two pseudohaps
combinePseudohaps = function(sample_list, lastz_bkpt_ALL){
  
    lastz_bkpt_ALL_pseudohap_combined = NULL
    for (i in sample_list){ # Loop through every sample
      
      print(i)
      
      df = lastz_bkpt_ALL[lastz_bkpt_ALL$sample==i,]
      df = df[order(df$ref_chr, df$ref_start),]
      
      dup_list = which(duplicated(df[,c("ref_chr", "ref_start", "ref_end")]))
      df$haplo[dup_list-1] = paste0(df$haplo[dup_list-1],";",df$haplo[dup_list])
      df = df[-dup_list,]
      
      df$validated = NA
      df$lower_bound = NA
      df$upper_bound = NA
      # BN validation <== this requires a set of rules
      # 1. if SVs are smaller than 1kb, size difference has to be <=200bp to be validated
      # 2. if SVs are bigger than 1kb, size difference has to be <=500bp to be considered validated
      # 3. If label_dist is greater than 40kb and if the SV sizes are discordant, it will not get penalized or validated
      
      if (all(is.na(df$BN_size))){
        df$BN_size = 0
        df$label_dist = 0
        # directly append to dataframe
        lastz_bkpt_ALL_pseudohap_combined = rbind.data.frame(lastz_bkpt_ALL_pseudohap_combined,df, stringsAsFactors = F)
        next
      }
      
      df$lower_bound[df$insert_size<1000] = df$insert_size[df$insert_size<1000] - 200
      df$upper_bound[df$insert_size<1000] = df$insert_size[df$insert_size<1000] + 200
      df$lower_bound[df$insert_size>=1000] = df$insert_size[df$insert_size>=1000] - 700
      df$upper_bound[df$insert_size>=1000] = df$insert_size[df$insert_size>=1000] + 700
      
      df$validated = ifelse(df$BN_size>=df$lower_bound & df$BN_size<=df$upper_bound, "yes", "no")
      df$validated[df$BN_size==0] = NA
      
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
        
        if (all(sub_df$validated=="no") & diff<=500){
          df[rows@from,"validated"]="yes"
        }
      }
      # If insert_size is smaller than 500bp, BN validation is not used
      df$validated[df$insert_size<=500] = NA
      
      # Finally, check the label_dist
      df$validated[df$label_dist>=40000 & df$validated=="no"] = NA
      
      lastz_bkpt_ALL_pseudohap_combined = rbind.data.frame(lastz_bkpt_ALL_pseudohap_combined,df, stringsAsFactors = F)
      
    }
    return(lastz_bkpt_ALL_pseudohap_combined)
}

findRepresentativeSeq = function(membership, components_count, lastz_combined, metadata){
  
  # Group the alns based on component membership and find the most representative alignment
  representative_seq = NULL
  for (i in 1:components_count){
    
    component = i
    member = names(which(membership==component))
    df = lastz_combined[as.numeric(member),]
    
    # Add ethnicity info to df and analyze populations
    df = merge(df,metadata[,c("sample", "population")], by = "sample", all.x = T)
    sample_df = unique(df[,c("sample", "population")])
    sample_count = nrow(sample_df)
    
    for (pop in unique(metadata$population)){
      assign(paste0("sample_", pop), nrow(sample_df[sample_df$population==pop,]))
      assign(paste0("percent_", pop), round(eval(parse(text=paste0("sample_",pop)))/length(which(metadata$population==pop)),3))
    }
    
    sample_record = unlist(paste0(sample_df$sample,"_",sample_df$population))
    sample_record = paste(sample_record, collapse = ';')
    
    sample_nonAFR = sum(sample_AMR, sample_EAS, sample_EUR, sample_SAS)
    percent_nonAFR = sample_nonAFR/nrow(metadata[metadata$population!="AFR",])
      
    # Always use sequences without ngap
    if (length(which(df$ngap==0)) >= 1){
        df = df[df$ngap==0,]
    } else { # If sequences have ngaps, sort based on shortest ngap
        df = df[order(df$validated, decreasing = T),]
        df = df[order(df$ngap, df$ngap_boundaries),]
    }
    
    # Choose the most represented breakpoints, then insert size
    start_end_size_df = as.data.frame(table(df[,c("ref_start", "ref_end", "insert_size")]))
    rep_bkpt_size = start_end_size_df[which.max(start_end_size_df$Freq),c("ref_start", "ref_end", "insert_size")]
    
    df = df[df$ref_start==rep_bkpt_size$ref_start & df$ref_end==rep_bkpt_size$ref_end & df$insert_size==rep_bkpt_size$insert_size, ]
    
    newline = df[1,c(2:22,1,23:27)]
    newline = cbind.data.frame(newline, component, sample_record, sample_count, sample_AFR, sample_AMR, sample_EAS, sample_EUR, sample_SAS, sample_nonAFR, percent_AFR, percent_AMR, percent_EAS, percent_EUR, percent_SAS, percent_nonAFR)
    representative_seq = rbind.data.frame(representative_seq, newline, stringsAsFactors = F)
    
  }
  return(representative_seq)
}


