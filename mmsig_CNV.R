round_to_half = function(x) {
  ifelse(abs(x - 0.5) <= 0.4, 0.5, 
         ifelse(abs(x - 0.5) <= 0.6, round(x), 
                round(x * 2) / 2))
}


collapse_segments = function(df) {
  # Sort the data frame by chromosome and start position
  df=df[order(df$chr, df$start), ] 
  
  # Initialize an empty list to store the collapsed segments
  collapsed=list()
  
  # Iterate through each chromosome
  for (chr in unique(df$chr)) {
    chr_df=df[df$chr == chr, ]
    
    # Initialize a list to store segments for the current chromosome
    chr_segments=list()
    
    # Iterate through each segment in the chromosome
    for (i in 1:nrow(chr_df)) {
      if (length(chr_segments) == 0) {
        # If the list is empty, add the first segment
        chr_segments[[1]]=chr_df[i, ]
      } else {
        # Check if the current segment is contiguous with the last segment in the list
        if (chr_segments[[length(chr_segments)]]$end + 1 >= chr_df$start[i]) {
          # If contiguous, update the end position of the last segment
          chr_segments[[length(chr_segments)]]$end=chr_df$end[i]
        } else {
          # If not contiguous, add a new segment to the list
          chr_segments[[length(chr_segments) + 1]]=chr_df[i, ]
        }
      }
    }
    
    # Convert the list of segments to a data frame
    chr_collapsed <- do.call(rbind, chr_segments)
    
    # Add the collapsed segments for the current chromosome to the final list
    collapsed <- c(collapsed, list(chr_collapsed))
  }
  
  # Combine the collapsed segments for all chromosomes into a single data frame
  collapsed_df <- do.call(rbind, collapsed)
  
  return(collapsed_df)
}


gg_color_hue=function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


generate_cnv_feature_matrix <- function(seg_file) {
  # Remove chrX, NA major, and sort
  df = seg_file
  df$Chrom = gsub("chr", "", df$Chrom)
  df = df[df$Chrom != "X", ]
  df = df[!is.na(df$major), ]
  df = df[order(df$sample, df$Chrom, df$start), ]
  
  # Remove IG loci
  IG_loci = (
    (df$Chrom == 14 & df$start >106032614 &  df$end< 108288051) |
      (df$Chrom == 22 & df$start >21080474 &  df$end< 26065085) |
      (df$Chrom == 2 & df$start >87090568 &  df$end< 93274235)
  )
  df = df[!IG_loci, ]
  
  # Remove segments < 50kb
  df = df[(df$end - df$start) > 50000, ]
  
  # Remove Battenberg artifacts
  bat_art = df[(df$end-df$start) < 1e6 & df$minor == 0 & df$major > 2, ]
  df = df[!paste(df$sample, df$Chrom, df$start) %in% paste(bat_art$sample, bat_art$Chrom, bat_art$start), ]
  
  # Collapse adjacent segments with identical copy number
  df = collapse_segments(df)
  df$start = as.numeric(as.character(df$start))
  df$end = as.numeric(as.character(df$end))
  
  # --- Define feature bins (Maclachlan et al. 2021) ---
  count_cnv_coordinate = data.frame(code=1:5, count=c(0,1,2,3,9), type="count_cnv")
  jump_coordinate_10mb = data.frame(code=1:3, count=c(1,3,8), type="jump")
  size_coordinate = data.frame(code=1:10, count=c(171900,615800,2178800,7348000,24320926,53631700,66240800,107971048,144264141,249137896), type="size")
  
  # --- 1. Segment size bins ---
  size_all = data.frame(sample=df$sample, seg_size=df$end-df$start)
  breakpoints = c(-Inf, size_coordinate$count)
  size_all$size_code = cut(size_all$seg_size, breaks=breakpoints, labels=size_coordinate$code, right=TRUE)
  size_tab = as.data.frame.matrix(table(size_all$sample, size_all$size_code))
  colnames(size_tab) = paste0("size_cnv_", colnames(size_tab))
  size_tab$sample = rownames(size_tab)
  
  # --- 2. Absolute CNV state bins ---
  count_cnv = data.frame(sample=df$sample, count=df$major)
  count_cnv$count = pmin(as.numeric(as.character(count_cnv$count)), 9)
  breakpoints_cn = c(-1, count_cnv_coordinate$count)
  count_cnv$cnv_count_code = cut(count_cnv$count, breaks=breakpoints_cn, labels=count_cnv_coordinate$code, right=TRUE)
  count_tab = as.data.frame.matrix(table(count_cnv$sample, count_cnv$cnv_count_code))
  colnames(count_tab) = paste0("count_cnv_", colnames(count_tab))
  count_tab$sample = rownames(count_tab)
  
  # --- 3. CN jumps between adjacent segments ---
  jump_list = list()
  for (s in unique(df$sample)) {
    dsub = df[df$sample == s, ]
    dsub = dsub[order(dsub$Chrom, dsub$start), ]
    dsub$jump = c(NA, abs(diff(as.numeric(as.character(dsub$major)))))
    jump_list[[s]] = data.frame(sample=dsub$sample[-1], jump=dsub$jump[-1])
  }
  jump_all = do.call(rbind, jump_list)
  breakpoints_jump = c(-Inf, jump_coordinate_10mb$count)
  jump_all$size_code = cut(jump_all$jump, breaks=breakpoints_jump, labels=jump_coordinate_10mb$code, right=TRUE)
  jump_tab = as.data.frame.matrix(table(jump_all$sample, jump_all$size_code))
  colnames(jump_tab) = paste0("jump_", colnames(jump_tab))
  jump_tab$sample = rownames(jump_tab)
  
  # --- Merge all features into a matrix ---
  merge_list = list(size_tab, count_tab, jump_tab)
  feature_mat = Reduce(function(x,y) merge(x, y, by="sample", all=TRUE), merge_list)
  rownames(feature_mat) = feature_mat$sample
  feature_mat2 = feature_mat[, -which(colnames(feature_mat)=="sample")]
  feature_mat2[is.na(feature_mat2)] = 0
  feature_mat2 = as.data.frame(lapply(feature_mat2, as.numeric))
  rownames(feature_mat2) = rownames(feature_mat)
  return(feature_mat2)
}




em_signatures_CNV=function(sigts.defn, mut.freqs, max.iter, dbg) {
  library(base)
  nof.sigts.defn=ncol(sigts.defn)
  alpha = stats::runif(nof.sigts.defn); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (i in 1:max.iter) {
    contr = t(array(alpha, dim=c(nof.sigts.defn,28))) * sigts.defn
    probs = contr/array(rowSums(contr), dim=dim(contr))
    probs[is.na(probs)] = 0
    probs = probs * mut.freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  # spit(dbg, "em: exit iteration: %d", i)
  return( alpha )
}

find_best_cutoff=function(time, event, value, cutoffs) {
  results=sapply(cutoffs, function(cutoff) {
    # Dichotomize the variable based on the cutoff
    dichotomized=ifelse(value > cutoff, 1, 0)
    
    # Fit the Cox model with the dichotomized variable
    model=coxph(Surv(time, event) ~ dichotomized)
    
    # Extract the log-likelihood or any other metric (e.g., p-value)
    logLik(model)
  })
  
  # Return the cutoff with the best (e.g., highest log-likelihood)
  best_cutoff=cutoffs[which.max(results)]
  list(best_cutoff = best_cutoff, results = results)
}
