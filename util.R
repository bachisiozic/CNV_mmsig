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
