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




collapse_segments2 <- function(df) {
  df = df[order(df$chr, df$start), ]
  collapsed = list()
  for (chr in unique(df$chr)) {
    chr_df = df[df$chr == chr, ]
    chr_segments = list()
    for (i in 1:nrow(chr_df)) {
      if (length(chr_segments) == 0) {
        chr_segments[[1]] = chr_df[i, ]
      } else {
        if (chr_segments[[length(chr_segments)]]$end + 1 >= chr_df$start[i]) {
          chr_segments[[length(chr_segments)]]$end = chr_df$end[i]
        } else {
          chr_segments[[length(chr_segments) + 1]] = chr_df[i, ]
        }
      }
    }
    chr_collapsed <- do.call(rbind, chr_segments)
    collapsed <- c(collapsed, list(chr_collapsed))
  }
  collapsed_df <- do.call(rbind, collapsed)
  return(collapsed_df)
}

generate_cn_feature_matrix <- function(
    segmentation_obj,     # path to CNV segmentation (tab-delimited: sample, Chrom, start, end, major, minor)
    chrom_length_file,     # path to chromosome sizes (tab-delimited, Offset=chrom number, Size=length)
    cytoband_file          # path to cytoband file (UCSC cytoBand, or similar)
) {
  # 1. Reference data: Chromosome sizes (assume hg38-like)
  snps <- read.table(chrom_length_file, header=TRUE, sep="\t")
  snps$V2 <- paste0("chr", snps$Offset)
  levels(snps$V2) <- c(levels(snps$V2), "chr23","chr24")
  snps$V2[snps$V2 == "chrX"] <- "chr23"
  snps$V2[snps$V2 == "chrY"] <- "chr24"
  
  # 2. Make 10Mb bins
  cyto_10mb <- list()
  for  (i in c(1:22,"X")) {
    chr_filter = i
    limit = snps[snps$Offset == chr_filter, ]
    post = limit$Size
    x = seq(from=0, to=(post-10000000), by=10000000)
    p = seq(from=10000000, to=post, by=10000000)
    alfa = length(p)
    tail = c(p[alfa], post)
    y = data.frame(x, p)
    y = rbind(y, tail)
    y$chr = rep(chr_filter, nrow(y))
    colnames(y) = c("start", "end", "chr")
    cyto_10mb[[chr_filter]] = y
  }
  cyto_10mb_all = do.call("rbind", cyto_10mb)
  gr1_10mb = with(cyto_10mb_all, GRanges(chr, IRanges(start=start, end=end)))
  
  # 3. Create centromere coordinates (collapsed_df)
  band = read.delim(cytoband_file, header=FALSE)
  band2 = rbind.data.frame(
    band[grep("var", band$V5),],
    band[grep("acen", band$V5),]
  )
  band2 = band2[order(band2$V1, band2$V2),]
  band2 = band2[band2$V2 != "chrY",]
  colnames(band2) <- c("chr","start","end","band","code")
  collapsed_df = collapse_segments2(band2)
  collapsed_df$size = collapsed_df$end - collapsed_df$start
  
  # 4. Read segmentation file (assume header: sample, Chrom, start, end, major, minor)
  #cnv_mmrf_sel = read.table(segmentation_file, header=TRUE, sep="\t")
  cnv_mmrf_sel = segmentation_obj
  colnames(cnv_mmrf_sel) <- c("sample","Chrom","start","end","major","minor")
  cnv_mmrf_sel$major[cnv_mmrf_sel$major < 0] <- 0
  cnv_mmrf_sel$minor[cnv_mmrf_sel$minor < 0] <- 0
  cnv_mmrf_sel$major[cnv_mmrf_sel$major == 0.5] <- 0
  cnv_mmrf_sel$major[cnv_mmrf_sel$major == 1.5] <- 1
  cnv_mmrf_sel$major[cnv_mmrf_sel$major == 2.5] <- 3
  cnv_mmrf_sel$major[cnv_mmrf_sel$major == 3.5] <- 4
  cnv_mmrf_sel$major[cnv_mmrf_sel$major == 4.5] <- 5
  cnv_mmrf_sel$major[cnv_mmrf_sel$major == 5.5] <- 5
  cnv_mmrf_sel$minor[cnv_mmrf_sel$minor == 0.5] <- 0
  cnv_mmrf_sel$minor[cnv_mmrf_sel$minor == 1.5] <- 1
  cnv_mmrf_sel$minor[cnv_mmrf_sel$minor == 2.5] <- 3
  cnv_mmrf_sel$minor[cnv_mmrf_sel$minor == 3.5] <- 4
  cnv_mmrf_sel$minor[cnv_mmrf_sel$minor == 4.5] <- 5
  cnv_mmrf_sel$minor[cnv_mmrf_sel$minor == 5.5] <- 5
  cnv_mmrf_sel$Chrom <- gsub("chr", "", cnv_mmrf_sel$Chrom)
  cnv_mmrf_sel = cnv_mmrf_sel[order(cnv_mmrf_sel$sample, cnv_mmrf_sel$Chrom, cnv_mmrf_sel$start),]
  cnv_mmrf_sel$code_row = 1:nrow(cnv_mmrf_sel)
  
  # Remove IG loci and chrX
  igh_cnv = cnv_mmrf_sel[
    (cnv_mmrf_sel$Chrom == 14 & cnv_mmrf_sel$start > 106032614 & cnv_mmrf_sel$end < 108288051) |
      (cnv_mmrf_sel$Chrom == 22 & cnv_mmrf_sel$start > 21080474 & cnv_mmrf_sel$end < 26065085) |
      (cnv_mmrf_sel$Chrom == 2  & cnv_mmrf_sel$start > 87090568  & cnv_mmrf_sel$end < 93274235),]
  cnv_mmrf_no_igh = cnv_mmrf_sel[!cnv_mmrf_sel$code_row %in% igh_cnv$code_row,]
  cnv_mmrf_no_igh_no_x = cnv_mmrf_no_igh[cnv_mmrf_no_igh$Chrom != "X",]
  cnv_mmrf_no_igh_no_x$major <- as.numeric(as.character(cnv_mmrf_no_igh_no_x$major))
  cnv_mmrf_no_igh_no_x$minor <- as.numeric(as.character(cnv_mmrf_no_igh_no_x$minor))
  
  # Collapse adjacent segments with same major
  collapse_adjacent <- function(df) {
    out = list()
    sample_list = unique(df$sample)
    for(j in seq_along(sample_list)) {
      cna_mmrf_single = df[df$sample == sample_list[j], ]
      sam_cnv_list = list()
      chr_list = unique(cna_mmrf_single$Chrom)
      for(i in seq_along(chr_list)) {
        cna_mmrf_single_chr = cna_mmrf_single[cna_mmrf_single$Chrom == chr_list[i], ]
        list_chr = list()
        vec = rle((paste(cna_mmrf_single_chr$major)))$length
        for(w in seq_along(vec)) {
          if(w == 1) {
            int = cna_mmrf_single_chr[1:vec[w], ]
            cna_mmrf_single_row = c(int$sample[1], int$Chrom[1], int$start[1], int$end[nrow(int)], int$major[1], max(int$minor))
          } else {
            int = cna_mmrf_single_chr[(sum(vec[1:(w-1)])+1):sum(vec[1:(w)]), ]
            cna_mmrf_single_row = c(int$sample[1], int$Chrom[1], int$start[1], int$end[nrow(int)], int$major[1], max(int$minor))
          }
          list_chr[[w]] = cna_mmrf_single_row
        }
        list_chr2 = do.call("rbind", list_chr)
        sam_cnv_list[[i]] = list_chr2
      }
      sam_cnv_list2 = do.call("rbind", sam_cnv_list)
      out[[j]] = sam_cnv_list2
    }
    result = do.call("rbind", out)
    result = as.data.frame(result)
    colnames(result) = c("sample", "Chrom", "start", "end", "major", "minor")
    result$end = as.numeric(as.character(result$end))
    result$start = as.numeric(as.character(result$start))
    return(result)
  }
  cnv_mmrf_no_igh_no_x2 = collapse_adjacent(cnv_mmrf_no_igh_no_x)
  cnv_mmrf_final2 = cnv_mmrf_no_igh_no_x2[(cnv_mmrf_no_igh_no_x2$end - cnv_mmrf_no_igh_no_x2$start) > 50000, ]
  bad_artifacts = cnv_mmrf_final2[
    (cnv_mmrf_final2$end-cnv_mmrf_final2$start) < 1000000 &
      cnv_mmrf_final2$minor == 0 & cnv_mmrf_final2$major > 2, ]
  cnv_mmrf_final = cnv_mmrf_final2[
    !paste(cnv_mmrf_final2$sample, cnv_mmrf_final2$Chrom, cnv_mmrf_final2$start) %in%
      paste(bad_artifacts$sample, bad_artifacts$Chrom, bad_artifacts$start), ]
  cnv_mmrf = collapse_adjacent(cnv_mmrf_final)
  cnv_mmrf$sample = as.character(as.character(cnv_mmrf$sample))
  cnv_mmrf$Chrom = as.character(as.character(cnv_mmrf$Chrom))
  cnv_mmrf[,3:ncol(cnv_mmrf)] = apply(cnv_mmrf[,3:ncol(cnv_mmrf)], 2, function(x){as.numeric(as.character(x))})
  cnv_mmrf = cnv_mmrf[cnv_mmrf$Chrom != "23",]
  cnv_mmrf = cnv_mmrf[!is.na(cnv_mmrf$major),]
  
  # ---- CN feature extraction ----
  
  # Coordinate reference cutoffs
  count_cnv_coordinate_10mb = data.frame(code = c(1,3,4), count = c(3,6,31), type="10mb")
  count_cnv_coordinate = data.frame(code = 1:5, count = c(0,1,2,3,9), type="count_cnv")
  jump_coordinate_10mb = data.frame(code = 1:3, count = c(1,3,8), type="jump")
  band_coordinate_new = data.frame(code = 1:3, count = c(5,17,60), type="band")
  osci_coordinate = data.frame(code = 1:4, count = c(1,4,9,38), type="osci")
  size_coordinate = data.frame(
    code = 1:10,
    count = c(171900,615800,2178800,7348000,24320926,53631700,66240800,107971048,144264141,249137896),
    type="size"
  )
  
  # ------------- Generate all features ---------------
  
  # Get sample list
  sample_list = unique(cnv_mmrf$sample)
  
  # Prepare result lists
  size_all = list()
  count_10mb_all = list()
  count_jump_all = list()
  count_cnv_all = list()
  band_rate_all = list()
  osci_all = list()
  
  for(w in seq_along(sample_list)) {
    cnv_mmrf2 = cnv_mmrf[cnv_mmrf$sample == sample_list[w], ]
    
    # Segment size
    cnv_mmrf2$seg_size = (cnv_mmrf2$end - cnv_mmrf2$start)
    size_all[[w]] = cnv_mmrf2[,c("sample","seg_size")]
    
    # Breaks per 10Mb
    cnv_temp_brk = cnv_mmrf2[,c(1,2,3,5,6,7)]
    cnv_mmrf2_second = cnv_mmrf2[,c(1,2,4,5,6,7)]
    int_dipl = as.data.frame.matrix(table(cnv_mmrf2_second$Chrom, cnv_mmrf2_second$major))
    diploid_chr = rownames(int_dipl[int_dipl$`2`==1 & rowSums(int_dipl)==1,])
    cnv_temp_brk = cnv_mmrf2_second[!cnv_mmrf2_second$Chrom %in% diploid_chr,]
    gr_cna_comm = with(cnv_temp_brk, GRanges(Chrom, IRanges(start=end, end=end)))
    values(gr_cna_comm) = DataFrame(sample = cnv_temp_brk$sample, major = cnv_temp_brk$major, minor= cnv_temp_brk$minor, seg_size= cnv_temp_brk$seg_size)
    range_10mb = merge(as.data.frame(gr_cna_comm),as.data.frame(gr1_10mb),by="seqnames",suffixes=c("A","B"))
    range_dri_10mb = range_10mb[with(range_10mb, startB <= startA & endB >= endA),]
    count_brk_10mb = as.data.frame(table(paste(range_dri_10mb$seqnames, range_dri_10mb$startB)))
    cyto_10mb_all$Var1 = paste(cyto_10mb_all$chr, cyto_10mb_all$start)
    count_10mb_file = join(cyto_10mb_all, count_brk_10mb, by="Var1")
    count_10mb_file[is.na(count_10mb_file)] = 0
    count_10mb_file$sample = sample_list[w]
    count_10mb_df = count_10mb_file[,c("sample","Freq")]
    colnames(count_10mb_df)[2] = "count"
    count_10mb_all[[w]] = count_10mb_df
    
    # Copy number jumps
    cnv_mmrf2_second$jump = NA
    cnv_mmrf2_second_jump = cnv_mmrf2_second
    if(length(unique(cnv_mmrf2_second_jump$Chrom)) != 0) {
      chr_list_jump = unique(cnv_mmrf2_second_jump$Chrom)
      cnv_mmrf2_second_jump$jump = NA
      all_chr_jump = list()
      for(jj in chr_list_jump) {
        cnv_mmrf2_second_jump_int = cnv_mmrf2_second_jump[cnv_mmrf2_second_jump$Chrom == jj, ]
        for(z in 1:nrow(cnv_mmrf2_second_jump_int)) {
          if(z == 1) {
            cnv_mmrf2_second_jump_int$jump[1] = NA
          } else {
            cnv_mmrf2_second_jump_int$jump[z] = abs((cnv_mmrf2_second_jump_int$major[z]) - (cnv_mmrf2_second_jump_int$major[z-1]))
          }
        }
        all_chr_jump[[jj]] = cnv_mmrf2_second_jump_int
      }
      all_chr_jump2 = do.call("rbind", all_chr_jump)
    } else {
      all_chr_jump2 = cnv_mmrf2_second[1,]
      all_chr_jump2$jump = 0
    }
    all_chr_jump2 = all_chr_jump2[! is.na(all_chr_jump2$jump), ]
    temp_jump = all_chr_jump2[,c("sample","jump")]
    colnames(temp_jump) = c("sample","count")
    count_jump_all[[w]] = temp_jump
    
    # Copy number of each segment
    count_cnv_final_df = cnv_mmrf2[,c("sample","major")]
    colnames(count_cnv_final_df)[2] = "count"
    count_cnv_all[[w]] = count_cnv_final_df[,c("sample","count")]
    
    # Breaks per chromosome arm
    chrom_arms = collapsed_df[collapsed_df$start != 0, ]
    chrom_arms$chrom = gsub("chr", "", chrom_arms$chr)
    cnv_mmrf2_second = cnv_mmrf2[,c(1,2,4,5,6,7)]
    cnv_temp_brk_arm = cnv_mmrf2_second
    if(nrow(cnv_temp_brk_arm) != 0) {
      gr_cna_comm = with(cnv_temp_brk_arm, GRanges(Chrom, IRanges(start=end, end=end)))
      values(gr_cna_comm) = DataFrame(sample = cnv_temp_brk_arm$sample,
                                      major = cnv_temp_brk_arm$major, minor= cnv_temp_brk_arm$minor, seg_size= cnv_temp_brk_arm$seg_size)
      gr_band = with(chrom_arms, GRanges(chrom, IRanges(start=start, end=end)))
      range_arm = merge(as.data.frame(gr_cna_comm),as.data.frame(gr_band),by="seqnames",suffixes=c("A","B"))
      range_arm$arm = NA
      range_arm$arm[range_arm$startA > range_arm$endB] = "q_arm"
      range_arm$arm[range_arm$startA < range_arm$startB] = "p_arm"
      range_arm$arm[range_arm$startB <= range_arm$startA & range_arm$endB >= range_arm$startA] = "centro"
      db_arm_counts = as.data.frame(table(paste(range_arm$seqnames, range_arm$arm)))
    } else {
      db_arm_counts = matrix(c("13 q_arm", 0), nrow=1)
      db_arm_counts = as.data.frame(db_arm_counts)
      colnames(db_arm_counts) = c("Var1","Freq")
    }
    file_int_band = as.data.frame(c(paste0(c(1:22), (" p_arm")), paste0(c(1:22), (" q_arm"))))
    colnames(file_int_band) = "Var1"
    band_rate = join(file_int_band, db_arm_counts, by="Var1")
    band_rate[is.na(band_rate)] = 0
    band_rate$sample = sample_list[w]
    colnames(band_rate)[2] = "count"
    band_rate_all[[w]] = band_rate[,c("sample","count")]
    
    # Oscillating CN
    out = c()
    chrs = unique(cnv_mmrf2$Chrom)
    cnv_mmrf2$tot = cnv_mmrf2$major
    oscCounts = c()
    for(c in chrs) {
      currseg = cnv_mmrf2[cnv_mmrf2$Chrom == c,"tot"]
      currseg = round(as.numeric(currseg))
      if(length(currseg) > 3) {
        prevval = currseg[1]
        count = 0
        for(j in 3:length(currseg)) {
          if(j == length(currseg)) {
            oscCounts = rbind(oscCounts, c(c,count))
            count = 0
          } else {
            if(abs(currseg[j] - prevval) <= 1 & currseg[j] != currseg[j-1]) {
              count = count + 1
            } else {
              oscCounts = rbind(oscCounts, c(c,count))
              count = 0
            }
          }
          prevval = currseg[j-1]
        }
      } else {
        oscCounts = rbind(oscCounts, c(c,0))
      }
    }
    oscCounts_df = as.data.frame(oscCounts)
    oscCounts_df$sample = sample_list[w]
    oscCounts_df$V2 = as.numeric(as.character(oscCounts_df$V2))
    osci_all[[w]] = oscCounts_df
  }
  
  # ----- Bin features and build final matrix -----
  
  size_all22 = do.call("rbind", size_all)
  size_all22 = as.data.frame(size_all22)
  size_all22$seg_size = as.numeric(as.character(size_all22$seg_size))
  breakpoints = c(-Inf, size_coordinate$count)
  values_classified = cut(size_all22$seg_size,
                          breaks = breakpoints,
                          labels = size_coordinate$code,
                          right = TRUE)
  size_all22$size_code = values_classified
  
  count_cnv_all2 = do.call("rbind", count_cnv_all)
  count_cnv_all2$count[count_cnv_all2$count > 9] = 9
  count_cnv_all2$cnv_count_code = NA
  count_cnv_all2$cnv_count_code[count_cnv_all2$count == 0] = 1
  count_cnv_all2$cnv_count_code[count_cnv_all2$count == 1] = 2
  count_cnv_all2$cnv_count_code[count_cnv_all2$count == 2] = 3
  count_cnv_all2$cnv_count_code[count_cnv_all2$count == 3] = 4
  count_cnv_all2$cnv_count_code[count_cnv_all2$count >= 4] = 5
  
  count_10mb_all2 = do.call("rbind", count_10mb_all)
  count_10mb_all2 = as.data.frame(count_10mb_all2)
  count_10mb_all2$count = as.numeric(as.character(count_10mb_all2$count))
  count_10mb_all2_alt = count_10mb_all2[count_10mb_all2$count != 0, ]
  breakpoints = c(-Inf, count_cnv_coordinate_10mb$count)
  values_classified = cut(count_10mb_all2_alt$count,
                          breaks = breakpoints,
                          labels = count_cnv_coordinate_10mb$code,
                          right = TRUE)
  count_10mb_all2_alt$size_code = values_classified
  
  count_jump_all2 = do.call("rbind", count_jump_all)
  count_jump_all2 = as.data.frame(count_jump_all2)
  count_jump_all2$count = as.numeric(as.character(count_jump_all2$count))
  count_jump_all2_alt = count_jump_all2
  breakpoints = c(-Inf, jump_coordinate_10mb$count)
  values_classified = cut(count_jump_all2_alt$count,
                          breaks = breakpoints,
                          labels = jump_coordinate_10mb$code,
                          right = TRUE)
  count_jump_all2_alt$size_code = values_classified
  
  band_rate_all2 = do.call("rbind", band_rate_all)
  band_rate_all2 = as.data.frame(band_rate_all2)
  band_rate_all2$Freq = as.numeric(as.character(band_rate_all2$count))
  band_rate_all2_alt = band_rate_all2[band_rate_all2$count != 0, ]
  breakpoints = c(-Inf, band_coordinate_new$count)
  values_classified = cut(band_rate_all2_alt$count,
                          breaks = breakpoints,
                          labels = band_coordinate_new$code,
                          right = TRUE)
  band_rate_all2_alt$size_code = values_classified
  band_rate_all2_alt = band_rate_all2_alt[, -3]
  
  osci_all2 = do.call("rbind", osci_all)
  osci_all2 = as.data.frame(osci_all2[,c(3,2)])
  colnames(osci_all2)[2] = "count"
  osci_all2$count = as.numeric(as.character(osci_all2$count))
  osci_all2_alt = osci_all2
  breakpoints = c(-Inf, osci_coordinate$count)
  values_classified = cut(osci_all2_alt$count,
                          breaks = breakpoints,
                          labels = osci_coordinate$code,
                          right = TRUE)
  osci_all2_alt$size_code = values_classified
  
  # Final tables (wide)
  osci_tab = as.data.frame.matrix(table(osci_all2_alt$sample, osci_all2_alt$size_code))
  colnames(osci_tab) = paste0("osci_", colnames(osci_tab))
  osci_tab$sample = rownames(osci_tab)
  band_tab = as.data.frame.matrix(table(band_rate_all2_alt$sample, band_rate_all2_alt$size_code))
  colnames(band_tab) = paste0("band_", colnames(band_tab))
  band_tab$sample = rownames(band_tab)
  jump_tab = as.data.frame.matrix(table(count_jump_all2_alt$sample, count_jump_all2_alt$size_code))
  colnames(jump_tab) = paste0("jump_", colnames(jump_tab))
  jump_tab$sample = rownames(jump_tab)
  mb_10_tab = as.data.frame.matrix(table(count_10mb_all2_alt$sample, count_10mb_all2_alt$size_code))
  colnames(mb_10_tab) = paste0("mb_10_", colnames(mb_10_tab))
  mb_10_tab$sample = rownames(mb_10_tab)
  count_tab = as.data.frame.matrix(table(count_cnv_all2$sample, count_cnv_all2$cnv_count_code))
  colnames(count_tab) = paste0("count_cnv_", colnames(count_tab))
  count_tab$sample = rownames(count_tab)
  size_tab = as.data.frame.matrix(table(size_all22$sample, size_all22$size_code))
  colnames(size_tab) = paste0("size_cnv_", colnames(size_tab))
  size_tab$sample = rownames(size_tab)
  
  hdp_final = Reduce(merge, list(mb_10_tab, count_tab, jump_tab, band_tab, osci_tab, size_tab))
  rownames(hdp_final) = hdp_final$sample
  hdp_final2 = hdp_final[,-which(colnames(hdp_final) == "sample")]
  hdp_final2[,1:ncol(hdp_final2)] = apply(hdp_final2[,1:ncol(hdp_final2)], 2, function(x){as.numeric(as.character(x))})
  
  return(hdp_final2)
}
                       
