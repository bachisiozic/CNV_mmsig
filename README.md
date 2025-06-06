# CNV_mmsig

Copy number variant signature fitting function for fitting based on “Copy number signatures predict chromothripsis and clinical outcomes in newly diagnosed multiple myeloma” by Maclachlan et al.

---

## Overview

- **Input**: Segmentation data, chromosome length, cytoband references, and a CN signature reference file.
- **Process**: 
  - Generate a CN feature matrix from segmentation data.
  - Fit known CN signatures to each sample using a custom EM function.
- **Output**: Relative CN signature contributions for each sample.

---

## Usage

### 1. Prepare Input Files

- Segmentation object: e.g., `cnv_mmrf_sel`
- Chromosome length file: `chrom_sizes_hg38.txt`
- Cytoband file: `cytoBand_hg38.txt`
- CN signature reference: e.g., `CNV_SIGNATURES_PROFILES.txt`

### 2. Run the Script

```r
library(GenomicRanges)
library(plyr)
library(dplyr)

source("~/path/to/CNV_mmsig.R")

cnv_seg=read.delim("~/path/to/your/cnv_seg.txt")
colnames(cnv_seg)=c("sample","Chrom","start","end","major","minor")
cnv_seg=cnv_seg[order(cnv_mmrf_sel$sample, cnv_mmrf_sel$Chrom, cnv_mmrf_sel$start),]
cnv_seg$Chrom=gsub("chr","", cnv_seg$Chrom)
cnv_seg$code_row=1:nrow(cnv_seg)
head(cnv_seg)
#sample   Chrom    start       end major minor code_row
#Sample_1     1   181508   1979102     2     1        1
#Sample_2     1   791595 248946322     2    NA        2
#Sample_3     1  1979102   1982826     3     0        3
#Sample_4     1  1982826  12834382     2     1        4
#Sample_5     1 12834382  13367518     2     0        5
#Sample_6     1 13367518  16546038     2     1        6

mat=generate_cn_feature_matrix(
  segmentation_obj = cnv_seg,
  chrom_length_file = "~/path/to/Ref/chrom_sizes_hg38.txt",
  cytoband_file = "~/path/to/Ref/cytoBand_hg38.txt"
)
head(mat)
#         mb_10_1 mb_10_3 mb_10_4 count_cnv_1 count_cnv_2 count_cnv_3 count_cnv_4 count_cnv_5 jump_1 jump_2 jump_3 band_1 band_2 band_3 osci_1 osci_2 osci_3 osci_4
#Sample_1      26       1       0           0           5          25          13           7     24      4      0     25      1      0     24      0      0      0
#Sample_2       1       0       0           1           1          22           0           0      1      1      0     22      0      0     22      0      0      0
#Sample_3      18       0       0           1           4          10          13           0      2      4      0     24      0      0     22      0      0      0
#Sample_4       8       0       0           1           0          15           4           3      0      1      0     22      0      0     22      0      0      0
#Sample_5      16       0       0           1           4          13          11           0      2      5      0     25      0      0     22      0      0      0
#Sample_6      26       0       0           0           0          14          13           9      3      9      2     24      1      0     21      0      1      0
#            size_cnv_1 size_cnv_2 size_cnv_3 size_cnv_4 size_cnv_5 size_cnv_6 size_cnv_7 size_cnv_8 size_cnv_9 size_cnv_10
#Sample_1          0          1          5          3          4          7          6         12          5           7
#Sample_2          1          0          1          0          0          2          2          6          4           8
#Sample_3          0          0          1          0          3          5          2          6          4           7
#Sample_4          0          0          1          0          0          2          2          6          4           8
#Sample_5          0          0          2          0          4          4          2          7          4           6
#Sample_6          1          0          0          5          6          3          3          6          6           6


ref = read.delim("~/path/to/Ref/CNV_SIGNATURES_PROFILES.txt", stringsAsFactors = FALSE)
colnames(ref) = paste0("CNV_SIG_", 1:5)
rownames(ref) = colnames(hdp_final)[2:29]

all_cnv_sig = list()
for (i in (1:nrow(mat))) {
  max.em.iter = 1000
  dbg = FALSE
  sample.consigts.defn = ref
  sample.mut.freqs = mat[i, ]
  alpha = em_signatures_CNV(
    sigts.defn = sample.consigts.defn,
    mut.freqs = as.numeric(as.character(sample.mut.freqs)),
    max.iter = max.em.iter,
    dbg = dbg
  )
  all_cnv_sig[[i]] = c(rownames(mat)[i], alpha)
}
all_cnv_sig2 = do.call(rbind.data.frame, all_cnv_sig)
```
