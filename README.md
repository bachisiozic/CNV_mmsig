# CNV_mmsig

Copy-number signature fitting function for fitting  based on “Copy number signatures predict chromothripsis and clinical outcomes in newly diagnosed multiple myeloma” by Maclachlan et al.

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
mat <- generate_cn_feature_matrix(
  segmentation_obj = cnv_mmrf_sel,
  chrom_length_file = "~/path/to/Ref/chrom_sizes_hg38.txt",
  cytoband_file = "~/path/to/Ref/cytoBand_hg38.txt"
)
head(mat)
head(hdp_final)

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
