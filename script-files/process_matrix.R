library(dplyr)
library(tidyverse)
library(jsonlite)

###################################################################################################

args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################

#Parse JSON header for metadata
FIRSTline = readr::read_lines(FILE, n_max = 1)
stopifnot(startsWith(FIRSTline, "@"))
META = jsonlite::fromJSON(sub("^@", "", FIRSTline))
  
SAMPLElabels = META$sample_labels               # e.g. c("siLuc_H3K9me3","siPiwi_H3K9me3")
SAMPLEbounds = META$sample_boundaries           # e.g. c(0, 1001, 2002)  (zero-based)
GROUPlabels = META$group_labels                # e.g. c("cluster_1","cluster_2")
GROUPbounds = META$group_boundaries            # zero-based row boundaries
  
BINsize    <- META[["bin size"]][1]            # 10
US    <- META$upstream[1]                 # 5000
DS  <- META$downstream[1]               # 5000
nBINS       <- US/BINsize + 1 + DS/BINsize  # 1001
BINlabels  <- - (US/BINsize) : (DS/BINsize)  # -500..500
  
#Read the numeric part (skip JSON line)
# Columns: chrom, start, end, name, score, strand, then bins...
RAW <- readr::read_tsv(
  FILE, skip = 1, col_names = FALSE, progress = FALSE,
  col_types = cols(
    X1 = col_character(),  # chrom
    X2 = col_double(),     # start
    X3 = col_double(),     # end
    X4 = col_character(),  # name
    X5 = col_double(),     # score
    X6 = col_character(),  # strand
    .default = col_double()
  )
)
  
#Locate the bin columns per sample (convert deepTools' 0-based boundaries to 1-based + offset for first 6 columns)
#For sample k: columns are (6 + samp_bounds[k] + 1) ... (6 + samp_bounds[k+1])
SAMPLEbins_cols <- map2(
  SAMPLEbounds[-length(SAMPLEbounds)],
  SAMPLEbounds[-1],
  ~ seq.int(6 + .x + 1, 6 + .y)
)
names(SAMPLEbins_cols) <- SAMPLElabels
  
#Build per-sample long tables with explicit bin labels
make_long_for_sample <- function(SAMPLEname) {
  cols <- SAMPLEbins_cols[[SAMPLEname]]
  stopifnot(length(cols) == nBINS)
  
  # set bin labels as column names for this sample’s bin matrix
  vals <- RAW[, cols, drop = FALSE]
  colnames(vals) <- as.character(BINlabels)
  
  tibble(
    CHR  = RAW$X1,
    START  = as.integer(RAW$X2),
    END    = as.integer(RAW$X3),
    NAME   = RAW$X4,
    STRAND = RAW$X6
  ) %>%
    bind_cols(vals) %>%
    pivot_longer(
      cols = all_of(as.character(BINlabels)),
      names_to = "BIN",
      values_to = SAMPLEname
    ) %>%
    mutate(BIN = as.integer(BIN))  # -500..500
}
  
LONGlist = lapply(SAMPLElabels, make_long_for_sample)
  
# Join the samples side-by-side by region keys + bin
TABLE = reduce(
  LONGlist,
  ~ full_join(.x, .y, by = c("CHR","START","END","NAME","STRAND","BIN"))
)
  
# --- 5) Attach k-means cluster per region using group_boundaries (deepTools uses 0-based row indexing)
row_zero_idx = seq_len(nrow(RAW)) - 1L
GROUPindex <- findInterval(row_zero_idx, GROUPbounds, rightmost.closed = TRUE)
REGIONcluster = tibble(
  CHR  = RAW$X1,
  START  = as.integer(RAW$X2),
  END    = as.integer(RAW$X3),
  NAME   = RAW$X4,
  STRAND = RAW$X6,
  CLUSTER = factor(GROUPlabels[pmax(1, pmin(GROUPindex, length(GROUPlabels)))],
                   levels = GROUPlabels)
)

write_tsv(
  REGIONcluster,
  file.path(locTMP, "TEs_H3K9me3_matrix_region_clusters.txt")
)
