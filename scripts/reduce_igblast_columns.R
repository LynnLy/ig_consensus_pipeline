# Title     : filter_and_prep_for_smolecule.R
# Objective : Find reads with enough full length sequences and cell barcode support
# Created by: lly

library(argparser)

p <- arg_parser("Subset IMGT output to fewer columns")
p <- add_argument(p, "--output", help=".tsv with only selected columns")
p <- add_argument(p, "--input", help="IMGT output in AIRR Rearrangement Format")
argv <- parse_args(p)

library(tidyverse)
library(data.table)
library(LaF)

reduced_cols <- c("sequence_id", "locus", "stop_codon", "vj_in_frame",
                  "v_frameshift","productive", "complete_vdj",
                  "v_call", "d_call", "j_call", "c_call", "v_identity", "c_identity",
                  "sequence", "v_sequence_start", "v_sequence_end", "d_sequence_start",
                  "d_sequence_end", "j_sequence_start", "j_sequence_end",
                  "c_sequence_start", "c_sequence_end", "cdr1_start", "cdr1_end",
                  "cdr2_start", "cdr2_end", "cdr3_start", "cdr3_end", "cell_barcode",
                  "cdr3", "cdr3_aa", "sequence"
)

model <- detect_dm_csv(argv$input, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=10000)
outfile <- argv$output

# Write the column names to the file
data.laf <- laf_open(model)
cols <- next_block(data.laf, nrows=1) %>% as.data.table()
write.table(cols[0] %>% select(reduced_cols), outfile, quote=FALSE, sep="\t", row.names=FALSE)

# Process data in batches (too large to load into memory) and write filtered results to outfile
data.laf <- laf_open(model)
process_data <- function(data, prev) {
    write.table(data %>% select(reduced_cols),
                outfile, na="", quote=FALSE, append=TRUE, col.names=FALSE, sep="\t", row.names=FALSE)
}
process_blocks(data.laf, process_data, nrows=5000000)
