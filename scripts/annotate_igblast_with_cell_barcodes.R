# Title     : filter_and_prep_for_smolecule.R
# Objective : Find reads with enough full length sequences and cell barcode support
# Created by: lly

library(argparser)

p <- arg_parser("Add cell barcodes to IGBlast output")
p <- add_argument(p, "--output", help="New AIRR Rearrangement File with cell_barcode column")
p <- add_argument(p, "--cell_barcodes", help= ".tsv containing 2 columns: readname and CB:Z: tag ")
p <- add_argument(p, "--input", help="IMGT output in AIRR Rearrangement Format")
argv <- parse_args(p)

library(tidyverse)
library(data.table)
library(LaF)

fread_tags <- function(path) {
  tags <- fread(path, fill=TRUE, header=FALSE, col.names=c("sequence_id", "cell_barcode")) %>%
    unique(by="sequence_id") %>%
    as.data.table()

  # Remove the type from the tags
  tags[, cell_barcode := gsub("CB:Z:", "", cell_barcode)]
  return(tags)
}

cell_barcodes <- fread_tags(argv$cell_barcodes)
model <- detect_dm_csv(argv$input, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=10000)
outfile <- argv$output

# Write the column names to the file
data.laf <- laf_open(model)
cols <- next_block(data.laf, nrows=1) %>% as.data.table()
cols[, cell_barcode := "placeholder"]

write.table(cols[0], outfile, quote=FALSE, sep="\t", row.names=FALSE)

# Process data in batches (too large to load into memory) and write filtered results to outfile
data.laf <- laf_open(model)
print("Processing Data")
process_data <- function(data, prev) {
    output <- cell_barcodes[as.data.table(data), on="sequence_id"]
    write.table(output %>% select(names(cols)),
                outfile, na="", quote=FALSE, append=TRUE, col.names=FALSE, sep="\t", row.names=FALSE)
    }
process_blocks(data.laf, process_data, nrows=1000000)
