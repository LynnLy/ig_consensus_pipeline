# Title     : filter_and_prep_for_smolecule.R
# Objective : Find reads with enough full length sequences and cell barcode support
# Created by: lly

library(argparser)

p <- arg_parser("Filter and prep reads for medaka smolecule")
# p <- add_argument(p, "--folder", help="name of folder")
p <- add_argument(p, "--output", help="tsv with 1) Old readname and 2) New readname: {cell_barcode}_{locus}_{number}")
p <- add_argument(p, "--output_depth", help="tsv with cell barcode and total number of reads")
p <- add_argument(p, "--cell_barcodes", help="tsv containing 2 columns: readname and CB:Z: tag ")
p <- add_argument(p, "--readnames", help="filter on additional readnames")
p <- add_argument(p, "--min_reads", help="Minimum number of reads for that cell", default=3)
p <- add_argument(p, "--locus", help="IGH, IGK, or IGL?")
p <- add_argument(p, "--max_reads", help="Maximum number of reads per consensus; too many reads requires too much RAM when running smolecule", default=150)
argv <- parse_args(p)

library(tidyverse)
library(data.table)

fread_tags <- function(path) {
  tags <- fread(path, fill=TRUE, header=FALSE, col.names=c("readname", "CB")) %>%
    unique(by="readname") %>%
    as.data.table()

  # Remove the type from the tags
  tags[, CB := gsub("CB:Z:", "", CB)]
  return(tags)
}

tags <- fread_tags(argv$cell_barcodes)
readnames <- fread(argv$readnames, header=FALSE)

filtered_tags <- tags[readname %in% readnames$V1]
sufficient_depth <- filtered_tags[, .N, CB][N > argv$min_reads]
fwrite(sufficient_depth, file=argv$output_depth, sep="\t")

renaming_guide <- filtered_tags[CB %in% sufficient_depth$CB][
    , index:=rowid(CB)][
    index <= argv$max_reads][
    , .(readname, newname=paste0(CB, "_", argv$locus, "_", index))]
fwrite(file=argv$output, renaming_guide, sep="\t")
