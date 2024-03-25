library(argparser)

p <- arg_parser("Filter and prep reads for medaka smolecule")
p <- add_argument(p, "--input", help="IGBlast output in AIRR format")
p <- add_argument(p, "--output", help="Filtered IGBLast output in AIRR format")
p <- add_argument(p, "--strictness", help="Strictness of filters [lenient, strict, consensus]")
p <- add_argument(p, "--c_lengths", help="File containing 2 cols: C gene name and length")
p <- add_argument(p, "--v_lengths", help="File containing 2 cols: V gene name and length")

argv <- parse_args(p)

library(tidyverse)
library(data.table)
library(LaF)

add_cols <- function(data) {
    # Make sure these columns were classified as numeric
   print("Converting to numeric")
   numeric_cols <- c(
      "v_support", "j_support", "c_support", "c_germline_end", "c_germline_start",
      "v_germline_start", "v_germline_end", "c_sequence_start", "c_sequence_end",
      "v_sequence_start", "v_sequence_end", "v_identity", "d_identity", "j_identity",
      "c_identity"
    )
   data[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]

    # How much "excess" sequence is there other than the VDJC?
    # There should be a (not well-defined) 5' and 3' UTR region
    data[, start_excess_seq := substring(sequence, 1, v_sequence_start)]
    data[, end_excess := nchar(sequence) - c_sequence_end]
    data[, end_excess_seq := substring(sequence, c_sequence_end, nchar(sequence))]
    data[, transcript_only := substring(sequence, v_sequence_start, c_sequence_end)]
    # For the purpose of finding V/C gene length, only use the first call in an ambiguous call
    data <- data[!is.na(c_call) & !is.na(v_call)]
    if (nrow(data) == 0) {
      return(data)
    }
    data[, c_call_ambi := tstrsplit(c_call, ",", fixed=TRUE)[[1]]]
    data[, v_call_ambi := tstrsplit(v_call, ",", fixed=TRUE)[[1]]]
    data <- merge(data[!is.na(c_call_ambi)], constant_lengths, by.x="c_call_ambi", by.y="c_call") %>% as.data.table()
    data <- merge(data[!is.na(v_call_ambi)], v_lengths, by.x="v_call_ambi", by.y="v_call") %>% as.data.table()
    return(data)
}

apply_pyir_filters <- function(data, strictness="lenient") {
    if (nrow(data) == 0) {return(data)}
    if(strictness == "lenient") {
        data[v_support < 1.0e-6 &
            c_support < 1.0e-6
        ] %>% return()
    } else if(strictness == "strict") {
      # Must also contain nearly-full length V and C genes
      # Additionally, attempts to remove chimeric sequences
      # Here, defined as having too many bp before the V gene or after the C gene
      # Require cdr3_aa length to be appropriate
        data[v_support < 1.0e-6 &
            c_support < 1.0e-6 &
            (v_germline_end - v_germline_start) > v_full_length * 0.90 &
            c_sequence_start - v_sequence_end < 250 &
            v_sequence_start < 400 &
            end_excess < 800 &
            nchar(cdr3_aa) > 3 & nchar(cdr3_aa) < 50
        ][
          (c_germline_end - c_germline_start) > c_full_length * 0.90 |
            (grepl("IGHD", c_call) & c_germline_end - c_germline_start > 400) | # Noted a large population of shorter IGHD seqs, not present for any other isotype
            (grepl("IGHG3", c_call) & c_germline_end - c_germline_start > 1017) # IGHG3 version used by default is the membrane form; 1017 is 0.9 * secreted length
        ] %>% return()
    } else if(strictness == "consensus") {
        # Filter on fields including in-frame codons and productivity
        data[v_support < 1.0e-6 &
            c_support < 1.0e-6 &
            (v_germline_end - v_germline_start) > v_full_length * 0.90 &
            c_sequence_start - v_sequence_end < 250 &
            v_sequence_start < 400 &
            end_excess < 800 &
            nchar(cdr3_aa) > 3 & nchar(cdr3_aa) < 50 &
            stop_codon == "F" &
            vj_in_frame == "T" &
            productive == "T" &
            complete_vdj == "T"
        ][
          (c_germline_end - c_germline_start) > c_full_length * 0.90 |
            (grepl("IGHD", c_call) & c_germline_end - c_germline_start > 400) |
            (grepl("IGHG3", c_call) & c_germline_end - c_germline_start > 1017)
        ] %>% return()
    }
}

constant_lengths <- fread(argv$c_lengths,
                          col.names=c("c_call", "c_full_length"),
                          sep="\t")
v_lengths <- fread(argv$v_lengths,
                   col.names=c("v_call", "v_full_length"),
                   sep="\t")

infile <- argv$input
outfile <- argv$output
model <- detect_dm_csv(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE, nrows=10000)

# Write the column names to the file
data.laf <- laf_open(model)
cols <- next_block(data.laf, nrows=1) %>% as.data.table()
write.table(cols[0], outfile, quote=FALSE, sep="\t", row.names=FALSE)

# Process data in batches (too large to load into memory) and write filtered results to outfile
data.laf <- laf_open(model)
print("Processing data")
process_data <- function(data, prev) {
    data <- apply_pyir_filters(data %>% as.data.table() %>% add_cols(), strictness=argv$strictness)
    write.table(data %>% select(names(cols)),
                outfile, quote=FALSE, append=TRUE,
                col.names=FALSE, sep="\t", row.names=FALSE,
                na=""
                )
}
process_blocks(data.laf, process_data, nrows=5000000)
