# Parse blastn alignments of IG genes to the assembly, and get the best functional hits
# Annotation accuracy critically assumes that the start of each gene is more than X bp away from the previous gene.
# If this is not true, then two distinct genes may be collapsed and only have 1 best hit between the two of them.
# The distance between start sites varies for the different genes:
# D mingap = 5 ; J mingap = 50 ; C mingap = 200 ; V mingap = 300

library(argparser)

p <- arg_parser("Parse blastn alignments of IG genes to the assembly, and get the best hits")
p <- add_argument(p, "--output_prefix", help="Files will be written to {output_prefix}.best_hits.tsv, {output_prefix}_V.bed, etc.")
p <- add_argument(p, "--input_prefix", help="ie. results/annotations/{sample}.{region}")
p <- add_argument(p, "--locus", help="locus [IGH, IGK, IGL]")
argv <- parse_args(p)

library(tidyverse)
library(data.table)

convert_tags <- function(sam) {
    # Convert some sam tags to numbers
    sam[, AS := gsub("AS:i:", "", AS) %>% as.numeric()]
    sam[, EV := gsub("EV:f:", "", EV) %>% as.numeric()]
    sam[, NM := gsub("NM:i:", "", NM) %>% as.numeric()]
    sam[, PI := gsub("PI:f:", "", PI) %>% as.numeric()]
    sam[, BS := gsub("BS:f:", "", BS) %>% as.numeric()]
}


fread_sam <- function(path) {
    cols <- c("qname", "flag", "rname", "pos",
    "mapq", "cigar", "rnext", "pnext",
     "tlen", "seq", "qual", "AS", "EV", "NM", "PI", "BS")

    skip_sam_header <- paste0("grep -v '^@' ", path)
    sam <- fread(cmd=skip_sam_header, col.names=cols)
    convert_tags(sam)
    return(sam)
}


annotate_anchor_alignments <- function(sam, min_gap=300) {
    # Alignments that begin a new potential IG V position will receive new_group==TRUE
    # If the new alignment is within min_gap bp of the start of the first (anchor) alignment, it will not be a new group
    sam[, prev_pos := shift(pos, 1), by="rname"]
    sam[, diff := pos - prev_pos]
    sam[, new_group := FALSE]
    sam[diff > min_gap | is.na(diff), new_group := TRUE]
    sam[new_group == TRUE, group_end := pos + min_gap]
}


process_sam <- function(sam, min_gap=300) {
    annotate_anchor_alignments(sam, min_gap)
}


roll_alignments <- function(sam) {
    # Annotate each entry with the closest anchor position, to group nearby alignments together
    anchor_alignments <- sam[new_group == TRUE][order(rname, pos)][, .(rname, flag, pos, prev_pos, diff, new_group, group_end)]

    anchor_alignments[, rollpos := pos]
    sam[, rollpos := pos]
    setkey(anchor_alignments, "rname", "rollpos")
    setkey(sam, "rname", "rollpos")

    rolled <- anchor_alignments[
        sam[, .(qname, rname, flag, rollpos, pos, cigar, AS, PI, diff)]
        , roll = TRUE][
        order(rname, pos)
    ]

    setnames(rolled, "pos", "anchor_pos") # Position of the first / "anchor" alignment for this group of alignments
    setnames(rolled, "i.pos", "entry_pos") # Position of the actual alignment

    rolled[, rollpos := NULL]
    return(rolled)
}

get_best_hits <- function(rolled) {
  #  top_readcov <- rolled[]

    top_AS_index <- rolled[, .I[AS == max(AS)], by=c("rname", "anchor_pos")]$V1
    top_AS <- rolled[top_AS_index]
    top_AS[, i.diff := NULL]

    # For positions with multiple hits with equal scores, label all hits except the first as "redundant_hit"
    top_AS[anchor_pos == shift(anchor_pos, 1), redundant_hit := TRUE]
    return(top_AS)
}

process_annotations <- function(fileprefix, segment) {
    if(segment == "D") {
        mingap = 5
    } else if(segment == "J") {
        mingap=50
    } else if(segment == "C") {
        mingap=200
    } else if(segment == "V") {
        mingap=300
    }

    seqkit <- fread(paste0(fileprefix, "_", segment, ".tsv"))
    readlens <- seqkit[, .(qname=Read, ReadLen)][, .N, by=c("qname", "ReadLen")]
    sam <- fread_sam(paste0(fileprefix, "_", segment, ".sam")) %>% process_sam(min_gap=mingap)
    rolled <- roll_alignments(sam)
    hits <- get_best_hits(rolled)
    hits <- merge(hits, readlens[, .(qname, ReadLen)], by="qname")
    hits[, ReadCov := round(AS/ReadLen*100, 1)]
    hits[, allele := tstrsplit(qname, "|", fixed=TRUE)[[2]]]

    result <- hits
    result[, fullgene := tstrsplit(qname, "|", fixed=TRUE)[2]]
    result[, c("gene", "allele") := tstrsplit(fullgene, "*", fixed=TRUE)]
    result[, ':='(prev_pos=NULL, new_group=NULL, group_end=NULL)]
    result[, functionality := tstrsplit(qname, "|", fixed=TRUE)[3]]
    return(result[, .(qname, rname, flag, anchor_pos, entry_pos, diff, cigar, AS, PI, ReadLen, ReadCov, gene, allele, redundant_hit, fullgene, functionality)])
}

extract_positions <- function(hits) {
    perfect_matches <- hits[ReadCov == 100 & PI == 100]
    perfect_bed <- perfect_matches[functionality == "F", .(rname, start=entry_pos-1, end=entry_pos+ReadLen-1, name=fullgene, score=flag)]
    print(paste0("Perfect matches: ", nrow(perfect_matches)))

    potentially_novel <- hits[functionality == "F"][ReadCov > 91][ReadCov < 100 | PI < 100][is.na(redundant_hit)]
    if(nrow(potentially_novel) > 0) {
         novel_bed <- potentially_novel[, .(rname, start=entry_pos-1, end=entry_pos+ReadLen-1, name=paste0(fullgene, "-like"), score=flag)]
         print(paste0("Potentially novel alleles: ", nrow(potentially_novel)))
         bed <- rbind(perfect_bed, novel_bed)[order(rname, start, end, name)]
    } else {
        bed <- perfect_bed
    }

    bed[score == 16, strand:="-"]
    bed[score == 0, strand:="+"]

    return(bed)
}

write_bed <- function(bed, segment) {
    filepath = paste0(argv$output_prefix, "_", segment, ".bed")
    fwrite(bed, file=filepath, sep="\t", col.names=FALSE)
}

v_hits <- process_annotations(argv$input_prefix, "V")
j_hits <- process_annotations(argv$input_prefix, "J")
c_hits <- process_annotations(argv$input_prefix, "C")

extract_positions(v_hits) %>% write_bed("V")
extract_positions(c_hits) %>% write_bed("C")
extract_positions(j_hits) %>% write_bed("J")

if(argv$locus == "IGH") {
    d_hits <- process_annotations(argv$input_prefix, "D")
    extract_positions(d_hits) %>% write_bed("D")

    annotations <- rbind(
        v_hits[functionality == "F"][, ':='(region = "V", locus=argv$locus)],
        d_hits[][, ':='(region = "D", locus=argv$locus)], #? Some D genes are identical in sequence yet differ in functionality - keep all to be safe
        j_hits[functionality == "F"][, ':='(region = "J", locus=argv$locus)],
        c_hits[, ':='(region = "C", locus=argv$locus)]
    )
} else {
    annotations <- rbind(
        v_hits[functionality == "F"][, ':='(region = "V", locus=argv$locus)],
        j_hits[functionality == "F"][, ':='(region = "J", locus=argv$locus)],
        c_hits[, ':='(region = "C", locus=argv$locus)]
    )
}

outfile=paste0(argv$output_prefix, ".best_hits.tsv")
filtered_annotations <- annotations[ReadCov > 91][order(rname, entry_pos)]

fwrite(filtered_annotations, file=outfile, col.names=TRUE, sep="\t")
