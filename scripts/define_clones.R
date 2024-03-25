library(argparser)

p <- arg_parser("Identify clone threshold and define clones")
p <- add_argument(p, "--output", help="table with clone_id column")
p <- add_argument(p, "--input", help="Output from changeo MakeDb")
p <- add_argument(p, "--plots", help="Output: plots")
p <- add_argument(p, "--sample_id", help="Label to append to the sequence_id")
p <- add_argument(p, "--threads", help="threads", default=8)
argv <- parse_args(p)

library(tidyverse)
library(data.table)
library(dowser)
library(scoper)
library(shazam)
library(alakazam)
library(patchwork)

exclude_cells_without_heavy_chains <- function(data) {
  heavy_cells <- data[locus == "IGH"]$cell_barcode
  light_cells <- data[locus %in% c("IGK", "IGL")]$cell_barcode
  no_heavy_cells <- setdiff(light_cells, heavy_cells)
  with_heavy <- data[!cell_barcode %in% no_heavy_cells]
  cat(paste(length(no_heavy_cells), "rows removed.\n"))
  cat(paste("There are", nrow(with_heavy), "rows in the data after filtering out cells without heavy chains.\n"))
  return(with_heavy)
}

data <- fread(argv$input)
data <- data[productive == "T"]
cat(paste(nrow(data), "productive sequences"))
data[, cell_barcode := tstrsplit(sequence_id, "_", fixed=TRUE)[[1]]]
data[, sample_id := argv$sample_id]
data[, sequence_id := paste(cell_barcode, locus, sample_id, sep=".")]
with_heavy <- exclude_cells_without_heavy_chains(data) %>% as.data.frame()

dist_nearest <- distToNearest(
  with_heavy,
  model="hh_s5f",
  nproc=10,
  cellIdColumn="cell_barcode",
  onlyHeavy=FALSE
)
output <- findThreshold(dist_nearest$dist_nearest, method="gmm", model="gamma-norm", cutoff="user", spc=0.995)
threshold <- output@threshold
threshold_plot <- plot(output)

sc_clones <- hierarchicalClones(
  with_heavy,
  cell_id = "cell_barcode",
  threshold=threshold,
  only_heavy=FALSE,
  split_light=TRUE,
  summarize_clones=FALSE,
  fields="sample_id",
  nproc=argv$threads
)

# Use both Heavy and Light chain

sc_hl_clones <- resolveLightChains(
  data=sc_clones,
  cell="cell_barcode",
  nproc=argv$threads
)

abund <- estimateAbundance(sc_clones %>% filter(locus == "IGH"), group="locus", nboot=100)
abundance_plot <- plot(abund)

clone_sizes <- countClones(sc_clones %>% filter(locus == "IGH"))
clone_sizes_plot <- clone_sizes %>% filter(seq_count > 1) %>%
  ggplot() +
  geom_histogram(aes(x = seq_count), binwidth=1) +
  ggtitle(paste0("Clone sizes (>1) ", argv$sample_id)) +
  theme_bw(14)

div <- alphaDiversity(sc_clones %>% filter(locus == "IGH"), group="locus", nboot=100)
div_plot <- plot(div)
all_plots <- threshold_plot / abundance_plot / clone_sizes_plot / div_plot + plot_layout(ncol=1)

ggsave(all_plots, file=argv$plots, height=15, width=6)
fwrite(sc_hl_clones, file=argv$output, sep="\t")
