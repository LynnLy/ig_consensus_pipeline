## Aims
- Identify near full length IG reads from single-cell transcriptome sequencing and create a consensus sequence for each locus for each cell barcode

## Analysis
##### Step 0. Run [wf-single-cell](https://github.com/epi2me-labs/wf-single-cell)
  - in order to get adapter trimmed, strand oriented reads.
  - These outputs from wf-single-cell will be required:
    1. `bams/{name}.tagged.sorted.bam`, which will be used to extract cell barcodes and a trimmed fasta
    2. `read_config.tsv` file from the work/temp directory, which will be used to identify full length reads (reads possessing the expected adapter structure).
      - Importantly, this is **not a default output from wf-single-cell**. Instructions from: https://github.com/epi2me-labs/wf-single-cell/issues/68
> You can locate this file by looking in <out_dir>/execution/trace.txt for the line containing pipeline:stranding:combine_adapter_tables which should look like this:
> - "21	8f/a6c024	484342	pipeline:stranding:combine_adapter_tables (1)	COMPLETED	0	2023-12-16 09:11:01.834	1.1s	28ms	93.5%	0	0	1 MB	553.2 K"
>
> The second column specifies the prefix to the relevant work directory. So in this case, you would look for the following file:  <work_dir>/8f/a6c024.../<alias>_read_config.tsv, which contains the adapter configurations for each read.


##### The [consensus generation](rules/consensus_generation.smk) module will:
   - Extract a trimmed fasta from the wf-single-cell output
   - Run [IgBlast](https://ncbi.github.io/igblast/) on the raw reads
   - Annotate the IgBlast output with each read's cell barcode
   - [Filter](scripts/filter_igblast_results.R) the IgBlast output to get likely full length, non-chimeric IG reads
   - Create a smaller IgBlast output with fewer columns, which can be loaded into memory where necessary
   - [Extract](scripts/filter_and_prep_for_smolecule.R) reads for cell barcodes that have a full-length adapter configuration and enough IG read depth, and then format the names of those reads for medaka smolecule
   - Run [medaka smolecule](https://github.com/nanoporetech/medaka/blob/master/medaka/smolecule.py) to create a consensus for each cell barcode for each locus (IGH, IGK, or IGL)
   - Run IgBlast on consensus sequences


##### The [assembly](rules/whole_genome_assembly.smk) module will:
   - Create a haploid genome assembly using [Flye](https://github.com/fenderglass/Flye)
   - Remap the raw reads to the assembly with minimap2
   - Run [HapDup](https://github.com/KolmogorovLab/hapdup) to create a diploid assembly
   - **This module requires additional prerequisites**.

##### The [annotation](rules/annotation.smk) module will:
   - Download and extract human IG alleles from the [IMGT/GENE-DB](https://www.imgt.org/download/GENE-DB/)
   - blastn IG reads against provided IG contigs
   - Identify the [best IG gene hit for each locus](scripts/get_best_ig_hits.R)


##### The [clonotyping](rules/immcantation.smk) / [Immcantation](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html) module will:
   - Run IgBlast on the consensus sequences and parse the results with [changeo](https://changeo.readthedocs.io/en/stable/)
   - Use [SHazaM](https://shazam.readthedocs.io/en/stable/) and [SCOPer](https://scoper.readthedocs.io/en/stable/) to identify thresholds for clonal clustering and to annotate clones using either only heavy chains or adding light chains where available
   - Use changeo to reconstruct germline sequences, masking the D region
   - **This module requires additional prerequisites**.


## Installation
This has only been tested on Linux.

Prerequisites:
1. [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Download the [IgBlast](https://ncbi.github.io/igblast/cook/How-to-set-up.html) binary for your system and extract it.
   1. Specify the location of the IgBlast folder in `config/config.yml`
   2. Make a symlink to IgBlast's `internal_data` folder in the working directory. ie. `ln -s bin/ncbi-igblast-1.22.0/internal_data internal_data`


Create the environment from which to call Snakemake
```bash
git clone https://github.com/LynnLy/ig_consensus_pipeline
cd ig_consensus_pipeline
conda env create -f environment.yml
conda activate ig_consensus_pipeline
```


##### Test with pre-configured test data
```bash
snakemake --use-conda -p -j {threads} consensus annotate_igh
```


## Usage
### Configuration
- Parameters in `config/config.yml` should verified and updated as appropriate to the dataset.


##### Run the consensus generation module
- Fill in `config/cDNA_basecall.csv` with the sample name and a link to the two wf-single-cell outputs described in [Analysis](#step-0-run-wf-single-cellhttpsgithubcomepi2me-labswf-single-cell)
- Output: .fasta of consensus sequences for each locus (IGH/IGK/IGL) and cell barcode `results/ig_consensus/medaka_smolecule/{cDNA_sample}_{locus}.{database}/consensus.fasta`
```bash
snakemake --use-conda -p -j {threads} consensus
```


##### Run the assembly module
- Fill in `config/gDNA_basecall.csv` with the sample name and link to the genomic fastq
- Create a conda environment named "hapdup" and [install hapdup from source as described here](https://github.com/KolmogorovLab/hapdup?tab=readme-ov-file#source-installation)
- Output: `results/global_assembly/{sample}/hapdup/`
```bash
snakemake --use-conda -p -j {threads} assembly
```

Optional: Extract 10kb+ contigs that mapped to an IG region in GRCh38 with >50%+ ReadCov
- Requires a GRCh38 reference genome to be placed at `results/refgenome/GRCh38.fa`
- Manual inspection recommended to make sure the right contigs were captured
- Output: `results/ig_contigs/{sample}/{locus}/hapdup_phased_both.fasta`
```bash
snakemake --use-conda -p -j {threads} ig_contigs
```


##### Run the annotation module
- Fill in `config/config.yml` with the sample name under "sample_to_annotate"
- Fill in `config/config.yml` with a path to IG contigs under the "contigs_to_annotate/IGK/IGL" fields
- Outputs:
  - `results/annotations/{sample}.{locus}.best_hits.tsv`: List of all the best alignments of functional genes
  - `results/annotations/{sample}.{locus}_{VDJ_segment}.fasta`
- Manual inspection recommended: IGH D segments may have spurious alignments outside of the expected IGH D region because they are so short.
```bash
snakemake --use-conda -p -j {threads} annotate_igh annotate_igl
```


##### Run the clonotyping / immcantation module
- Uses the consensus sequences generated above
- Enter the immcantation docker image, and copy the directory
`/usr/local/share/germlines/imgt/human` from the docker image to
`results/refgenome/immcantation_gapped/human` on the local machine (If someone knows the original link to the gapped genes database, please let me know!)
- Additionally, create a conda environment named "immcantation" and install immcantation packages through R
- Output `results/ig_consensus/clonotyping/{cDNA_sample}_all.{database}.fmt7_db-pass_clone-pass_germ-pass-*.tsv`
```bash
conda create --name immcantation -c conda-forge -c bioconda r-essentials r-tidyverse r-data.table r-patchwork r-argparser
conda activate immcantation
R
install.packages(c("dowser", "scoper", "shazam", "alakazam"))
```

```bash
snakemake --use-conda -p -j {threads} clonotyping
```
