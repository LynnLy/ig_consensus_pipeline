### wf-sc-reliant method
def get_read_config(wildcards):
    read_config = cDNA_basecall_df.query("cDNA_sample == @wildcards.cDNA_sample")[
        "read_config"
    ]
    if len(read_config) > 1:
        raise ValueError(
            "cDNA_basecall.csv has multiple rows with the same cDNA_sample"
        )
    elif len(read_config) == 0:
        raise ValueError("Could not find cDNA_sample in cDNA_basecall.csv")
    return read_config[0]


def get_tagged_bam(wildcards):
    tagged_bam = cDNA_basecall_df.query("cDNA_sample == @wildcards.cDNA_sample")[
        "tagged_bam"
    ]
    if len(tagged_bam) > 1:
        raise ValueError(
            "cDNA_basecall.csv has multiple rows with the same cDNA_sample"
        )
    elif len(tagged_bam) == 0:
        raise ValueError("Could not find cDNA_sample in cDNA_basecall.csv")
    return tagged_bam[0]


rule get_full_length_readnames:
    # Input file is from the wf-single-cell work/temp files https://github.com/epi2me-labs/wf-single-cell/issues/68
    input:
        table=get_read_config,
    output:
        "results/wf-single-cell/{cDNA_sample}_output/full_length.readnames.txt",
    shell:
        """
        grep "full_len" {input.table} | tail -n +2 | cut -f 2 -d '\t' > {output}
        """


rule extract_cell_barcodes:
    input:
        bam=get_tagged_bam,
    output:
        cellbarcodes="results/wf-single-cell/{cDNA_sample}_output/{cDNA_sample}/bams/{cDNA_sample}.CBs.tsv",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools view {input.bam} | cut -f 1,12- \
             | grep "CB:Z:" \
             | awk '{{ printf "%s%s", $1, FS; for (i=2; i<=NF; i++) if ($i ~ /^CB/) printf "%s%s\\n", $i, (i==NF ? RS : FS) }}' \
             > {output.cellbarcodes}
        """


rule extract_trimmed_fastas:
    # Get a .fasta of all reads from wf-sc bam output
    input:
        bam=get_tagged_bam,
    output:
        fa="results/ig_consensus/igblast_raw/{cDNA_sample}.fa",
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        samtools bam2fq {input.bam} | seqkit fq2fa | seqkit rmdup -n > {output.fa}
        """


rule split_fastas:
    # Split fastas into 20 parts to decrease memory usage during IgBlast
    input:
        reads="results/ig_consensus/igblast_raw/{cDNA_sample}.fa",
    output:
        expand(
            "results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.fa",
            chunk=[f"part_{str(i).zfill(3)}" for i in range(1, 21)],
            allow_missing=True,
        ),
    params:
        outdir="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}",
    conda:
        "../envs/seqkit.yml"
    threads: 5
    shell:
        """
        seqkit split2 -p 20 -j {threads} -O {params.outdir} {input.reads}
        """


# Each part is then IgBlasted : See igblast.smk


rule annotate_with_cell_barcodes:
    # Add cell barcodes as a field
    input:
        igblast="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.tsv",
        cell_barcodes="results/wf-single-cell/{cDNA_sample}_output/{cDNA_sample}/bams/{cDNA_sample}.CBs.tsv",
    output:
        "results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.CBs.tsv",
    conda:
        "../envs/R.yml"
    threads: 25
    params:
        script="scripts/annotate_igblast_with_cell_barcodes.R",
    shell:
        """
        Rscript {params.script} --output {output} \
            --cell_barcodes {input.cell_barcodes} \
            --input {input.igblast}
        """


rule filter_igblast:
    input:
        igblast="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.CBs.tsv",
        c_lengths="databases/imgt/ncbi_human_c_genes.lengths",
        v_lengths="databases/imgt/IG_V.lengths",
    output:
        "results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.filtered.tsv",
    params:
        strictness="strict",
    threads: 20
    conda:
        "../envs/R.yml"
    shell:
        """
        Rscript scripts/filter_igblast_results.R \
            --input {input.igblast} --output {output} \
            --strictness {params.strictness} \
            --c_lengths {input.c_lengths} \
            --v_lengths {input.v_lengths}
        """


rule reduce_fields:
    # Limit to specific columns so we can load the whole file where necessary, such as for counting barcodes
    input:
        igblast="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.filtered.tsv",
    output:
        "results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.filtered.reduced.tsv",
    conda:
        "../envs/R.yml"
    threads: 20
    params:
        script="scripts/reduce_igblast_columns.R",
    shell:
        """
        Rscript {params.script} --output {output} \
            --input {input.igblast}
        """


use rule combine_files_with_header as combine_chunked_filtered_igblast with:
    input:
        expand(
            "results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.filtered.tsv",
            chunk=[f"part_{str(i).zfill(3)}" for i in range(1, 21)],
            allow_missing=True,
        ),
    output:
        "results/ig_consensus/igblast_raw/{cDNA_sample}.{database}.filtered.tsv",


use rule combine_files_with_header as combine_chunked_reduced_igblast with:
    input:
        expand(
            "results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.filtered.reduced.tsv",
            chunk=[f"part_{str(i).zfill(3)}" for i in range(1, 21)],
            allow_missing=True,
        ),
    output:
        "results/ig_consensus/igblast_raw/{cDNA_sample}.{database}.filtered.reduced.tsv",


rule get_ig_readnames:
    input:
        "results/ig_consensus/igblast_raw/{cDNA_sample}.{database}.filtered.tsv",
    output:
        "results/ig_consensus/igblast_raw/{cDNA_sample}.{roi_name}.{database}.filtered.readnames.txt",
    shell:
        """
        cut -f 1,4 {input} | tail -n +2 | grep {wildcards.roi_name} > {output}
        """


### Consensus generation with medaka smolecule
rule get_reads_and_cellbarcodes:
    # Script figures out which cell barcodes have enough IG reads
    # and how reads should be renamed for smolecule
    input:
        readnames="results/ig_consensus/igblast_raw/{cDNA_sample}.{roi_name}.{database}.filtered.readnames.txt",
        cell_barcodes="results/wf-single-cell/{cDNA_sample}_output/{cDNA_sample}/bams/{cDNA_sample}.CBs.tsv",
    output:
        guide="results/ig_consensus/medaka_input/{cDNA_sample}/{roi_name}.{database}.renaming_guide.tsv",
        depth="results/ig_consensus/medaka_input/{cDNA_sample}/{roi_name}.{database}.depths.tsv",
        readnames="results/ig_consensus/medaka_input/{cDNA_sample}/{roi_name}.{database}.smolecule_readnames.txt",
    params:
        script="scripts/filter_and_prep_for_smolecule.R",
        minreads=3,
        maxreads=config["medaka"]["max_reads"],
    threads: 40
    conda:
        "../envs/R.yml"
    shell:
        """
        Rscript {params.script} --cell_barcodes {input.cell_barcodes} \
            --readnames {input.readnames} --min_reads {params.minreads} \
            --max_reads {params.maxreads} --locus {wildcards.roi_name} \
            --output {output.guide} --output_depth {output.depth}
        cut -f 1 {output.guide} > {output.readnames}
        """


rule prep_for_smolecule:
    # Collect full length reads and rename them according to the renaming_guide
    # Further filter for full length adapter config from wf-single-cell
    input:
        guide="results/ig_consensus/medaka_input/{cDNA_sample}/{roi_name}.{database}.renaming_guide.tsv",
        readnames="results/ig_consensus/medaka_input/{cDNA_sample}/{roi_name}.{database}.smolecule_readnames.txt",
        fasta="results/ig_consensus/igblast_raw/{cDNA_sample}.fa",
        full_length_readnames="results/wf-single-cell/{cDNA_sample}_output/full_length.readnames.txt",
    output:
        filtered="results/ig_consensus/medaka_input/{cDNA_sample}_{roi_name}.{database}.candidates.fasta",
        renamed="results/ig_consensus/medaka_input/{cDNA_sample}_{roi_name}.{database}.pre_smolecule.fasta",
    conda:
        "../envs/seqkit.yml"
    threads: 40
    shell:
        """
        seqkit grep -f {input.readnames} {input.fasta} \
            | seqkit grep -f {input.full_length_readnames} \
            | seqkit fq2fa > {output.filtered}
        cat {output.filtered} | seqkit replace -p "(\S+.+?)$" -r '{{kv}}' -k {input.guide} \
            | seqkit sort --by-name \
            > {output.renamed}
        """


rule smolecule:
    input:
        fastq="results/ig_consensus/medaka_input/{cDNA_sample}_{roi_name}.{database}.pre_smolecule.fasta",
    output:
        "results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database}/consensus.fasta",
    params:
        outdir="results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database}/",
        depth=config["medaka"]["depth"],  # Minimum number of subreads in the group
        min_depth=config["medaka"]["min_depth"],  # Minimum depth at a particular base to be polished
    conda:
        "../envs/medaka.yml"
    threads: 60
    resources:
        mem_mb=120000,
        disk_mb=120000,
    shell:
        """
        medaka smolecule --threads 40 --batch_size 10 \
            --min_depth {params.min_depth} --depth {params.depth} \
            --model r941_e81_sup_g514 {params.outdir} \
            {input.fastq}
        """


use rule igblast_imgt as igblast_on_target_consensus with:
    input:
        reads="results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database1}/consensus.fasta",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="databases/imgt/IG_V.rmdup",
        D="databases/imgt/IG_D.rmdup",
        J="databases/imgt/IG_J.rmdup",
        Vi="databases/imgt/IG_V.rmdup.ndb",
        Di="databases/imgt/IG_D.rmdup.ndb",
        Ji="databases/imgt/IG_J.rmdup.ndb",
    output:
        "results/ig_consensus/igblast_consensus.{database1}/{cDNA_sample}_{roi_name}.{database2}.tsv",
    wildcard_constraints:
        database2="IMGT",


use rule igblast_on_target_consensus as igblast_on_target_consensus_non_dedup with:
    input:
        reads="results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database1}/consensus.fasta",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="databases/imgt/IG_V",
        D="databases/imgt/IG_D",
        J="databases/imgt/IG_J",
        Vi="databases/imgt/IG_V.ndb",
        Di="databases/imgt/IG_D.ndb",
        Ji="databases/imgt/IG_J.ndb",
    output:
        "results/ig_consensus/igblast_consensus.{database1}/{cDNA_sample}_{roi_name}.{database2}.tsv",
    wildcard_constraints:
        database2="IMGT_full",


use rule igblast_imgt as igblast_curated_consensus with:
    input:
        reads="results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database1}/consensus.fasta",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="databases/ogrdb/ogrdb_human_ig.V",
        D="databases/ogrdb/ogrdb_human_igh.D",
        J="databases/ogrdb/ogrdb_human_ig.J",
    output:
        "results/ig_consensus/igblast_consensus.{database1}/{cDNA_sample}_{roi_name}.{database2}.tsv",
    wildcard_constraints:
        database2="curated",


use rule igblast_imgt as igblast_personalized_consensus with:
    input:
        reads="results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database1}/consensus.fasta",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V=config['personalized_db'] + ".all_V_segments.fasta",
        D=config['personalized_db'] + ".all_D_segments.fasta",
        J=config['personalized_db'] + ".all_J_segments.fasta",
        Vi=config['personalized_db'] + ".all_V_segments.fasta.ndb",
        Di=config['personalized_db'] + ".all_D_segments.fasta.ndb",
        Ji=config['personalized_db'] + ".all_J_segments.fasta.ndb",
    wildcard_constraints:
        database2="personalized",
    output:
        "results/ig_consensus/igblast_consensus.{database1}/{cDNA_sample}_{roi_name}.{database2}.tsv",


use rule combine_files_with_header as gather_igblast_by_cDNA_sample with:
    input:
        expand(
            "results/ig_consensus/igblast_consensus.{database1}/{cDNA_sample}_{roi_name}.{database2}.tsv",
            roi_name=["IGL", "IGH", "IGK"],
            allow_missing=True,
        ),
    output:
        "results/ig_consensus/igblast_consensus.{database1}/{cDNA_sample}_all.{database2}.tsv",
