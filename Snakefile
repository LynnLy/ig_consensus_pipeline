import pandas as pd
import yaml
from pathlib import Path
import os
import re

from snakemake.utils import validate, min_version

min_version("8.0.0")

CWD = Path(os.getcwd())


configfile: "config/config.yml"


cDNA_basecall_df = pd.read_table(
    "config/cDNA_basecall.csv", comment="#", sep=","
).set_index(["cDNA_sample"], drop=False)

gDNA_basecall_df = pd.read_table(
    "config/gDNA_basecall.csv", comment="#", sep=","
).set_index(["gDNA_sample"], drop=False)

roi_df = pd.read_table("config/roi_coordinates.csv", comment="#", sep=",").set_index(
    ["roi_name"], drop=False
)

cDNA_samples = cDNA_basecall_df["cDNA_sample"]
gDNA_samples = gDNA_basecall_df["gDNA_sample"]
assemblers = config["assemblers"]
haplotypes = ["1", "2"]
segment_types = config["segment_types"]
chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


include: "rules/common.smk"  # python helper functions
include: "rules/igblast.smk"
include: "rules/consensus_generation.smk"
include: "rules/immcantation.smk"
include: "rules/whole_genome_assembly.smk"
include: "rules/annotation.smk"
include: "rules/methylation.smk"

wildcard_constraints:
    #gDNA_sample='|'.join([re.escape(x) for x in basecall_df["run_id"]]),
    chrom=r"(chr\d+)|(chrX)|(chrY)",
    roi_name="|".join([re.escape(x) for x in ["IGH", "IGK", "IGL"]]),
    haplotype="(1)|(2)|(none)",
    segment_type="|".join([re.escape(x) for x in segment_types]),
    chunk=r"part_\d+",
    cDNA_sample="|".join([re.escape(x) for x in cDNA_samples]),
    database="|".join(["curated", "IMGT", "personalized", "IMGT_full"]),


rule consensus:
    input:
        medaka_input=expand(
            "results/ig_consensus/medaka_input/{cDNA_sample}_{roi_name}.{database}.pre_smolecule.fasta",
            cDNA_sample=cDNA_samples,
            roi_name=["IGH", "IGL", "IGK"],
            database="IMGT"
        ),
        medaka=expand(
            "results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database}/consensus.fasta",
            cDNA_sample=cDNA_samples,
            roi_name=["IGH", "IGL", "IGK"],
            database="IMGT"
        ),
        igblast_consensus=expand(
            "results/ig_consensus/igblast_consensus.{database1}/{cDNA_sample}_{roi_name}.{database2}.tsv",
            cDNA_sample=cDNA_samples,
            roi_name=["IGH", "IGL", "IGK"],
            database1="IMGT",
            database2=config['database']
        )


rule clonotyping:
    input:
        germlines_heavy=expand(
            "results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass_clone-pass_germ-pass.tsv",
            cDNA_sample=cDNA_samples,
            database1="IMGT",
            database2=config["database"],
        ),
        germlines_heavy_and_light=expand(
            "results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass_clone-pass_germ-pass-hl.tsv",
            cDNA_sample=cDNA_samples,
            database1="IMGT",
            database2=config["database"],
        ),


rule assembly:
    input:
        assembled=expand(
            "results/global_assembly/{gDNA_sample}/hapdup/hapdup_{mode}_both.fasta",
            gDNA_sample=gDNA_samples,
            mode=["phased", "dual"],
        ),


rule ig_contigs:
    input:
        regional_asm=expand(
            "results/ig_contigs/{sample}/{roi_name}/hapdup_phased_both.fasta",
            sample=gDNA_samples,
            roi_name=["IGH", "IGK", "IGL"],
        ),


rule annotate_igh:
    input:
        heavy_blast_results=expand(
            "results/annotations/{gDNA_sample}.{roi_name}_{segment_type}.fasta",
            gDNA_sample=config["sample_to_annotate"],
            segment_type=["V", "D", "J", "C"],
            roi_name="IGH",
        ),


rule annotate_igl:
    input:
        light_blast_results=expand(
            "results/annotations/{gDNA_sample}.{roi_name}_{segment_type}.fasta",
            gDNA_sample=config["sample_to_annotate"],
            segment_type=["V", "J", "C"],
            roi_name="IGL",
        ),


rule annotate_igk:
    input:
        light_blast_results=expand(
            "results/annotations/{gDNA_sample}.{roi_name}_{segment_type}.fasta",
            gDNA_sample=config["sample_to_annotate"],
            segment_type=["V", "J", "C"],
            roi_name="IGK",
        ),

rule methylation:
    input:
        mapped_to_diploid_assembly=expand(
            "results/ig_contigs/modkit/{gDNA_sample}.{cov}cov.{snp}snp.bed",
            gDNA_sample=config["sample_to_annotate"],
            cov=[80, 90],
            snp=3
        ),
        entropy=expand(
            "results/ig_contigs/modkit/{gDNA_sample}.{cov}cov.{snp}snp.entropy",
            gDNA_sample=config["sample_to_annotate"],
            cov=[80, 90],
            snp=3
        )