# Part 1. Generate the assembly
# Part 2. Identify IG contigs by comparing to GRCh38
def get_fastq(wildcards):
    return gDNA_basecall_df.query("gDNA_sample == @wildcards.sample")["fastq"]


rule flye_haploid:
    input:
        fastq=get_fastq,
    output:
        asm="results/global_assembly/{sample}/assembly.fasta",
    params:
        prefix="results/global_assembly/{sample}",
        mode=config["flye"]["mode"],
    threads: 50
    conda:
        "../envs/flye.yml"
    shell:
        """
        flye {params.mode} {input.fastq} \
            -t {threads} -o {params.prefix}
        """


rule remap_to_flye:
    input:
        fastq=get_fastq,
        asm="results/global_assembly/{sample}/assembly.fasta",
    output:
        bam="results/global_assembly/{sample}/remapped_to_asm.bam",
        bai="results/global_assembly/{sample}/remapped_to_asm.bam.bai",
    conda:
        "../envs/minimap2.yml"
    threads: 50
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.asm} {input.fastq} \
            | samtools sort -O bam -o {output.bam}
        samtools index {output.bam}
        """


rule do_hapdup:
    input:
        asm="results/global_assembly/{sample}/assembly.fasta",
        mapping="results/global_assembly/{sample}/remapped_to_asm.bam",
    output:
        "results/global_assembly/{sample}/hapdup/hapdup_phased_1.fasta",
        "results/global_assembly/{sample}/hapdup/hapdup_phased_2.fasta",
        "results/global_assembly/{sample}/hapdup/hapdup_dual_1.fasta",
        "results/global_assembly/{sample}/hapdup/hapdup_dual_2.fasta",
        "results/global_assembly/{sample}/hapdup/margin/MARGIN_PHASED.haplotagged.bam",
    params:
        conda=activate_conda_env("hapdup"),
        prefix="results/global_assembly/{sample}/hapdup",
    threads: 50
    shell:
        """
        {params.conda}
        ./bin/hapdup/hapdup.py \
            --assembly {input.asm} \
            --bam {input.mapping} \
            --out-dir {params.prefix} \
            --rtype ont \
            --min-aligned-length 5000 \
            -t {threads}
        """


rule combined_haplotyped_fastas_global:
    input:
        hap1="results/global_assembly/{sample}/hapdup/hapdup_{mode}_1.fasta",
        hap2="results/global_assembly/{sample}/hapdup/hapdup_{mode}_2.fasta",
    output:
        both="results/global_assembly/{sample}/hapdup/hapdup_{mode}_both.fasta",
    wildcard_constraints:
        mode="(dual)|(phased)",
    shell:
        """
        sed 's/^>/>hap1./' {input.hap1} > {output.both}
        sed 's/^>/>hap2./' {input.hap2} >> {output.both}
        """


# Part 2. Identify IG contigs by comparing to GRCh38
rule map_assembly_to_GRCh38:
    input:
        asm="results/global_assembly/{sample}/hapdup/hapdup_phased_both.fasta",
        ref="results/refgenome/GRCh38.fa",
    output:
        bam="results/global_assembly/{sample}/hapdup/hapdup_phased_both.GRCh38.bam",
    threads: 10
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -t {threads} -ax asm10 {input.ref} {input.asm} \
            | samtools sort -O bam -o {output.bam}
        samtools index {output.bam}
        """


rule extract_IG_contigs:
    input:
        asm="results/global_assembly/{sample}/hapdup/hapdup_phased_both.fasta",
        bam="results/global_assembly/{sample}/hapdup/hapdup_phased_both.GRCh38.bam",
        refgenome="results/refgenome/GRCh38.fa",
    output:
        regional_asm="results/ig_contigs/{sample}/{roi_name}/hapdup_phased_both.fasta",
        regional_bam="results/ig_contigs/{sample}/{roi_name}/hapdup_phased_both.bam",
        regional_tsv="results/ig_contigs/{sample}/{roi_name}/hapdup_phased_both.tsv",
        contigs="results/ig_contigs/{sample}/{roi_name}/hapdup_phased_both.contigs.txt",
        paf="results/ig_contigs/{sample}/{roi_name}/hapdup_phased_both.paf",
    params:
        chrom=lambda w: roi_df.query(
            "roi_name == @w.roi_name and refgenome == 'GRCh38'"
        )["chrom"][0],
        start=lambda w: roi_df.query(
            "roi_name == @w.roi_name and refgenome == 'GRCh38'"
        )["start"][0],
        end=lambda w: roi_df.query(
            "roi_name == @w.roi_name and refgenome == 'GRCh38'"
        )["end"][0],
    conda:
        "../envs/samtools.yml"
    threads: 10
    shell:
        """
        samtools view -h -O bam {input.bam} \
            {params.chrom}:{params.start}-{params.end} \
            > {output.regional_bam} ;
        samtools index {output.regional_bam} ;
        seqkit bam {output.regional_bam} 2> {output.regional_tsv} ;
        cat {output.regional_tsv} \
            | awk '$7 >= 10000' \
            | awk '$12 >= 50' \
            | grep -v Read \
            | cut -f 1 | sort | uniq > {output.contigs}
        seqkit grep -f {output.contigs} {input.asm} > {output.regional_asm}
        minimap2 -t {threads} -x asm10 {input.refgenome} {output.regional_asm} > {output.paf}
        """
