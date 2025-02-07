def get_ubam(wildcards):
    return gDNA_basecall_df.query("gDNA_sample == @wildcards.gDNA_sample")["ubam"]


rule map_to_diploid_assembly:
    input:
        reads = get_ubam,
        IGH = config["contigs_to_annotate_IGH"],
        IGL = config["contigs_to_annotate_IGL"]
    output:
        final_contigs = "results/ig_contigs/{gDNA_sample}.final_contigs.fa",
        bam = "results/ig_contigs/{gDNA_sample}.CpG.bam",
        bai = "results/ig_contigs/{gDNA_sample}.CpG.bam.bai"
    threads: 30
    conda: "../envs/minimap2.yml"
    shell:
        """
        cat {input.IGH} {input.IGL} > {output.final_contigs}
        samtools fastq {input.reads} -T '*' \
            | minimap2 -y -ax map-ont -t {threads} {output.final_contigs} - \
            | samtools sort -@ 4 -O bam -o {output.bam}
        samtools index {output.bam}
        """
        

rule filter_alignments:
    # Get alignments that are >80% coverage and < 5% SNPs
    input:
        bam = "results/ig_contigs/{gDNA_sample}.CpG.bam"
    output:
        snp_filt = "results/ig_contigs/{gDNA_sample}.filtered.{cov}cov.{snp}snp.bam"
    params:
        snp=lambda w: int(w.snp) / 100,
        cov_filt = lambda w: f"results/ig_contigs/{w.gDNA_sample}.filtered.{w.cov}cov.bam",
    conda: "../envs/pysam.yml"
    shell:
        """
        python scripts/filter_by_read_coverage.py \
            --input {input.bam} \
            --output {params.cov_filt} \
            -c {wildcards.cov}
        samtools index {params.cov_filt}
        samtools view -e "[NM]/(qlen+sclen) <= {params.snp}" \
            -O bam -o {output.snp_filt} {params.cov_filt}
        samtools index {output.snp_filt}
        """


rule modkit_pileup:
    input:
        bam="results/ig_contigs/{gDNA_sample}.filtered.{cov}cov.{snp}snp.bam",
        ref="results/ig_contigs/{gDNA_sample}.final_contigs.fa"
    output:
        bed="results/ig_contigs/modkit/{gDNA_sample}.{cov}cov.{snp}snp.bed",
        summary="results/ig_contigs/modkit/{gDNA_sample}.{cov}cov.{snp}snp.summary"
    threads: 20
    shell:
        """
        bin/dist/modkit summary {input.bam} -t {threads} --tsv > {output.summary} 
        bin/dist/modkit pileup {input.bam} {output.bed} \
            --cpg \
            --ref {input.ref} \
            --ignore h \
            --combine-strands 
        """


rule entropy:
    input:
        bam="results/ig_contigs/{gDNA_sample}.filtered.{cov}cov.{snp}snp.bam",
        ref="results/ig_contigs/{gDNA_sample}.final_contigs.fa"
    output:
        entropy = "results/ig_contigs/modkit/{gDNA_sample}.{cov}cov.{snp}snp.entropy",
    threads: 15
    shell:
        """
        bin/dist/modkit entropy --in-bam {input.bam} --ref {input.ref} \
            --out-bed {output.entropy} --threads {threads} --cpg \
            --header
        """


# rule modkit_extract:
#     # Read-level mod info parsed into a tsv


rule bam_to_tsv:
    input:
        bam = "results/ig_contigs/{gDNA_sample}.CpG.bam"
    output:
        "results/ig_contigs/{gDNA_sample}.CpG.tsv"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit bam {input.bam} 2>{output}
        """
    

# rule map_to_haploid_assembly:
#     input:
#         reads = get_ubam,
#         IGH = "results/manual/curation/IGH_final/hap1.contig_17_IGH.fa",
#         IGL = "results/manual/curation/IMGT_IGL/hapdup_phased_1.fasta"
#     output:
#         hap1_contigs = "results/ig_contigs/{gDNA_sample}.hap1_contigs.fa",
#         bam = "results/ig_contigs/{gDNA_sample}.hap1.CpG.bam",
#         bai = "results/ig_contigs/{gDNA_sample}.hap1.CpG.bam.bai"
#     threads: 50
#     conda: "../envs/minimap2.yml"
#     shell:
#         """
#         cat {input.IGH} {input.IGL} > {output.final_contigs}
#         samtools fastq {input.reads} -T '*' \
#             | minimap2 -y -ax map-ont -t {threads} {output.final_contigs} - \
#             | samtools sort -@ 4 -O bam -o {output.bam}
#         samtools index {output.bam}
#         """


# rule call_variants:
#     input:
#         bam = "results/ig_contigs/{gDNA_sample}.hap1.filtered.{cov}cov.{snp}snp.bam",
#     output:
#         vcf = "results/ig_contigs/{gDNA_sample}.hap1.filtered.{cov}cov.{snp}snp.vcf"
#     conda: "../envs/"
#     shell:
#         """
#         
#         """
#         
#         
# rule haplotag:
#     input:
#         bam = "results/ig_contigs/{gDNA_sample}.hap1.filtered.{cov}cov.{snp}snp.bam"
#         vcf =
#     output:
#         haplotagged = "results/ig_contigs/{gDNA_sample}.hap1.haplotagged.{cov}cov.{snp}snp.bam"
#     conda: "../envs/whatshap.yml"
#     shell:
#         """
#         whatshap haplotag
#         """