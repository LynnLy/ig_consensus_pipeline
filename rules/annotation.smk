rule download_IG_segments:
    output:
        "results/annotations/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.fasta",
    shell:
        """
        curl https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP > {output}
        """


rule extract_IG_segments:
    input:
        "results/annotations/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.fasta",
    output:
        "results/annotations/{roi_name}_{segment_type}.fasta",
    conda:
        "../envs/seqkit.yml"
    wildcard_constraints:
        segment_type="V|J|D",
    shell:
        """
        cat {input} \
            | seqkit grep -n -r -p "\|Homo sapiens\|" \
            | seqkit grep -n -r -p {wildcards.roi_name} \
            | seqkit grep -n -r \
                -p "{wildcards.segment_type}-REGION" \
                -p "{wildcards.segment_type}-LIKE" \
            | seqkit replace -p '^([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|.+' \
                 -r "\$1|\$2|\$4|\$5|" > {output}
        """


use rule extract_IG_segments as extract_IGK_IGL_Constants with:
    wildcard_constraints:
        segment_type="C",
        roi_name="(IGL)|(IGK)",


rule extract_IGH_Constants:
    # Requires a slightly different search because each IGHC gene has
    # at least 3 chain segments and don't have the "-REGION" tags
    input:
        "results/annotations/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.fasta",
    output:
        "results/annotations/{roi_name}_{segment_type}.fasta",
    conda:
        "../envs/seqkit.yml"
    wildcard_constraints:
        segment_type="C",
        roi_name="IGH",
    shell:
        """
        cat {input} \
            | seqkit grep -n -r -p "\|Homo sapiens\|" \
            | seqkit grep -n -r -p {wildcards.roi_name} \
            | seqkit grep -n -r -v -p "REGION" -p "LIKE" \
            | seqkit replace -p '^([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|.+' \
                 -r "\$1|\$2|\$4|\$5|" \
            > {output}
        """


def get_assembled_IG_contigs(wildcards):
    if wildcards.roi_name == "IGH":
        return config["contigs_to_annotate_IGH"]
    if wildcards.roi_name == "IGK":
        return config["contigs_to_annotate_IGK"]
    if wildcards.roi_name == "IGL":
        return config["contigs_to_annotate_IGL"]


rule blastn:
    input:
        asm=get_assembled_IG_contigs,
        db="results/annotations/{roi_name}_{segment_type}.fasta",
        index="results/annotations/{roi_name}_{segment_type}.fasta.ndb",
    output:
        bam="results/annotations/{sample}.{roi_name}_{segment_type}.bam",
        sam="results/annotations/{sample}.{roi_name}_{segment_type}.sam",
    wildcard_constraints:
        segment_type="V|C",
    conda:
        "../envs/blast.yml"
    shell:
        """
        blastn -parse_deflines -query {input.asm} \
            -db {input.db} -outfmt "17 SQ" -penalty -5 \
            | samtools sort -O bam > {output.bam} ;
        samtools index {output.bam}
        samtools view -h {output.bam} > {output.sam}
        """


rule blastn_short:
    input:
        asm=get_assembled_IG_contigs,
        db="results/annotations/{roi_name}_{segment_type}.fasta",
        index="results/annotations/{roi_name}_{segment_type}.fasta.ndb",
    output:
        bam="results/annotations/{sample}.{roi_name}_{segment_type}.bam",
        sam="results/annotations/{sample}.{roi_name}_{segment_type}.sam",
    wildcard_constraints:
        segment_type="J|D",
    conda:
        "../envs/blast.yml"
    shell:
        """
        blastn -parse_deflines -query {input.asm} \
            -db {input.db} -outfmt "17 SQ" -penalty -4 \
            -task blastn-short \
            | samtools sort -O bam > {output.bam} ;
        samtools index {output.bam}
        samtools view -h {output.bam} > {output.sam}
        """


rule blastn_ncbi_c_genes:
    # Alternative C gene annotation where the C genes are not broken up
    input:
        asm=get_assembled_IG_contigs,
        db="results/transcriptome/ncbi-igblast-1.22.0/reference/ncbi_human_c_genes",
    output:
        bam="results/annotations/{sample}.{roi_name}.ncbi_{segment_type}.bam",
        sam="results/annotations/{sample}.{roi_name}.ncbi_{segment_type}.sam",
    wildcard_constraints:
        segment_type="C",
    conda:
        "../envs/blast.yml"
    shell:
        """
        blastn -parse_deflines -query {input.asm} \
            -db {input.db} -outfmt "17 SQ" -penalty -5 \
            | samtools sort -O bam > {output.bam} ;
        samtools index {output.bam}
        samtools view -h {output.bam} > {output.sam}
        """


rule seqkit_bam:
    input:
        bam="results/annotations/{sample}.{roi_name}_{segment_type}.bam",
    output:
        tsv="results/annotations/{sample}.{roi_name}_{segment_type}.tsv",
    conda:
        "../envs/seqkit.yml"
    wildcard_constraints:
        segment_type="V|D|J|C",
    shell:
        """
        seqkit bam {input.bam} 2> {output.tsv}
        """


rule get_best_light_hits:
    input:
        expand(
            "results/annotations/{sample}.{roi_name}_{segment_type}.{suffix}",
            suffix=["tsv", "sam"],
            segment_type=["V", "J", "C"],
            allow_missing=True,
        ),
    output:
        best_hits="results/annotations/{sample}.{roi_name}.best_hits.tsv",
        V="results/annotations/{sample}.{roi_name}_V.bed",
        J="results/annotations/{sample}.{roi_name}_J.bed",
        C="results/annotations/{sample}.{roi_name}_C.bed",
    params:
        prefix="results/annotations/{sample}.{roi_name}",
    wildcard_constraints:
        roi_name="(IGK)|(IGL)",
    conda:
        "../envs/R.yml"
    shell:
        """
        Rscript scripts/get_best_ig_hits.R \
            --output_prefix {params.prefix} \
            --input_prefix {params.prefix} \
            --locus {wildcards.roi_name}
        """


rule get_best_heavy_hits:
    input:
        expand(
            "results/annotations/{sample}.{roi_name}_{segment_type}.{suffix}",
            suffix=["tsv", "sam"],
            segment_type=["V", "D", "J", "C"],
            allow_missing=True,
        ),
    output:
        best_hits="results/annotations/{sample}.{roi_name}.best_hits.tsv",
        V="results/annotations/{sample}.{roi_name}_V.bed",
        D="results/annotations/{sample}.{roi_name}_D.bed",
        J="results/annotations/{sample}.{roi_name}_J.bed",
        C="results/annotations/{sample}.{roi_name}_C.bed",
    params:
        prefix="results/annotations/{sample}.{roi_name}",
    wildcard_constraints:
        roi_name="IGH",
    conda:
        "../envs/R.yml"
    shell:
        """
        Rscript scripts/get_best_ig_hits.R \
            --output_prefix {params.prefix} \
            --input_prefix {params.prefix} \
            --locus {wildcards.roi_name}
        """


rule extract_alleles:
    # Extract alleles and remove the (-) and (+) part of the names
    input:
        asm=get_assembled_IG_contigs,
        bed="results/annotations/{sample}.{roi_name}_{segment}.bed",
    output:
        "results/annotations/{sample}.{roi_name}_{segment}.fasta",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools getfasta -fi {input.asm} -bed {input.bed} -nameOnly -s |
            sed 's/(-)//' | sed 's/(+)//' > {output}
        """


# The remaining rules use annotations from only IGH and IGL, because we did not have a full IGK assembly.
# Instead, IGK genes are taken from the database
rule combine_VJ:
    input:
        IGH="results/annotations/{sample}.IGH_{segment}.fasta",
        IGL="results/annotations/{sample}.IGL_{segment}.fasta",
        IGK="databases/imgt/IGK{segment}.fasta",
    output:
        "results/annotations/{sample}.all_{segment}.fasta",
    wildcard_constraints:
        segment="V|J",
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        cat {input.IGH} {input.IGL} {input.IGK} | seqkit rmdup --by-seq > {output}
        """


rule combine_D:
    input:
        IGH="results/annotations/{sample}.IGH_{segment}.fasta",
    output:
        IGH="results/annotations/{sample}.all_{segment}.fasta",
    wildcard_constraints:
        segment="D",
    shell:
        """
        cat {input.IGH} | seqkit rmdup --by-seq > {output}
        """
