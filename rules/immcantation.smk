rule changeo_assign_genes:
    # AssignGenes.py is a wrapper around IgBLAST
    input:
        reads="results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database1}/consensus.fasta",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="databases/imgt/IG_V.rmdup",
        D="databases/imgt/IG_D.rmdup",
        J="databases/imgt/IG_J.rmdup",
        Vi="databases/imgt/IG_V.rmdup.ndb",
        Di="databases/imgt/IG_D.rmdup.ndb",
        Ji="databases/imgt/IG_J.rmdup.ndb",
        C="databases/imgt/ncbi_human_c_genes.ndb",
    output:
        fmt7="results/ig_consensus/clonotyping/{cDNA_sample}_{roi_name}.{database1}.{database2}.fmt7.tsv",
    params:
        aux=config["igblast_folder"] + "optional_file/human_gl.aux",
        C="databases/imgt/ncbi_human_c_genes",
        internal=config["igblast_folder"] + "internal_data",
    wildcard_constraints:
        database2="IMGT",
    threads: 10
    shell:
        """
        {input.igblastn} \
            -germline_db_V {input.V} \
            -germline_db_D {input.D} \
            -germline_db_J {input.J} \
            -c_region_db {params.C} \
            -auxiliary_data {params.aux} \
            -domain_system imgt -ig_seqtype Ig -organism human \
            -outfmt '7 std qseq sseq btop' \
            -query {input.reads} \
            -out {output.fmt7}
        """


use rule changeo_assign_genes as changeo_assign_genes_personalized with:
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


rule changeo_makedb:
    input:
        fmt7="results/ig_consensus/clonotyping/{cDNA_sample}_{roi_name}.{database1}.{database2}.fmt7.tsv",
        reads="results/ig_consensus/medaka_smolecule/{cDNA_sample}_{roi_name}.{database1}/consensus.fasta",
        gapped="results/refgenome/immcantation_gapped/human/vdj/",
    output:
        db="results/ig_consensus/clonotyping/{cDNA_sample}_{roi_name}.{database1}.{database2}.fmt7_db-pass.tsv",
    conda:
        "../envs/changeo.yml"
    shell:
        """
        MakeDb.py igblast -i {input.fmt7} -s {input.reads} \
            -r {input.gapped}  \
            --extended
        """


use rule combine_files_with_header as gather_changeo with:
    input:
        expand(
            "results/ig_consensus/clonotyping/{cDNA_sample}_{roi_name}.{database1}.{database2}.fmt7_db-pass.tsv",
            roi_name=["IGL", "IGH", "IGK"],
            allow_missing=True,
        ),
    output:
        "results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass.tsv",


rule defineClones:
    input:
        db="results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass.tsv",
    output:
        plots="results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.clones.png",
        clones="results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass_clone-pass.tsv",
    params:
        conda=activate_conda_env("immcantation"),
        script="scripts/define_clones.R",
    shell:
        """
        {params.conda}
        Rscript {params.script} \
            --input {input.db} \
            --sample_id {wildcards.cDNA_sample} \
            --output {output.clones} \
            --plots {output.plots}
        """


rule createGermlines_only_heavy:
    input:
        clones="results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass_clone-pass.tsv",
        v="results/refgenome/immcantation_gapped/human/vdj/imgt_human_IGHV.fasta",
        d="results/refgenome/immcantation_gapped/human/vdj/imgt_human_IGHD.fasta",
        j="results/refgenome/immcantation_gapped/human/vdj/imgt_human_IGHJ.fasta",
    output:
        "results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass_clone-pass_germ-pass.tsv",
    conda:
        "../envs/changeo.yml"
    shell:
        """
        CreateGermlines.py -d {input.clones} -g dmask \
            --cloned \
            -r {input.v} {input.d} {input.j} \
            -o {output}
        """


rule createGermlines_heavy_and_light:
    input:
        clones="results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass_clone-pass.tsv",
        #gapped="results/refgenome/immcantation_gapped/human/vdj/"
        v="results/refgenome/immcantation_gapped/human/vdj/imgt_human_IGHV.fasta",
        d="results/refgenome/immcantation_gapped/human/vdj/imgt_human_IGHD.fasta",
        j="results/refgenome/immcantation_gapped/human/vdj/imgt_human_IGHJ.fasta",
    output:
        "results/ig_consensus/clonotyping/{cDNA_sample}_all.{database1}.{database2}.fmt7_db-pass_clone-pass_germ-pass-hl.tsv",
    conda:
        "../envs/changeo.yml"
    shell:
        """
        CreateGermlines.py -d {input.clones} -g dmask \
            --cloned --cf clone_subgroup_id \
            -r {input.v} {input.d} {input.j} \
            -o {output}
        """
