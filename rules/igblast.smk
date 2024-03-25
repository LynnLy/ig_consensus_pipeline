# https://www.imgt.org/vquest/refseqh.html#V-D-J-C-sets


rule download_igblast:
    # Downloads version 1.22.0 for linux. If using a different system or version, download and
    # extract it yourself and specify the bin/ncbi-igblast-1.22.0/ folder equivalent in `config/config.yml`
    output:
        "bin/ncbi-igblast-1.22.0/bin/igblastn",
        "bin/ncbi-igblast-1.22.0/bin/edit_imgt_file.pl",
        "bin/ncbi-igblast-1.22.0/bin/makeblastdb",
        "bin/ncbi-igblast-1.22.0/optional_file/human_gl.aux",
        directory("bin/ncbi-igblast-1.22.0/internal_data"),
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz
        tar -C bin/ -xvf ncbi-igblast-1.22.0-x64-linux.tar.gz
        """


# rule symlink_internal_db_to_working_directory:
#     # IgBlast does not detect the IGDATA environment variable and expects
#     # the internal data folder to be in the working directory
#     # This doesn't work and makes snakemake think the data's been updated constantly
#     input:
#         config['igblast_folder'] + "internal_data"
#     output:
#         internal=directory("internal_data")
#     shell:
#         """
#         ln -s {input} internal_data
#         """


rule download_imgt_databases:
    input:
        edit_imgt_file=config["igblast_folder"] + "bin/edit_imgt_file.pl",
        makeblastdb=config["igblast_folder"] + "bin/makeblastdb",
    output:
        formatted_V="databases/imgt/IG_V",
        formatted_D="databases/imgt/IG_D",
        formatted_J="databases/imgt/IG_J",
        indexed_V="databases/imgt/IG_V.ndb",
        indexed_D="databases/imgt/IG_D.ndb",
        indexed_J="databases/imgt/IG_J.ndb",
        C="databases/imgt/ncbi_human_c_genes.ndb",
    shell:
        """
        wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar
        wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta
        wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta
        wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHJ.fasta
        wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKV.fasta
        wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKJ.fasta
        wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLV.fasta
        wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLJ.fasta

        cat IGHV.fasta IGKV.fasta IGLV.fasta > IG_V.fasta
        cat IGHJ.fasta IGKJ.fasta IGLJ.fasta > IG_J.fasta

        {input.edit_imgt_file} IG_V.fasta > {output.formatted_V}
        {input.edit_imgt_file} IGHD.fasta > {output.formatted_D}
        {input.edit_imgt_file} IG_J.fasta > {output.formatted_J}

        {input.makeblastdb} -parse_seqids -dbtype nucl -in {output.formatted_V}
        {input.makeblastdb} -parse_seqids -dbtype nucl -in {output.formatted_D}
        {input.makeblastdb} -parse_seqids -dbtype nucl -in {output.formatted_J}

        tar -C databases/imgt/ -xf ncbi_human_c_genes.tar
        mv ncbi_human_c_genes.tar databases/imgt/
        mv IGHV.fasta IGHD.fasta IGHJ.fasta IGKV.fasta IGKJ.fasta IGLV.fasta IGLJ.fasta databases/imgt/
        rm IG_V.fasta IG_J.fasta
        """


rule get_C_fasta_from_blastdb:
    input:
        C="databases/imgt/ncbi_human_c_genes.ndb",
    output:
        "databases/imgt/ncbi_human_c_genes",
    conda:
        "../envs/blast.yml"
    shell:
        """
        blastdbcmd -entry all -db databases/imgt/ncbi_human_c_genes -out {output}
        """


rule download_ogrdb_databases:
    output:
        "databases/ogrdb/ogrdb_human_igh.V",
        "databases/ogrdb/ogrdb_human_igh.D",
        "databases/ogrdb/ogrdb_human_igh.J",
    shell:
        """
        wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/airr/ogrdb_human_blastdb.tar.gz
        tar -C databases/ogrdb -xvf ogrdb_human_blastdb.tar.gz
        """


rule rmdup_database:
    # Remove duplicates by sequence from the IMGT database
    input:
        "databases/imgt/IG_{segment}",
    output:
        "databases/imgt/IG_{segment}.rmdup",
    conda:
        "../envs/seqkit.yml"
    threads: 25
    shell:
        """
        cat {input} | seqkit rmdup --by-seq > {output}
        """


rule igblast_imgt:
    input:
        reads="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.fa",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="databases/imgt/IG_V",
        D="databases/imgt/IG_D",
        J="databases/imgt/IG_J",
        C="databases/imgt/ncbi_human_c_genes.ndb",
    params:
        aux=config["igblast_folder"] + "optional_file/human_gl.aux",
        C="databases/imgt/ncbi_human_c_genes",
        internal=config["igblast_folder"] + "internal_data",
    wildcard_constraints:
        database="IMGT_full",
    threads: 35
    resources:
        mem_mb=150000,
        mem_mib=150000,
    output:
        "results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.{database}.tsv",
    shell:
        """
        {input.igblastn} \
            -germline_db_V {input.V} \
            -germline_db_D {input.D} \
            -germline_db_J {input.J} \
            -c_region_db {params.C} \
            -organism human \
            -query {input.reads} \
            -auxiliary_data {params.aux} \
            -outfmt 19 \
            -num_threads {threads} \
            > {output}
        """


use rule igblast_imgt as igblast_imgt_dedup with:
    input:
        reads="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.fa",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="databases/imgt/IG_V.rmdup",
        D="databases/imgt/IG_D.rmdup",
        J="databases/imgt/IG_J.rmdup",
        Vi="databases/imgt/IG_V.rmdup.ndb",
        Di="databases/imgt/IG_D.rmdup.ndb",
        Ji="databases/imgt/IG_J.rmdup.ndb",
    wildcard_constraints:
        database="IMGT",


use rule igblast_imgt as igblast_curated with:
    input:
        reads="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.fa",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="databases/ogrdb/ogrdb_human_ig.V",
        D="databases/ogrdb/ogrdb_human_igh.D",
        J="databases/ogrdb/ogrdb_human_ig.J",
    wildcard_constraints:
        database="curated",


use rule igblast_imgt as igblast_personalized with:
    input:
        reads="results/ig_consensus/igblast_raw/chunked/{cDNA_sample}/{cDNA_sample}.{chunk}.fa",
        igblastn=config["igblast_folder"] + "bin/igblastn",
        V="results/manual/personalized_genome/20231025_wgs_duplex.all_V_segments.fasta",
        D="results/manual/personalized_genome/20231025_wgs_duplex.all_D_segments.fasta",
        J="results/manual/personalized_genome/20231025_wgs_duplex.all_J_segments.fasta",
        Vi="results/manual/personalized_genome/20231025_wgs_duplex.all_V_segments.fasta.ndb",
        Di="results/manual/personalized_genome/20231025_wgs_duplex.all_D_segments.fasta.ndb",
        Ji="results/manual/personalized_genome/20231025_wgs_duplex.all_J_segments.fasta.ndb",
        Ci="results/manual/personalized_genome/20231025_wgs_duplex.all_C_segments.fasta.ndb",
    wildcard_constraints:
        database="personalized",
