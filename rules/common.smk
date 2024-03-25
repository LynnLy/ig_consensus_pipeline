outdir = "results"


def activate_conda_env(name: str) -> str:
    """Manually activate a conda environment"""
    return f"set +u; source ~/miniconda3/etc/profile.d/conda.sh ; conda activate ; conda activate {name}; "


def to_log(path: str) -> str:
    """Log file location based on output file"""
    return str(outdir / "logs" / path) + ".log"


def to_benchmark(path: str) -> str:
    """Benchmark file location based on output file"""
    return str(outdir / "benchmarks" / path) + ".bench.txt"


rule combine_files_with_header:
    input:
        expand(
            "results/example/{roi_name}.tsv",
            roi_name=["IGL", "IGH", "IGK"],
            allow_missing=True,
        ),
    output:
        "results/example/all.tsv",
    shell:
        """
        single_header_cat () {{
            output=$1
            shift
            head -n 1 $1 > $output
            for FILE in "$@"
            do
                tail -n +2 $FILE >> $output
            done
        }}
        single_header_cat {output} {input}
        """


rule makeblastdb:
    input:
        db="{anything}",
    output:
        index="{anything}.ndb",
    conda:
        "../envs/blast.yml"
    shell:
        """
        makeblastdb -in {input.db} -parse_seqids -dbtype nucl
        """


rule get_fasta_lengths:
    input:
        "{anything}",
    output:
        "{anything}.lengths",
    conda:
        "../envs/bioawk.yml"
    shell:
        """
        bioawk -c fastx '{{ print $name, length($seq) }}' {input} > {output}
        """
