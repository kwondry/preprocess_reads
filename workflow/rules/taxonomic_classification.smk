# Taxonomic classification using Kraken2

rule kraken2:
    """Taxonomic classification of host-filtered reads using Kraken2."""
    input:
        r1 = "results/hisat2/{sample}_1.human_filtered.fastq.gz",
        r2 = "results/hisat2/{sample}_2.human_filtered.fastq.gz"
    output:
        report = "results/kraken/{sample}_kraken_report.txt"
    params:
        kraken_db = config["kraken_db_path"],
        confidence = config["kraken"]["confidence"]
    threads: config["kraken"]["threads"]
    resources:
        cpus_per_task=config["kraken"]["threads"],
        mem_mb=800000,
        runtime="4h",
        slurm_partition="highmem"
    conda:
        "../envs/kraken2.yaml"
    log:
        "results/logs/kraken/{sample}_kraken2.log"
    shell:
        """
        kraken2 --db {params.kraken_db} \\
            --paired \\
            --gzip-compressed \\
            --threads {threads} \\
            --confidence {params.confidence} \\
            --report {output.report} \\
            {input.r1} {input.r2} \\
            2> {log}
        """

rule bracken:
    """Abundance estimation from Kraken2 report using Bracken."""
    input:
        report = "results/kraken/{sample}_kraken_report.txt"
    output:
        bracken_report = "results/bracken/{sample}_bracken_report.txt",
        bracken_output = "results/bracken/{sample}_bracken_output.txt"
    params:
        kraken_db = config["kraken_db_path"],
        read_len = config["bracken"]["read_len"],  # e.g., 150
        level = config["bracken"]["level"]  # e.g., "S" for species
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=16000,
        runtime="12h",
        slurm_partition="short"
    conda:
        "../envs/kraken2.yaml"  # Bracken is typically included with kraken2 environment
    log:
        "results/logs/bracken/{sample}_bracken.log"
    shell:
        """
        bracken -d {params.kraken_db} \\
            -i {input.report} \\
            -o {output.bracken_output} \\
            -w {output.bracken_report} \\
            -r {params.read_len} \\
            -l {params.level} \\
            2> {log}
        """
