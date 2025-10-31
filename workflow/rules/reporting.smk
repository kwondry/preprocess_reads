rule multiqc:
    input:
        fastqc_raw = expand("results/fastqc_raw/{sample}_{read}_fastqc.html",
                           sample=SAMPLES, read=["R1", "R2"]),
        fastqc_filtered = expand("results/fastqc_filtered/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["R1", "R2"]),
        fastp = expand("results/fastp/{sample}.fastp.json", sample=SAMPLES),
        hisat2 = expand("results/hisat2/{sample}_hisat2_report.txt", sample=SAMPLES),
        kraken = expand("results/kraken/{sample}_kraken_report.txt", sample=SAMPLES)
    output:
        report = "results/multiqc/multiqc_report.html",
        data = directory("results/multiqc/multiqc_data")
    params:
        config_file = "config/multiqc_config.yaml",
        output_dir = "results/multiqc",
        search_dirs = f"results"
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=250000,
        runtime="12h",
        slurm_partition="short"
    conda:
        "../envs/multiqc.yaml"
    log:
        "results/logs/multiqc.log"
    shell:
        """
        # Run MultiQC
        multiqc \\
            --config {params.config_file} \\
            --outdir {params.output_dir} \\
            --force \\
            --interactive \\
            {params.search_dirs} \\
            2> {log}
        """

rule aggregate_read_counts:
    """Create a summary table of read counts at each processing step."""
    input:
        fastp = expand("results/fastp/{sample}.fastp.json", sample=SAMPLES),
        hisat2 = expand("results/hisat2/{sample}_hisat2_report.txt", sample=SAMPLES),
        kraken = expand("results/kraken/{sample}_kraken_report.txt", sample=SAMPLES)
    output:
        summary = "results/read_count_summary.tsv"
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=2000,
        runtime="1h",
        slurm_partition="short"
    run:
        import json
        import pandas as pd
        import re

        data = []
        for sample in SAMPLES:
            row = {"sample": sample}

            # Parse fastp JSON
            fastp_file = f"results/fastp/{sample}.fastp.json"
            with open(fastp_file) as f:
                fastp_data = json.load(f)
                row["raw_reads"] = fastp_data["summary"]["before_filtering"]["total_reads"]
                row["trimmed_reads"] = fastp_data["summary"]["after_filtering"]["total_reads"]

            # Parse HISAT2 report
            hisat2_file = f"results/hisat2/{sample}_hisat2_report.txt"
            with open(hisat2_file) as f:
                content = f.read()

                # Extract total reads
                total_match = re.search(r"(\d+) reads; of these:", content)
                if total_match:
                    total_reads = int(total_match.group(1))

                # Extract reads that aligned 0 times (non-host reads)
                # These are read pairs that did not align to the host genome
                unaligned_match = re.search(r"(\d+) \([\d.]+%\) aligned concordantly 0 times", content)
                if unaligned_match:
                    # This gives us pairs, multiply by 2 for total reads
                    non_host_pairs = int(unaligned_match.group(1))
                    row["non_host_reads"] = non_host_pairs * 2

                # Also extract overall alignment rate for reference
                align_rate_match = re.search(r"([\d.]+)% overall alignment rate", content)
                if align_rate_match:
                    row["host_alignment_rate"] = float(align_rate_match.group(1))

            # Parse Kraken report
            kraken_file = f"results/kraken/{sample}_kraken_report.txt"
            with open(kraken_file) as f:
                lines = f.readlines()
                if lines:
                    # First line contains unclassified reads
                    unclassified = lines[0].strip().split("\t")
                    if len(unclassified) > 1:
                        row["unclassified_percent"] = float(unclassified[0])
                        row["classified_reads"] = int(row.get("non_host_reads", 0) * (1 - float(unclassified[0])/100))

            data.append(row)

        # Create DataFrame and save
        df = pd.DataFrame(data)
        df.to_csv(output.summary, sep="\t", index=False)