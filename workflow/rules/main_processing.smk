def get_fastq_for_sample(wildcards):
    """Get input FASTQ files for a sample from the sample sheet."""
    row = samples_df[samples_df["sample"] == wildcards.sample].iloc[0]
    return {
        "r1": row["fastq_1"],
        "r2": row["fastq_2"]
    }

rule fastqc_raw:
    """Run FastQC on raw reads."""
    input:
        unpack(get_fastq_for_sample)
    output:
        html_r1 = "results/fastqc_raw/{sample}_R1_fastqc.html",
        html_r2 = "results/fastqc_raw/{sample}_R2_fastqc.html",
        zip_r1 = "results/fastqc_raw/{sample}_R1_fastqc.zip",
        zip_r2 = "results/fastqc_raw/{sample}_R2_fastqc.zip"
    params:
        outdir = "results/fastqc_raw",
        rename_r1 = "temp/fastqc_raw/{sample}_R1.fastq.gz",
        rename_r2 = "temp/fastqc_raw/{sample}_R2.fastq.gz"
    threads: 2
    resources:
        cpus_per_task=2,
        mem_mb=4000,
        runtime="2h",
        slurm_partition="short"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p temp/fastqc_raw/;
        # Create symbolic links with sample names for better reporting
        ln -sf $(realpath {input.r1}) {params.rename_r1}
        ln -sf $(realpath {input.r2}) {params.rename_r2}

        # Run FastQC
        fastqc -t {threads} -o {params.outdir} {params.rename_r1} {params.rename_r2}

        # Clean up symbolic links
        rm -f {params.rename_r1} {params.rename_r2}
        """

rule fastp:
    """Adapter trimming and quality filtering using fastp."""
    input:
        sample = lambda wildcards: [
            samples_df[samples_df["sample"] == wildcards.sample].iloc[0]["fastq_1"],
            samples_df[samples_df["sample"] == wildcards.sample].iloc[0]["fastq_2"]
        ]
    output:
        trimmed = [
            "results/fastp/{sample}_1.trimmed.fastq.gz",
            "results/fastp/{sample}_2.trimmed.fastq.gz"
        ],
        unpaired = "results/fastp/{sample}.unpaired.fastq.gz",
        json = "results/fastp/{sample}.fastp.json",
        html = "results/fastp/{sample}.fastp.html"
    log:
        "results/logs/fastp/{sample}_fastp.log"
    params:
        adapters = lambda wildcards: (
            f"--adapter_sequence {config['fastp']['adapter_r1']}" if config['fastp']['adapter_r1'] else ""
        ) + (
            f" --adapter_sequence_r2 {config['fastp']['adapter_r2']}" if config['fastp']['adapter_r2'] else ""
        ),
        extra = lambda wildcards:
        (
            " --trim_poly_g" if config["fastp"]["trim_poly_g"] else ""
        ) + (
            " --trim_poly_x" if config["fastp"]["trim_poly_x"] else ""
        ) + (
            f" --complexity_threshold {config['fastp']['complexity_threshold']}"
            if config["fastp"]["complexity_threshold"] else ""
        ) + (
            f" --qualified_quality_phred {config['fastp']['qualified_quality_phred']}"
        ) + (
            f" --length_required {config['fastp']['length_required']}"
        )
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=8000,
        runtime="4h",
        slurm_partition="short"
    wrapper:
        "v7.5.0/bio/fastp"

rule build_hisat2_index:
    """Build HISAT2 index from reference genome."""
    input:
        reference = config["T2T_reference"]
    output:
        done = touch(f"{config['hisat2_index']}.index_done")
    params:
        prefix = config["hisat2_index"]
    threads: 8
    resources:
        cpus_per_task=8,
        mem_mb=32000,
        runtime="12h",
        slurm_partition="short"
    conda:
        "../envs/hisat2.yaml"
    log:
        "results/logs/hisat2_index_build.log"
    shell:
        """
        # Check if index already exists
        if [ ! -f "{params.prefix}.1.ht2" ]; then
            hisat2-build -p {threads} {input.reference} {params.prefix} 2> {log}
        else
            echo "Index already exists, skipping build" > {log}
        fi
        """

rule hisat2_host_removal:
    """Remove human reads using HISAT2 alignment."""
    input:
        r1 = "results/fastp/{sample}_1.trimmed.fastq.gz",
        r2 = "results/fastp/{sample}_2.trimmed.fastq.gz",
        index_done = f"{config['hisat2_index']}.index_done"
    output:
        r1 = "results/hisat2/{sample}_1.human_filtered.fastq.gz",
        r2 = "results/hisat2/{sample}_2.human_filtered.fastq.gz",
        report = "results/hisat2/{sample}_hisat2_report.txt"
    params:
        hisat2_index = config["hisat2_index"],
        score_min = config["hisat2"]["score_min"]
    threads: config["hisat2"]["threads"]
    resources:
        cpus_per_task=config["hisat2"]["threads"],
        mem_mb=32000,
        runtime="36h",
        slurm_partition="medium"
    conda:
        "../envs/hisat2.yaml"
    log:
        "results/logs/hisat2/{sample}_hisat2_host_removal.log"
    shell:
        """
        hisat2 -x {params.hisat2_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            --no-spliced-alignment \
            --score-min {params.score_min} \
            --un-conc-gz results/hisat2/{wildcards.sample}_%.human_filtered.fastq.gz \
            -S /dev/null \
            2> {output.report}
        """


rule fastqc_filtered:
    """Run FastQC on host-filtered reads."""
    input:
        r1 = "results/hisat2/{sample}_1.human_filtered.fastq.gz",
        r2 = "results/hisat2/{sample}_2.human_filtered.fastq.gz"
    output:
        html_r1 = "results/fastqc_filtered/{sample}_R1_fastqc.html",
        html_r2 = "results/fastqc_filtered/{sample}_R2_fastqc.html",
        zip_r1 = "results/fastqc_filtered/{sample}_R1_fastqc.zip",
        zip_r2 = "results/fastqc_filtered/{sample}_R2_fastqc.zip"
    params:
        outdir = "results/fastqc_filtered",
        rename_r1 = "temp/fastqc_filtered/{sample}_R1.fastq.gz",
        rename_r2 = "temp/fastqc_filtered/{sample}_R2.fastq.gz"
    threads: 2
    resources:
        cpus_per_task=2,
        mem_mb=4000,
        runtime="2h",
        slurm_partition="short"
    conda:
        "../envs/fastqc.yaml"
    log:
        "results/logs/fastqc_filtered/{sample}_fastqc.log"
    shell:
        """
        # Create temp directory if it doesn't exist
        mkdir -p temp/fastqc_filtered

        # Create symbolic links with standardized names for better MultiQC reporting
        ln -sf $(realpath {input.r1}) {params.rename_r1}
        ln -sf $(realpath {input.r2}) {params.rename_r2}

        # Run FastQC
        fastqc -t {threads} -o {params.outdir} {params.rename_r1} {params.rename_r2} 2> {log}

        # Clean up symbolic links
        rm -f {params.rename_r1} {params.rename_r2}
        """
