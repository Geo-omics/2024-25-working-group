import os
import re
import subprocess
import snakemake.io
from glob import glob
import pandas as pd
import time

# Get a list of sample names from fastq files
samples = glob_wildcards("data/fastqs/{SampleID}_{direction}.fastq.gz").SampleID

# Get a list of directions from fastq files
directions = glob_wildcards("data/fastqs/{SampleID}_{direction}.fastq.gz").direction

# Create target rule for running fastqc
rule run_fastqc:
    input: expand("data/fastqc/{SampleID}_{direction}_fastqc.html", SampleID = samples, direction = directions)

rule run_fastp:
    input: expand("data/fastqs/{SampleID}_{direction}_qcd.fastq.gz", SampleID = samples, direction = directions)


# Rule for running fastqc
rule fastqc:
    input:
        reads = "data/fastqs/{SampleID}_{direction}.fastq.gz"
    output:
        report = "data/fastqc/{SampleID}_{direction}_fastqc.html"
    params:
        out_dir = "data/fastqc"
    conda:
        "config/conda/fastqc.yaml"
    log: "logs/fastqc/{SampleID}_{direction}.log"
    resources: time_min = 1000, cpus = 8, mem_mb = 50000
    shell:
        """
        fastqc -o {params.out_dir} -t {resources.cpus} {input.reads} | tee {log}
        """

rule fastp:
    input:
        fwd_reads = "data/fastqs/{SampleID}_fwd.fastq.gz",
        rev_reads = "data/fastqs/{SampleID}_rev.fastq.gz"
    output:
        tmp_fwd = temp("data/fastqs/{SampleID}_fwd_dedup_tmp.fastq.gz"),
        tmp_rev = temp("data/fastqs/{SampleID}_rev_dedup_tmp.fastq.gz"),
        rev_reads = "data/fastqs/{SampleID}_rev_qcd.fastq.gz",
        fwd_reads = "data/fastqs/{SampleID}_fwd_qcd.fastq.gz",
        json_dedup = "data/fastp/{SampleID}_dedup.json",
        html_dedup = "data/fastp/{SampleID}_dedup.html",
        html = "data/fastp/{SampleID}.html",
        json = "data/fastp/{SampleID}.json"
    conda: "config/conda/fastp.yaml"
    log: "logs/fastp/{SampleID}.log"
    resources: cpus = 16, mem_mb = 60000, tim_min=2880
    shell: 
        """
        # First deduplicate
        fastp \
            -i {input.fwd_reads} -I {input.rev_reads} \
            -o {output.tmp_fwd} -O {output.tmp_rev} \
            -h {output.html_dedup} -j {output.json_dedup} \
            --thread {resources.cpus} \
            -z 3 \
            --dedup \
            --dup_calc_accuracy 6  2>&1 | tee {log}

        # Trim and filter reads, remove adapters
        fastp \
            -i {output.tmp_fwd} -I {output.tmp_rev} \
            -o {output.fwd_reads} -O {output.rev_reads} \
            -h {output.html} -j {output.json} \
            --thread {resources.cpus} \
            -z 9 \
            --length_required 50 \
            --n_base_limit 5 \
            --low_complexity_filter --complexity_threshold 7 \
            --detect_adapter_for_pe \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --overrepresentation_analysis 2>&1 | tee -a {log}
        """

rule bb_index:
    input:
        human_genome = ancient(rules.get_contaminants.output.human_genome)
    output:
        "data/reference/contaminants/ref/genome/1/summary.txt",
        index = directory("data/reference/contaminants/ref/")
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/bbmap_index.log"
    benchmark:
        "benchmarks/bb_index.txt"
    resources: cpus = 8, mem_mb = 50000
    shell:
        """        
        bbmap.sh \
            ref={input.human_genome} \
            path={params.bbmap_index_path} \
            t={resources.cpus} \
            2>&1 | tee {log}
        """

rule remove_contaminants:
    input:
        rev_reads = "data/fastqs/{SampleID}_rev_qcd.fastq.gz",
        fwd_reads = "data/fastqs/{SampleID}_fwd_qcd.fastq.gz",
        human_genome = "data/references/contaminants/human.fa.gz",
        spike_ins = "data/references/contaminants/spike-ins.fa",
        adapters = "data/references/contaminants/adapters.fa"
        #bbmap_index = "data/references/contaminants/ref"
    output:
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda:
    log:
    resources: mem_mb = 120000
    shell:
        """
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
        echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"
        
        # Remove PhiX reads
        bbduk.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={input.fwd_reads} \
            in2={input.rev_reads} \
            out1={output.phix_rm_fwd} \
            out2={output.phix_rm_rev} \
            t={resources.cpus} k=31 hdist=1 \
            ref={input.skip_ins} \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}

        # Remove Human reads
        echo "\n\n***doing remove contaminants***\n\n" >> {log}




        """