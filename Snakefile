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
    input:  expand("data/fastqc/{SampleID}_{direction}/report.html", SampleID = samples, direction = directions)

# Rule for running fastqc
rule fastqc:
    input:
        reads = "data/fastqs/{SampleID}_{direction}.fastq.gz"
    output:
        report = "data/fastqc/{SampleID}_{direction}.html"
    params:
        out_dir = "data/fastqc"
    conda:
        "config/conda/fastqc.yaml"
    resources: time_min = 1000, cpus = 8, mem_mb = 50000
    shell:
        """
        fastqc -o {params.out_dir} -t {resources.cpus} {input.reads}
        """