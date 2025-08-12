import os
import re
import subprocess
import snakemake.io
from glob import glob
import pandas as pd
import time

# Get a list of sample names from fastq files
#samples = glob_wildcards("data/fastqs/{SampleID}_fwd.fastq.gz").SampleID

samples="samp_447"

genomes = glob_wildcards("data/metabat2/{SampleID}/{genome}.fa").genome

# Get a list of metatrancritome samples
metaT_samples = glob_wildcards("data/transcriptomes/{sample}_fwd.fastq.gz").sample

# Get a list of directions from fastq files
directions = glob_wildcards("data/fastqs/{SampleID}_{direction}.fastq.gz").direction

# Create target rule for running fastqc
rule run_fastqc:
    input: expand("data/fastqc/{SampleID}_{direction}_fastqc.html", SampleID = samples, direction = directions)

rule run_fastp:
    input: expand("data/fastqs/{SampleID}_qcd_{direction}.fastq.gz", SampleID = samples, direction = directions)

rule run_human_read_removal:
    input: 
        expand("data/fastqs/{SampleID}_decon_fwd.fastq.gz", SampleID = samples, direction = directions),
        expand("data/fastqs/{SampleID}_read_count_fastp.tsv", SampleID = samples)

rule run_megahit:
    input: expand("data/assembly/megahit/{SampleID}", SampleID = samples)

rule run_contig_coverage:
    input: expand("data/contig_coverage/{SampleID}-metabat_style_contig_coverage.tsv", SampleID = samples)

rule run_metabat2:
    input: expand("data/metabat2/{SampleID}/.done", SampleID = samples)

rule run_checkm:
    input: expand("data/checkm/{SampleID}.txt", SampleID = samples)

rule run_gtdb:
    input: expand("data/GTDB_{database_version}/{SampleID}/.done_GTDB", SampleID = samples, database_version = "release220")

rule run_bin_abundance:
    input: expand("data/bin_abundance/{SampleID}.txt", SampleID = samples)

rule run_drep:
    input: expand("data/drep/{SampleID}/.drep_done", SampleID = samples)

rule run_kofamscan:
    input: expand("data/metabat2/{SampleID}/kofamscan/{genome}_kofam_results.txt", SampleID = samples, genome = genomes)

rule run_kraken2_reads:
    input: expand("data/kraken2_reads/{SampleID}_kraken2_report.txt", SampleID = samples)

rule run_bracken_reads:
    input: expand("data/kraken2_reads/{SampleID}_bracken_mpa.tsv", SampleID = samples)

rule run_map_metaT_reads:
    input: expand("data/transcriptome_bams/GCF_002095975/{sample}.bam", sample = metaT_samples) # GCF_002095975 is a specific genome, could replace with wildcard to map to multiple genomes

rule run_featurecounts:
    input: expand("data/transcriptome_counts/GCF_002095975/{sample}.counts.txt", sample = metaT_samples) # GCF_002095975 is a specific genome, could replace with wildcard to map to multiple genomes


# Make a graph of all our rules
rule make_rulegraph:
    output:
        "rulegraph.pdf",
        "rulegraph.png"
    shell:
        """
        snakemake data/checkm/samp_447.txt --rulegraph --dry-run | dot -Tpdf > rulegraph.pdf
        """

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
        tmp_fwd = temp("data/fastqs/{SampleID}_dedup_tmp_fwd.fastq.gz"),
        tmp_rev = temp("data/fastqs/{SampleID}_dedup_tmp_rev.fastq.gz"),
        rev_reads = "data/fastqs/{SampleID}_qcd_rev.fastq.gz",
        fwd_reads = "data/fastqs/{SampleID}_qcd_fwd.fastq.gz",
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
        "data/references/contaminants/human.fa.gz",
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
        rev_reads = "data/fastqs/{SampleID}_qcd_rev.fastq.gz",
        fwd_reads = "data/fastqs/{SampleID}_qcd_fwd.fastq.gz",
        human_genome = "data/references/contaminants/human.fa.gz",
        spike_ins = "data/references/contaminants/spike-ins.fa",
        adapters = "data/references/contaminants/adapters.fa"
        #bbmap_index = "data/references/contaminants/ref"
    output:
        phix_rm_fwd = "data/fastqs/{SampleID}_phix_rm_fwd.fastq.gz",
        phix_rm_rev = "data/fastqs/{SampleID}_phix_rm_rev.fastq.gz",
        decon_fwd = "data/fastqs/{SampleID}_decon_fwd.fastq.gz",
        decon_rev = "data/fastqs/{SampleID}_decon_rev.fastq.gz"
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda/bbmap.yaml"
    log: "logs/bbmap/{SampleID}.log"
    resources: mem_mb = 120000, cpus = 24, time_min = 2880
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
            ref={input.spike_ins} \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}

        # Remove Human reads
        echo "\n\n***doing remove contaminants***\n\n" >> {log}

        bbmap.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={output.phix_rm_fwd} \
            in2={output.phix_rm_rev} \
            outu1={output.decon_fwd} \
            outu2={output.decon_rev} \
            ref={input.human_genome} \
            t={resources.cpus} fast=t \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}
        """


rule count_reads_fastp:
    input:
        raw_reads_fwd = "data/fastqs/{SampleID}_fwd.fastq.gz",
        raw_reads_rev = "data/fastqs/{SampleID}_rev.fastq.gz",
        deduped_reads_fwd = "data/fastqs/{SampleID}_dedup_tmp_fwd.fastq.gz",
        deduped_reads_rev = "data/fastqs/{SampleID}_dedup_tmp_rev.fastq.gz",
        qual_filt_and_trimmed_fwd = "data/fastqs/{SampleID}_qcd_fwd.fastq.gz",
        qual_filt_and_trimmed_rev = "data/fastqs/{SampleID}_qcd_rev.fastq.gz",
        decon_reads_fwd = "data/fastqs/{SampleID}_decon_fwd.fastq.gz",
        decon_reads_rev = "data/fastqs/{SampleID}_decon_rev.fastq.gz"
    output:
        "data/fastqs/{SampleID}_read_count_fastp.tsv"
    resources: cpus=4
    shell:
        """
        printf "read_state\tfwd_read_count\trev_read_count\n" > {output} &&
        printf "raw_reads\t$(($(pigz -dc -p {resources.cpus} {input.raw_reads_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.raw_reads_rev} | wc -l) / 4 ))\n" >> {output} &&
        printf "deduped_reads\t$(($(pigz -dc -p {resources.cpus} {input.deduped_reads_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.deduped_reads_rev} | wc -l) / 4 ))\n" >> {output} &&
        printf "filt_and_trimmed_reads\t$(($(pigz -dc -p {resources.cpus} {input.qual_filt_and_trimmed_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.qual_filt_and_trimmed_rev} | wc -l) / 4 ))\n" >> {output} &&
        printf "decon_reads\t$(($(pigz -dc -p {resources.cpus} {input.decon_reads_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.decon_reads_rev} | wc -l) / 4 ))\n" >> {output}
        """

rule megahit:
    input:
        decon_reads_fwd = "data/fastqs/{SampleID}_decon_fwd.fastq.gz",
        decon_reads_rev = "data/fastqs/{SampleID}_decon_rev.fastq.gz"
    output:
        assembly_directory = directory("data/assembly/megahit/{SampleID}")
    log: "logs/megahit/{SampleID}.log"
    conda: "config/conda/megahit.yaml"
    resources: cpus = 24, mem_mb = 500000, time_min = 7200
    shell:
        """
        rm -rf {output.assembly_directory} # for re-running, megahit doesn't overwrite automatically

        megahit \
            -1 {input.decon_reads_fwd} \
            -2 {input.decon_reads_rev} \
            -t {resources.cpus} \
            --presets meta-sensitive \
            -m 0.5 \
            -o {output.assembly_directory} 2>&1 | tee -a {log}
        """

rule run_quast:
    input:
        expand("data/assembly/{SampleID}/quast/report.tsv", SampleID = samples)

rule quast_megahit:
    input:
        megahit_contigs = "data/assembly/megahit/{SampleID}/final.contigs.fa"
    output:
        report = "data/assembly/{SampleID}/quast/report.tsv"
    params:
        out_dir = "data/assembly/{SampleID}/quast"
    log: "log/megahit/quast_megahit/{SampleID}.log"
    benchmark:
        "benchmarks/assembly/quast_megahit/{SampleID}.txt"
    conda:
        "config/conda/quast.yaml"
    resources:
        cpus = 1, mem_mb = 20000
    shell:
        """
        quast.py {input.megahit_contigs} -o {params.out_dir} 2>&1 | tee {log}
        """

rule contig_coverage:
    input:
        decon_reads_fwd = "data/fastqs/{SampleID}_decon_fwd.fastq.gz",
        decon_reads_rev = "data/fastqs/{SampleID}_decon_rev.fastq.gz",
        assembly_directory = "data/assembly/megahit/{SampleID}"
    output: 
        coverage_metabat = "data/contig_coverage/{SampleID}-metabat_style_contig_coverage.tsv"
    params:
        tmpdir = "tmp/coverm_contig_coverage/{SampleID}"
    benchmark: "benchmarks/contig_coverage/{SampleID}.txt"
    conda: "config/conda/coverm.yaml"
    resources: cpus=24, mem_mb=120000, time_min=2880 # standard assemblies
    #resources: cpus=24, mem_mb=1000000, time_min=2880, partition = "largemem" # coassembly
    priority: 2
    shell:
        """
        export TMPDIR={params.tmpdir}
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER
        mkdir -p $TMPDIR

        coverm contig \
            -c data/fastqs/*_decon_fwd.fastq.gz \
            -r {input.assembly_directory}/final.contigs.fa \
            --bam-file-cache-directory $TMPDIR/{wildcards.SampleID} \
            --discard-unmapped \
            -t {resources.cpus} \
            --mapper minimap2-sr \
            --methods metabat \
            --output-file {output.coverage_metabat}

        rm -r $TMPDIR/{wildcards.SampleID}
        """

rule metabat2:
    input:
        assembly_directory = "data/assembly/megahit/{SampleID}",
        coverm_depth = "data/contig_coverage/{SampleID}-metabat_style_contig_coverage.tsv"
    output:
        done = touch("data/metabat2/{SampleID}/.done")
    params:
        bin_name = "data/metabat2/{SampleID}/metabat2"
    benchmark: "benchmarks/metabat2/{SampleID}.txt"
    singularity: "docker://metabat/metabat"
    resources: cpus=16, mem_mb=20000, time_min=2880 # standard samples
    shell:
        """
        metabat2 -i {input.assembly_directory}/final.contigs.fa \
            -a {input.coverm_depth} \
            -o {params.bin_name} \
            -m 2000 \
            -t {resources.cpus} \
            --unbinned
        """

rule checkm:
    input: "data/metabat2/{SampleID}/.done"
    output:
        results = "data/checkm/{SampleID}.txt"
    params:
        in_dir = "data/metabat2/{SampleID}",
        out_dir = "data/checkm/{SampleID}_full_results"
    conda: "config/conda/checkm.yaml"
    resources: cpus=16, mem_mb=80000, time_min=2880
    shell:
        """
        # Move the files corresponding to unbinned contigs to a subfolder so checkM isn't run on them
        mkdir -p data/metabat2/{wildcards.SampleID}/unbinned 
        mv data/metabat2/{wildcards.SampleID}/metabat2.[!0-9]*.fa data/metabat2/{wildcards.SampleID}/unbinned || true

        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.in_dir} {params.out_dir}
        """

rule GTDB_versioned:
    input:
        metabat_complete = "data/metabat2/{SampleID}/.done",
        refs = "/geomicro/data2/kiledal/GLAMR/data/reference/GTDBtk/{database_version}"
    params:
        input_bin_dir = "data/metabat2/{SampleID}",
        out_dir = "data/GTDB_{database_version}/{SampleID}",
        pplacer_cpus = 1
    output:
        done = touch("data/GTDB_{database_version}/{SampleID}/.done_GTDB")
    conda: "config/conda/gtdbtk_2.4.0.yaml"
    benchmark: "benchmarks/GTDB/{SampleID}_database-{database_version}.txt"
    log: "logs/GTDB/{SampleID}_database-{database_version}.log"
    resources: cpus=16, mem_mb=100000, time_min=2880
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf \
            --extension fa \
            --genome_dir {params.input_bin_dir} \
            --out_dir {params.out_dir} \
            --cpus {resources.cpus} \
            --pplacer_cpus {params.pplacer_cpus} \
            --mash_db $GTDBTK_DATA_PATH/mash_db
        """

rule bin_abundance:
    input:
        "data/metabat2/{SampleID}/.done",
        decon_reads_fwd = "data/fastqs/{SampleID}_decon_fwd.fastq.gz",
        decon_reads_rev = "data/fastqs/{SampleID}_decon_rev.fastq.gz"
    output:
        "data/bin_abundance/{SampleID}.txt",
    params:
        bin_dir = "data/metabat2/{SampleID}"
    benchmark: "benchmarks/bin_abundance/{SampleID}.txt"
    conda: "config/conda/coverm.yaml"
    resources: cpus=24, mem_mb=120000, time_min=2880 # standard assemblies
    shell:
        """
        coverm genome \
            -t {resources.cpus} \
            --methods relative_abundance mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm \
            --output-format sparse \
            --min-covered-fraction 0 \
            -1 {input.decon_reads_fwd} \
            -2 {input.decon_reads_rev} \
            --genome-fasta-files {params.bin_dir}/*.fa \
            -o {output}
        """

rule drep_new:
    input: 
        input_bins = "data/metabat2/{SampleID}"
    output:
        touch("data/drep/{SampleID}/.drep_done")
    params:
        main_dir = directory("data/drep/{SampleID}/"),
        input_bins = "data/metabat2/{SampleID}/*.fa"
    conda: "config/conda/drep.yaml"
    log: "logs/drep/{SampleID}.log"
    benchmark: "benchmarks/drep/{SampleID}.txt"
    resources: cpus=8, mem_mb=150000, time_min=2880
    shell:
        """
        rm -rf {params.main_dir} # Clear any old drep output
        
        dRep dereplicate \
            {params.main_dir} \
            -p {resources.cpus} \
            --contamination 50 \
            --completeness 30 \
            -pa 0.9 \
            -sa 0.99 \
            --length 10000 \
            -g {params.input_bins}
        """

rule prodigal:
    input:
        genome = "data/metabat2/{SampleID}/{genome}.fa"
    output: 
        proteins = "data/metabat2/{SampleID}/prodigal/{genome}.faa",
        genes = "data/metabat2/{SampleID}/prodigal/{genome}.gff"
    conda: "config/conda/prodigal.yaml"
    resources: cpus = 1, mem_mb = 10000
    shell:
        """
        prodigal -i {input.genome} -a {output.proteins} -d {output.genes} #1>{log} 2>&1
        """

rule kofam_scan:
    input:
        proteins = "data/metabat2/{SampleID}/prodigal/{genome}.faa",
        profile = "/geomicro/data2/kiledal/GLAMR/data/reference/kegg/kofamscan/profiles",
        ko_list = "/geomicro/data2/kiledal/GLAMR/data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/metabat2/{SampleID}/kofamscan/{genome}_kofam_results.txt"
    conda: "config/conda/kofamscan.yaml"
    #shadow: "shallow"
    benchmark: "benchmarks/kofamscan/{SampleID}-{genome}.txt"
    log: "logs/kofamscan/{SampleID}-{genome}.log"
    resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --format=detail-tsv \
            --cpu={resources.cpus}  \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.SampleID}_kofamscan \
            --ko-list {input.ko_list} {input.proteins} | tee {log}
        """

rule make_blast_db_PD:
    input: 
        "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/SagBay/blast/database/nucl/{blast_name}.fasta"
    output: 
        "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/SagBay/blast/database/nucl/{blast_name}.fasta.nin"
    conda: 
        "config/conda/blast.yaml"
    resources: 
        cpus=1, mem_mb=5000, time_min=10000
    shell: 
        """
        makeblastdb -in {input} -dbtype nucl -logfile {log}
        """

rule blastn_assemblies_UNIVERSAL:
    input: 
        blast_db = "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/SagBay/blast/database/nucl/{blast_name}.fasta",
        blast_db_index = "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/SagBay/blast/database/nucl/{blast_name}.fasta.nin",
        assembly = "/geomicro/data2/pdenuyl2/metagenomics_wg24/data/extract_refine_seqs_gvp/assembly/megahit_noNORM/{sample}_final.contigs.renamed.fa"
    output:
         "data/extract_refine_seqs_gvp/output/blastn/{sample}_{blast_name}_contigs_blastn.tsv"
    conda: 
        "config/conda/blast.yaml"
    resources:
        cpus=4, mem_mb=16000, time_min=10000
    shell:
        """
        blastn -db {input.blast_db} -query {input.assembly} -out {output} -outfmt '6 std qcovs stitle qlen' -num_threads {resources.cpus} -evalue 1e-2
        """

rule run_blastn_assemblies_gvpAC:
    input: 
        expand("data/extract_refine_seqs_gvp/output/blastn/{sample}_{blast_name}_contigs_blastn.tsv", sample = glob_wildcards("/geomicro/data2/pdenuyl2/metagenomics_wg24/2024-25-working-group/data/extract_refine_seqs_gvp/assembly/megahit_noNORM/{sample}_final.contigs.renamed.fa",followlinks=True).sample, blast_name = "blastdbbuoyancygvpAC")

rule kraken2_gvpAC_contigs_pluspf:
    input: 
        contigs = "data/extract_refine_seqs_gvp/contigs/kraken/gvp_contigs/{contig_name}.fa",
        db = "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2024/SagBay/kraken/database_ln/k2_pluspf_20240904"
    output:
        annot = "data/extract_refine_seqs_gvp/output/kraken/gvp_contigs/gvpAC_contigs_pluspf/{contig_name}_kraken2_output.txt",
        report = "data/extract_refine_seqs_gvp/output/kraken/gvp_contigs/gvpAC_contigs_pluspf/{contig_name}_kraken2_report.txt"
    conda: 
        "config/conda/kraken.yaml"
    resources:
        cpus = 12, 
        mem_mb = 64000, 
        time_min = 10000
    shell:
        """
        kraken2 --threads {resources.cpus} \
        --output {output.annot} \
        --report {output.report} \
        --report-minimizer-data \
        --use-names \
        --minimum-hit-groups 3 \gvpAC_contigs_pluspf/
        --db {input.db} \
        {input.contigs}
        """

rule run_kraken2_gvpAC_contigs_pluspf:
    input: 
        expand("data/extract_refine_seqs_gvp/output/kraken/gvp_contigs/gvpAC_contigs_pluspf/{contig_name}_kraken2_output.txt", contig_name = glob_wildcards("data/extract_refine_seqs_gvp/contigs/kraken/gvp_contigs/{sample}.fa",followlinks=True).sample)

rule bakta_mcy_buoyancy_contigs:
    input: 
        fasta = "../output/clinker/gvp_contigs/bakta/input/{contig_name}.fa",
        db = "/geomicro/data2/pdenuyl2/databases/bakta/db"
    output: 
        folder = directory("../output/clinker/gvp_contigs/bakta/output/{contig_name}")
    conda: 
        "config/conda/bakta.yaml"
    resources:
        cpus = 8, 
        mem_mb = 16000, 
        time_min = 10000
    shell:
        """
        bakta \
        --db {input.db} \
        --threads {resources.cpus} \
        --output {output.folder} \
        {input.fasta}

        cp {output.folder}/*.gbff ../output/clinker/gvp_contigs/bakta/output/gbff/ # consilodate gbff files from bakta
        """

rule run_bakta_mcy_buoyancy_contigs:
    input:
        expand("../output/clinker/gvp_contigs/bakta/output/{contig_name}", contig_name = glob_wildcards("../output/clinker/gvp_contigs/bakta/input/{contig}.fa",followlinks=True).contig)

rule clinker_buoyancy_gvpAC_mcy_contigs:
    conda: 
        "config/conda/clinker.yaml"
    shell:
        """
        clinker ../output/clinker/gvp_contigs/bakta/output/gbff/*.gbff \
            -p ../output/clinker/gvp_contigs/buoyancy_gvpAC_mcy_contigs_set41_workinggroup.html
        """

rule run_blastn_assemblies_samp_4400_818442_Lake_Huron_gvpAC:
    input: 
        expand("../output/blastn/{sample}_{blast_name}_contigs_blastn.tsv", sample = glob_wildcards("../glamr_data/assembly/{sample}_final.contigs.renamed.fa", followlinks=True).sample, run = "assemblies_vs_buoyancy", blast_name = "blastdbbuoyancysamp4400818442LH")

rule kraken2_microcystis_set41_contigs_pluspf:
    input: 
        contigs = "kraken/set_41/megahit_assemblies/{contig_name}.fa",
        db = "/geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2024/SagBay/kraken/database_ln/k2_pluspf_20240904"
    output:
        annot = "kraken/set_41/output/set41_contigs_pluspf/{contig_name}_kraken2_output.txt",
        report = "kraken/set_41/output/set41_contigs_pluspf/{contig_name}_kraken2_report.txt"
    conda: 
        "config/kraken.yml"
    resources:
        cpus = 12, 
        mem_mb = 64000, 
        time_min = 10000
    shell:
        """
        kraken2 --threads {resources.cpus} \
        --output {output.annot} \
        --report {output.report} \
        --report-minimizer-data \
        --use-names \
        --minimum-hit-groups 3 \
        --db {input.db} \
        {input.contigs}
        """

rule run_kraken2_microcystis_set41_contigs_pluspf:
    input: 
        expand("kraken/set_41/output/set41_contigs_pluspf/{contig_name}_kraken2_output.txt", contig_name = glob_wildcards("kraken/set_41/megahit_assemblies/{sample}.fa",followlinks=True).sample)



rule kraken2_refseq:
    input:
        fwd_reads = "data/fastqs/{sample}_decon_fwd.fastq.gz",
        rev_reads = "data/fastqs/{sample}_decon_rev.fastq.gz",
        db = "/geomicro/data2/kiledal/GLAMR/data/reference/kraken_databases/core_nt_202505"
    output: 
        kraken_res = "data/kraken2_reads/{sample}_kraken2_output.txt",
        kraken_report = "data/kraken2_reads/{sample}_kraken2_report.txt"
    conda: "config/conda/kraken2.yaml"
    resources:
        cpus = 48, 
        mem_mb = 640000, 
        time_min = 10000
    log: "logs/kraken2/{sample}.log"
    benchmark: "benchmarks/kraken2/{sample}.log"
    shell:
        """
        printf "Running sample: {wildcards.sample} \n"

        kraken2 \
            --threads {resources.cpus} \
            --output {output.kraken_res} \
            --report {output.kraken_report} \
            --report-minimizer-data \
            --minimum-hit-groups 3 \
            --db {input.db} \
            --paired {input.fwd_reads} {input.rev_reads} 2>&1 | tee {log}
        """

rule bracken:
    input:
        report = "data/kraken2_reads/{sample}_kraken2_report.txt",
        db = "/geomicro/data2/kiledal/GLAMR/data/reference/kraken_databases/core_nt_202505",
        kreport2mpa = "code/kreport2mpa.py"
    output:
        bracken_input = "data/kraken2_reads/{sample}_bracken_input.tsv",
        bracken = "data/kraken2_reads/{sample}_bracken.tsv",
        bracken_report = "data/kraken2_reads/{sample}_bracken_report.txt",
        bracken_mpa = "data/kraken2_reads/{sample}_bracken_mpa.tsv"
    params:
        uniq_minimizer_threshold = 150
    conda: "config/conda/kraken2.yaml"
    resources:
        cpus = 1, 
        mem_mb = 24000, 
        time_min = 1000
    log: "logs/bracken/{sample}.log"
    benchmark: "benchmarks/bracken/{sample}.log"
    shell:
        """
        awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {input.report} | cut --complement -f4,5 > {output.bracken_input}

        bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

        ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages
        """


rule map_metaT_reads:
    input:
        ref = "data/reference/genomes/{genome}/{genome}.fna",
        fwd_reads = "data/transcriptomes/{sample}_fwd.fastq.gz",
        rev_reads = "data/transcriptomes/{sample}_rev.fastq.gz"
    output:
        sam = temp("data/transcriptome_bams/{genome}/{sample}.sam"),
        temp_bam = temp("data/transcriptome_bams/{genome}/{sample}.temp.bam"),
        unsorted_bam = temp("data/transcriptome_bams/{genome}/{sample}.unsorted.bam"),
        bam = "data/transcriptome_bams/{genome}/{sample}.bam"
    log: "logs/map_metaT_reads/{genome}__{sample}.log"
    benchmark: "benchmarks/map_metaT_reads/{genome}__{sample}.log"
    conda: "config/conda/minimap2.yaml"
    resources: cpus = 16, mem_mb = 4000, time_min = 180
    shell:
        """
        minimap2 \
            -ax splice:sr \
            {input.ref} \
            {input.fwd_reads} {input.rev_reads} > {output.sam}        # paired-end

        samtools view -bS {output.sam} > {output.temp_bam}

        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 90

        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """


rule featureCounts:
    input:
        bam = "data/transcriptome_bams/{genome}/{sample}.bam",
        gtf = "data/reference/genomes/{genome}/{genome}.gff3"
    output:
        counts = "data/transcriptome_counts/{genome}/{sample}.counts.txt"
    log: "logs/featureCounts/{genome}__{sample}.log"
    conda: "config/conda/featurecounts.yaml"
    resources: cpus = 16, mem_mb = 4000, time_min = 180
    shell:
        """
        featureCounts \
           -a {input.gtf} \ 
           -g ID \
           -o {output.counts} \
           {input.bam}
        """