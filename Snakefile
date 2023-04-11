import os
import re
import snakemake.io
from glob import glob

sample_list = glob_wildcards("data/reads/{sample}__fwd.fastq.gz").sample

binning_samples = glob_wildcards("data/reads_for_binning/{sample}_{read_dir_coverm}.fastq.gz").sample

binning_example_sample = "samp_447" # I've only been binning one example sample for the working group, which I'll use below for the appropriate binning steps rather than all of the samples from  binning_samples or sample_list

rule all:
    input:
        #expand("data/sourmash/outputs/{sample}_tax.csv", sample = sample_list), # I commented this out to focus on binning
        expand("data/binning/{sample}/METABAT2/.done", sample = binning_example_sample),
        expand("data/binning/{sample}/bin_coverage.tsv", sample = binning_example_sample),
        expand("data/binning/{sample}/checkm.txt", sample = binning_example_sample),
        expand("data/binning/{sample}/gtdbtk", sample = binning_example_sample)


rule run_sourmash:
    input: expand("data/sourmash/outputs/{sample}_tax.csv", sample=sample_list)

rule sourmash_sketch:
    input:
        fwd = "data/reads/{sample}__fwd.fastq.gz",
        rev = "data/reads/{sample}__rev.fastq.gz"
    output: "data/sourmash/signatures/{sample}.sig"
    conda: "config/sourmash.yml"
    resources: cpus=1, mem_mb=20000, time_min=5000
    shell:
        """
        sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge {wildcards.sample} -o {output} {input.fwd} {input.rev}
        """

rule sourmash_gather:
    input:
        sig = rules.sourmash_sketch.output,
        gtdb_refDB = "data/references/sourmash/gtdb-rs207.dna.k31.zip",
        gtdb_taxDB = "data/references/sourmash/gtdb-rs207.taxonomy.sqldb"
    output:
        reps = "data/sourmash/outputs/{sample}_reps.csv",
        tax = "data/sourmash/outputs/{sample}_tax.csv"
    conda: "config/sourmash.yml"
    resources: cpus=1, mem_mb=20000, time_min=5000
    shell:
        """
        sourmash gather {input.sig} {input.gtdb_refDB} -o {output.reps}

        sourmash tax annotate -g {output.reps} -t {input.gtdb_taxDB}

        mv {wildcards.sample}_reps.with-lineages.csv {output.tax}
        """

rule coverm_contig_coverage:
    input:
        reads = expand("data/reads_for_binning/{sample}_{dir}.fastq.gz", sample = binning_samples, dir = ["R1", "R2"]),
        assembly = "data/assemblies/samp_447.fasta"
    output:
        coverage_metabat = "data/binning/{sample}/coverage_metabat.tsv"
    conda: "config/coverm.yml"
    resources: cpus = 16, mem_mb = 100000, time_min = 20000
    shell:
        """
        coverm contig \
            -c {input.reads} \
            -r {input.assembly} \
            -t {resources.cpus} \
            --methods metabat \
            --output-file {output.coverage_metabat}
        """

# I was having some issues with Docker containers not being able to find the correct path, so this finds our current working directory which we 'cd' into in the metabat rule
current_dir = os.getcwd()

rule metabat2:
    input:
        assembly = "data/assemblies/{sample}.fasta",
        coverm_depth = "data/binning/{sample}/coverage_metabat.tsv"
    output:
        #depth = "data/omics/metagenomes/{sample}/bins/jgi_depth_summary.txt",
        done = touch("data/binning/{sample}/METABAT2/.done") # this writes an empty file that marks that this rule was run, can be more stable than tracking outputs that are directories
    params:
        bin_name = directory("data/binning/{sample}/METABAT2/metabat2")
    singularity: "docker://metabat/metabat" # Notice we are using a docker container to 'install' metabat instead of conda. You should be automatically downloaded from this link without any other work on your part
    resources: cpus=16, mem_mb=20000, time_min=2880 # for standard samples
    #resources: cpus=36, mem_mb=150000, time_min=5880 # for coassembly
    shell:
        """
        pwd # print the current directory
        cd {current_dir} # change into the project root dir
        pwd # print the current directory

        metabat2 \
            -i {input.assembly} \
            -a {input.coverm_depth} \
            -o {params.bin_name} \
            -m 2000 \
            -t {resources.cpus} \
            --unbinned
        """

rule checkM:
    input: "data/binning/{sample}/METABAT2/.done"
    output:
        dir = temp(directory("data/binning/{sample}/checkm")),
        results = "data/binning/{sample}/checkm.txt"
    params:
        bin_dir = "data/binning/{sample}/METABAT2"
    conda: "config/checkm.yaml"
    resources: cpus=8, mem_mb=80000, time_min=2880, mem_gb = 80
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.bin_dir} {output.dir}
        """


rule gtdbtk:
    input:
        metabat_done = "data/binning/{sample}/METABAT2/.done",
        refs = "/geomicro/data2/kiledal/references/gtdbtk/release207_v2"
    output: directory("data/binning/{sample}/gtdbtk")
    params:
        bin_dir = "data/binning/{sample}/METABAT2"
    conda: "config/gtdbtk.yaml"
    resources: cpus=1, mem_mb=500000, time_min=2880, mem_gb = 500
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf \
            --extension fa \
            --genome_dir {params.bin_dir} \
            --out_dir {output} \
            --skip_ani_screen \
            --cpus {resources.cpus}
        """

rule coverm_bin_coverage:
    input:
        reads = expand("data/reads_for_binning/{sample}_{dir}.fastq.gz", sample = binning_samples, dir = ["R1", "R2"]),
        metabat_done = "data/binning/{sample}/METABAT2/.done"
    output:
        coverage_metabat = "data/binning/{sample}/bin_coverage.tsv"
    params:
        bin_dir = "data/binning/{sample}/METABAT2"
    conda: "config/coverm.yml"
    resources: cpus = 16, mem_mb = 100000, time_min = 20000
    shell:
        """
        coverm genome \
            -t {resources.cpus} \
            -m relative_abundance mean covered_bases variance length \
            --min-covered-fraction 0 \
            -c {input.reads} \
            --genome-fasta-files {params.bin_dir}/*.fa \
            -o {output}
        """


# rule ref_read_mapping:
#     input:
#         f_reads = "data/reads/{sample}__fwd.fastq.gz",
#         r_reads = "data/reads/{sample}__rev.fastq.gz"
#         ref = "data/reference/blast_queries/{ref_seqs}.fasta"
#     output:
#         temp_bam = temp("data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped_temp.bam"),
#         sam = temp("data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped.sam"),
#         unsorted_bam = temp("data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped_unsorted.bam"),
#         bam = "data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped.bam"
#     conda: "config/conda_yaml/minimap2.yaml"
#     log: "logs/ref_read_mapping/{sample_type}-{sample}.{ref_seqs}.log"
#     benchmark: "benchmarks/ref_read_mapping/{sample_type}-{sample}.{ref_seqs}.tsv"
#     resources: cpus=8
#     shell:
#         """
#         minimap2 \
#             -ax sr \
#             -t {resources.cpus} \
#             --secondary=yes \
#             {input.ref} \
#             {input.f_reads} {input.r_reads} > {output.sam}

#         samtools view -bS {output.sam} > {output.temp_bam}
        
#         filterBam \
#             --in {output.temp_bam} \
#             --out {output.unsorted_bam} \
#             --minCover 50 \
#             --minId 80
        
#         samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
#         samtools index -@ {resources.cpus} {output.bam}
#         """

# rule ref_read_mapping_pileup:
#     input:
#         bam = "data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped.bam",
#         ref = "data/reference/blast_queries/{ref_seqs}.fasta"
#     output:
#         pileup = "data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_pileup.txt"
#     conda: "config/conda_yaml/minimap2.yaml"
#     log: "logs/ref_read_mapping_pileup/{sample_type}-{sample}.{ref_seqs}.log"
#     benchmark: "benchmarks/ref_read_mapping_pileup/{sample_type}-{sample}.{ref_seqs}.tsv"
#     resources: cpus=1
#     shell:
#         """
#         samtools mpileup -f {input.ref} -o {output.pileup} {input.bam}
#         """

# rule run_toxin_gene_read_mapping:
#     input: 
#         #expand("data/omics/metagenomes/{sample}/ref_read_mapping/toxin-genes_mapped.bam", sample = qcd_samples),
#         expand("data/omics/metagenomes/{sample}/ref_read_mapping/toxin-genes_pileup.txt", sample = qcd_samples),
#         expand("data/omics/metatranscriptomes/{sample}/ref_read_mapping/toxin-genes_pileup.txt", sample = qcd_transcript_samples)


rule prodigal_mags:
    input: 
        bin = "data/binning/{sample}/METABAT2/{bin}.fa"
    output:
        genes = "data/prodigal_mags/{sample}/{bin}.gbk",
        proteins = "data/prodigal_mags/{sample}/{bin}.faa"
    conda: "config/prodigal.yaml"
    log: "logs/progdigal_mags/{sample}__{bin}.log"
    shell:
        """
        prodigal \
            -i {input.bin} \
            -a {output.proteins} \
            -d {output.genes} \
            1>{log} 2>&1
        """

rule kofam_scan:
    input:
        genes = rules.prodigal_mags.output.genes,
        profile = "/geomicro/data2/kiledal/GLAMR/data/reference/kegg/kofamscan/profiles",
        ko_list = "/geomicro/data2/kiledal/GLAMR/data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/kofamscan/{sample}/{bin}_kofam_results.txt"
    conda: "config/kofamscan.yaml"
    resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --cpu={resources.cpus}  \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.sample}_kofamscan \
            --ko-list {input.ko_list} {input.genes}
        """

rule run_kofamscan_bins:
    input: expand("data/kofamscan/{sample}/{bin}_kofam_results.txt", bin = glob_wildcards("data/binning/{sample}/METABAT2/{bin}.fa").bin, sample = binning_example_sample)