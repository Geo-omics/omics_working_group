import os
import re
import snakemake.io
from glob import glob

sample_list = glob_wildcards("data/reads/{sample}__fwd.fastq.gz").sample

binning_samples = glob_wildcards("data/reads_for_binning/{sample}_{read_dir_coverm}.fastq.gz").sample

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