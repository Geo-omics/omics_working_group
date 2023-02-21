import os
import re
import snakemake.io
from glob import glob

sample_list = glob_wildcards("data/reads/{sample}__fwd.fastq.gz").sample

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

