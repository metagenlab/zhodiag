# zhodiag: Ze shotgun metagenomics pipeline for diagnostics

## How to use:

### Pre-requisites
A conda environment with nextflow version 24.10.5 (other versions not tested.). You can find a working environment called `vcs_nextflow_24.10.5` in asterix:

```
conda activate vcs_nextflow_24.10.5
```

### Usage
1. Clone the repo
2. Prepare input sample table. This is a comma-separated table with the following columns:
    `sample,fastq_1,fastq_2,group`. You can see an example in `data/example_groups.csv`

* sample: sample name. This will prefix all files.
* fastq_1 and fastq_2: full path to R1 and R2 reads.
* group: variable used in final plots. You will need at least two groups for the plots to work. A "control" group is useful.

3. Edit `nextflow.config` with parameters of choice. The current config works, so you should have **no need to change anything**. 
Minimap2 using a reference fasta of your choice. Currently, the pandb fasta needs fixing so this one will not work (due to duplicated headers). Therefore, for now, I have set the ref fasta to the fpv one. If a reference annotation is provided, plots will be generated.

4. Run (with -resume if re-launching):

```
nextflow run main.nf -profile singularity --input data/example_groups.csv --outdir OUTPUT -resume
```

## Output
The output will be organised by software, for example:

* fastqc: output of qc control step.
* fastp: output of trimming step, including cleaned-reads fastq files and log files.
* minimap2: output of all mapping steps. The "host" subdirectory contains output files after removing human-mapping reads. The rest of the folder contains output of mapping to reference fasta. All mapping steps include flagstat log files and bam files. Reference-mapping step include paf files. 
* kraken2: output of kraken2 classification, including output and report files.
* mash: output of classification with mash.
* plots_and_tables: including summary tables derived from kraken2 (conf subfolders) and minimap2 (mapq/cov subfolders), as well as relevant result plots.
* pipeline_info: nextflow run logs.

## Diagram (not up-to-date)

![Diagram in progress](misc/diagram.png)
