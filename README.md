# :milky_way: zhodiag: Ze shotgun metagenomics pipeline for diagnostics


## :telescope: Pipeline overview
![Diagram in progress](misc/zhodiag_vertical_db_legend.png)


## :microscope: How to use:

### :stars: Pre-requisites
A conda environment with nextflow version 24.10.5 (other versions not tested.). You can find a working environment in asterix and obelix:

```
conda activate vcs_nextflow_24.10.5
```

### :rocket: Usage
1. Clone the repo.
2. Prepare input sample table. This is a comma-separated table with the following columns:
    `sample,fastq_R1,fastq_R2,group`. You can see an example in `data/samples.csv`

* sample: sample name. This will prefix all files.
* fastq_R1 and fastq_R2: full path to R1 and R2 reads.
* group: variable used for faceting in output plots. For control samples, use "control" as group.

3. Edit `nextflow.config` with parameters of choice. 

* you can choose between minimap2 and bowtie2 for mapping steps (option mapper).
* You can optionally run taxonomy classification with kraken2/mash and/or direct mapping by turning on/off (true/false) the corresponding processes.
* Minimap2 uses a reference fasta or mmi index of your choice. Kraken2 used a db of your choice. Bowtie2 requires the reference index and fasta.

4. Run (with -resume if re-launching):

```
nextflow run main.nf -profile singularity --input data/example_groups.csv --outdir OUTPUT -resume
```

## :fireworks: Output
The output will be organised by software, for example:

* fastqc: output of qc control step.
* fastp: output of trimming step, including cleaned-reads fastq files and log files.
* minimap2/bowtie2: output of all mapping steps. The "host" subdirectory contains output files after removing human-mapping reads. The rest of the folder contains output of mapping to the provided reference fasta. All mapping steps include flagstat log files and bam files.
* kraken2: output of kraken2 classification, including output and report files.
* mash: output of classification with mash.
* plots_and_tables: includes summary tables derived from kraken2 and minimap2 processes, as well as relevant result plots. The tables are the input to generate all plots and can be used to manually genereate any downstream plot of interest.
* pipeline_info: nextflow run logs.


