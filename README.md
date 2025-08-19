# zhodiag: Ze shotgun metagenomics pipeline for diagnostics

## How to use:

### Pre-requisites
A conda environment with nextflow version 24.10.5 (other versions not tested.). You can find a working environment called `vcs_nextflow_24.10.5` in asterix:
`conda activate vcs_nextflow_24.10.5`

### Usage
1. Clone the repo
2. Prepare input sample table. This is a comma-separated table with the following columns.
    sample,fastq_1,fastq_2,group. You can see an example in `data/example_groups.csv`

* sample: sample name. This will prefix all files.
* fastq_1 and fastq_2: full path to R1 and R2 reads.
* group: variable used in final plots. You will need at least two groups for the plots to work. A "control" group is useful.

3. Run (with -resume if re-launching):
`nextflow run main.nf -profile singularity --input data/example_groups.csv --outdir OUTPUT -resume`


## Diagram (not up-to-date)

![Diagram in progress](misc/diagram.png)
