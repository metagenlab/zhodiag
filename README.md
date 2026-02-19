# :milky_way: zhodiag: Ze shotgun metagenomics pipeline for diagnostics


## :telescope: Pipeline overview

![Pipeline overview](misc/pipeline.png)


## :microscope: How to use:

### :stars: Pre-requisites
A conda environment with nextflow version 24.10.5 (other versions not tested.). You can find a working environment in asterix and obelix:

```
conda activate vcs_nextflow_24.10.5
```

### :rocket: Usage
1. Clone the repo.
2. Prepare input sample table. This is a comma-separated table with the following columns:
    `sample,fastq_R1,fastq_R2,group`. You can see an example in `data/simdata1.tsv`

* sample: sample name. Unique. This will prefix all output files.
* fastq_R1 and fastq_R2: full path to R1 and R2 reads.
* group: variable used for faceting in output plots. *For control samples, use "control" as group*.

3. Edit `nextflow.config` with parameters of choice. Follow the recomended tools for full functionality (other tools not fully implemented yet)

* Recommended mapper is bowtie2 (option `mapper`). Bowtie2 requires the reference index (option `host_bowtie2_index`) and reference fasta (option `host_fasta`). 
* You can optionally run taxonomy classification with kraken2 and/or krakenuniq by turning on/off (true/false) the corresponding processes. 
* Currently, kraken2 can run with our custom fulldb or with a downloaded db, whereas krakenuniq runs only with the downloaded [MicrobialDB](https://benlangmead.github.io/aws-indexes/k2) (we couldn't compile our custom db into a krakenuniq index).
* After taxonomy assignation, the classified reads can be extracted (option `map_classified` true). You need to choose which classified reads to extract (kraken2 or krakenuniq, depending on which classifier you used, option `which_classified`)
* You can prepare a table file containing known contaminants that you want removed from plots (eg. Bradyrhyzobiums, Cutibacterium acnes, etc). The file is a header-less tsv and has two columns: first column is the taxonomy. This can be a species name, a genus name, or a taxID. The second column species whether is a "species", a "genus", or a "taxid". The species named exactly as writen in the file will be removed. All species starting by the genus as writen in the file will be removed. Taxa specified by the taxID will also be removed. You can pass this file with the parameter `contaminants`. See an example file in `data/contaminants.tsv`.


4. Run (with -resume if re-launching):

```
nextflow run main.nf -profile singularity --input data/simdata1.tsv -resume
```

## :fireworks: Output
The output will be organised by software, for example:

* fastqc: output of qc control step.
* fastp: output of trimming step, including cleaned-reads fastq files and log files.
* bowtie2: output of all mapping steps (host, and optionally, post-classification mapping).
* kraken2: output of kraken2 classification, including output and report files.
* krakenuniq: output of krakenuniq classification, including output and report files.
* plots_and_tables: includes summary tables derived from kraken2, krakenuniq and post-classification mapping processes, as well as relevant result plots. 
A stats_report table is also generated, containing relevant statistics on number of reads mapped, classified, and filtered at each step.
* pipeline_info: nextflow run logs.


