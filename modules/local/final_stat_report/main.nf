process FINAL_STAT_REPORT {

    tag "Statistics Report"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    // path(samples)
    path(multiqc_data)
    val host_removal
    path host_fasta
    val krakenuniq
    path krakenuniq_report
    path krakenuniq_kingdoms
    path krakenuniq_removal
    val kraken2
    path kraken2_report
    path kraken2_kingdoms
    path kraken2_removal
    val mapping
    path db_fasta
    path bowtie2_kingdoms
    path bowtie2_removal

    output:
    path '*.tsv'
    path '*.pdf'

    script:
    def custom_report_script = workflow.projectDir.resolve("bin/final_stat_report.r")
    def host_name = host_removal ? host_fasta.getBaseName().replaceFirst(/(\.fna|\.fa|\.fasta)?(\.gz)?$/, '') : 'NA'
    def db_name = mapping ? db_fasta.getBaseName().replaceFirst(/(\.fna|\.fa|\.fasta)?(\.gz)?$/, '') : 'NA'

    """
    Rscript $custom_report_script \
            $multiqc_data \
            $host_removal \
            $host_name \
            $krakenuniq \
            ${krakenuniq ? krakenuniq_report : 'NA'} \
            ${krakenuniq ? krakenuniq_kingdoms : 'NA'} \
            ${krakenuniq ? krakenuniq_removal : 'NA'} \
            $kraken2 \
            ${kraken2 ? kraken2_report : 'NA'} \
            ${kraken2 ? kraken2_kingdoms : 'NA'} \
            ${kraken2 ? kraken2_removal : 'NA'} \
            $mapping \
            $db_name \
            ${mapping ? bowtie2_kingdoms : 'NA'} \
            ${mapping ? bowtie2_removal : 'NA'}
    """

}

