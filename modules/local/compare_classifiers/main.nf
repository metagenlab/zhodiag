process COMPARE_CLASSIFIERS {

    tag "Comparing classifiers"
    container 'docker://rocker/tidyverse:latest'

    input:
    val krakenuniq
    path krakenuniq_report
    val kraken2
    path kraken2_report
    val bowtie2
    path bowtie2_report

    output:
    path '*.tsv'
    path '*.pdf'

    script:
    def custom_report_script = workflow.projectDir.resolve("bin/compare_classifiers.r")

    """
    Rscript $custom_report_script \
            $krakenuniq \
            ${krakenuniq ? krakenuniq_report : 'NA'} \
            $kraken2 \
            ${kraken2 ? kraken2_report : 'NA'} \
            $bowtie2 \
            ${bowtie2 ? bowtie2_report : 'NA'}
    """

}

