process KRAKEN2_COMBINE_REPORTS2 {
    tag "combining kraken2 reports"
    container 'docker://lbwang/rocker-genome:latest'

    input:
    path k2reports
    path(metadata)

    output:
    path "*/*/kraken2_combined*.tsv", emit: combine_long

    script:
    def kraken_merge_script_long = workflow.projectDir.resolve("bin/combine_kraken_reports.r")
    def db   = k2reports[0].getFileName().toString().split('_')[1]
    def conf = k2reports[0].getFileName().toString().split('_')[2]
    def out_dir = "${db}/${conf}"

    """
    mkdir -p ${out_dir}
    Rscript $kraken_merge_script_long ${metadata} ${out_dir}/kraken2_combined_${db}_${conf} ${k2reports.join(' ')} 
    """


}

