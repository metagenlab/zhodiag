process MINIMAP2_TAXONOMY {

    tag "Adding full taxonomy"

    conda "${moduleDir}/environment.yml"
    // container './modules/local/kraken2_taxonomy/ete3_custom.sif'
    container "docker://metagenlab/ete:1.0"

    input:
    path minimap2_output

    output:
    path("minimap2_report_taxonomy.tsv") , emit: taxonomy

    script:
    def script_dir = file("bin")
    def minimap_taxonomy_script = script_dir.resolve("minimap2taxonomy.py").toString()

    """
    python3 $minimap_taxonomy_script -i $minimap2_output -o minimap2_report_taxonomy.tsv
    """

}
