process KRAKEN2_TAXONOMY {

    conda "${moduleDir}/environment.yml"
    container './modules/local/kraken2_taxonomy/ete3_custom.sif'

    input:
    tuple val(meta), path(infile)
    val threshold
    val kraken2_db_name
    path db_file

    output:
    tuple val(meta), path("*_taxonomy.tsv") , emit: taxonomy

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${kraken2_db_name}"
    def script_dir = file("bin")
    def kraken_taxonomy_script = script_dir.resolve("kraken2_linear_taxonomy.py").toString()
    def infileName = infile.getName().replaceFirst(/(\.[^.]+)?$/, '')  // remove extension
    def outfile = "${infileName}_taxonomy.tsv"
    
    """
    python3 $kraken_taxonomy_script -i $infile -o $outfile --db $db_file
    """

}
