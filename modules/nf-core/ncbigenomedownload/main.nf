process NCBIGENOMEDOWNLOAD {
    // tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ (workflow.containerEngine ?: '') == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.3--pyh7cba7a3_0' :
    'biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0' }"

    input:
    // path(taxids)
    // path accessions
    val taxids
    val fna_cat
    // val groups

    output:
    path("*_genomic.gbff.gz")        , emit: gbk     , optional: true
    path("*_genomic.fna.gz")         , emit: fna     , optional: true
    path("*_rm.out.gz")              , emit: rm      , optional: true
    path("*_feature_table.txt.gz")   , emit: features, optional: true
    path("*_genomic.gff.gz")         , emit: gff     , optional: true
    path("*_protein.faa.gz")         , emit: faa     , optional: true
    path("*_protein.gpff.gz")        , emit: gpff    , optional: true
    path("*_wgsmaster.gbff.gz")      , emit: wgs_gbk , optional: true
    path("*_cds_from_genomic.fna.gz"), emit: cds     , optional: true
    path("*_rna.fna.gz")             , emit: rna     , optional: true
    path("*_rna_from_genomic.fna.gz"), emit: rna_fna , optional: true
    path("*_assembly_report.txt")    , emit: report  , optional: true
    path("*_assembly_stats.txt")     , emit: stats   , optional: true
    path("*_genomes.fna.gz")           , emit: cat_fna
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    // def prefix         = task.ext.prefix ?: "${meta.id}"
    // def accessions_opt = accessions ? "-A ${accessions}" : ""
    def taxids_opt     = taxids ? "-t ${taxids}" : ""
    def fna_output = " && cat *.fna.gz > ${fna_cat}_genomes.fna.gz"
    """
    export XDG_CACHE_HOME=\$PWD/.cache
    
    ncbi-genome-download \\
        $args \\
        $taxids_opt \\
        --flat-output \\
        --parallel $task.cpus \\
        --formats fasta \\
        all \\
        $fna_output

    rm -rf .cache
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$( ncbi-genome-download --version )
    END_VERSIONS
    """
}

        // # $accessions_opt \\
        // --output-folder ./ \\
