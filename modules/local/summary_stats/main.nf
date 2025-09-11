process SUMMARY_STATS_CANDIDATES {
    tag "$meta.id"
    label 'process_medium'

    container "docker://metagenlab/taxonomizr:1.0"

    input:
    tuple val(meta), path(depth_file), path(map_file)

    output:
    tuple val(meta), path("*by_accession.tsv"), emit: accession_table
    tuple val(meta), path("*by_taxid.tsv"), emit: taxid_table
    tuple val(meta), path("*.pdf")

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def analysis_script = workflow.projectDir.resolve("bin/analysis_by_sample.r")


    """
    awk -F'\\t' '
    {
        split(\$1, a, "|")
        accession = a[1]
        taxid     = a[2]

        key = accession "|" taxid

        # max position
        if (\$2 > max_pos[key]) max_pos[key] = \$2

        # coverage
        if (\$3 > 0) {
            nBases[key]++
            sumDepth[key] += \$3
            countDepth[key]++
        }

        if (!(key in taxid_map)) taxid_map[key] = taxid
    }
    END {
        print "accession\ttaxid\tlength\tnBases_covered\tfraction_covered\tmean_depth"
        for (key in max_pos) {
            split(key, b, "|")
            accession = b[1]
            taxid     = b[2]

            len  = max_pos[key]
            nb   = nBases[key]
            frac = (len > 0 ? nb/len : 0)
            mean = (countDepth[key] > 0 ? sumDepth[key]/countDepth[key] : 0)

            print accession "\\t" taxid "\\t" len "\\t" nb "\\t" frac "\\t" mean
        }
    }' "$depth_file" > ${prefix}_summary_by_accession.tsv

    export TAXONOMY_DB=${workflow.projectDir}/data/accessionTaxa.sql
    Rscript $analysis_script ${prefix}_summary_by_accession.tsv $map_file ${prefix}

    """
}
