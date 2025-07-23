#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { COLLECT_REPORTS } from './modules/local/multiqc/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main'
include { BBMAP_BBDUK } from './modules/nf-core/bbmap/bbduk/main'
include { FASTP } from './modules/nf-core/fastp/main'
include { BBMAP_INDEX } from './modules/nf-core/bbmap/index/main'
// include { BOWTIE2_BUILD } from './modules/nf-core/bowtie2/build/main'
include { BBMAP_ALIGN } from './modules/nf-core/bbmap/align/main'
// include { BOWTIE2_ALIGN } from './modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN as MINIMAP_HOST } from './modules/nf-core/minimap2/align/main'
include { SORT_INDEX_BAM } from './modules/local/bam_sort_index/main'
include { KRAKEN2_BUILD } from './modules/nf-core/kraken2/build/main'
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_COMBINEKREPORTS } from './modules/nf-core/krakentools/combinekreports/main'
include { KRAKEN2_COMBINE_REPORTS } from './modules/local/kraken2_combineReports/main'
include { KRAKEN2_KREPORT2KRONA } from './modules/nf-core/krakentools/kreport2krona/main'
include { PLOTS_KRONA } from './modules/nf-core/krona/ktimporttext/main'
include { PLOTS_KRAKEN2 } from './modules/local/plots_kraken2/main'
include { KRAKEN2_PARSE } from './modules/local/kraken2_parse/main'
include { KRAKEN2_TAXONOMY as KRAKEN2_TAXONOMY_PRE } from './modules/local/kraken2_taxonomy/main'
include { KRAKEN2_TAXONOMY as KRAKEN2_TAXONOMY_POST } from './modules/local/kraken2_taxonomy/main'
include { COUNTS_MERGE as COUNTS_MERGE_PRE } from './modules/local/counts_merge/main'
include { COUNTS_MERGE as COUNTS_MERGE_POST } from './modules/local/counts_merge/main'
include { PLOTS_HEATMAP as PLOTS_HEATMAP_PRE } from './modules/local/plots_heatmap/main'
include { PLOTS_HEATMAP as PLOTS_HEATMAP_POST } from './modules/local/plots_heatmap/main'
include { MASH_SCREEN } from './modules/nf-core/mash/screen/main'
include { MASH_SKETCH } from './modules/nf-core/mash/sketch/main'
include { NCBIGENOMEDOWNLOAD } from './modules/nf-core/ncbigenomedownload/main'
include { MINIMAP2_ALIGN as MINIMAP_CANDIDATES } from './modules/nf-core/minimap2/align/main'

workflow {
    // Convert sample sheet CSV into a channel of tuples (meta, reads)
    samples = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    samples = samples.map { row ->
        def sample_meta = row[0]  // e.g. [id: 'sample1']
        def fastq_R1 = row[1]
        def fastq_R2 = row[2]
        def group_label = (row.size() > 3 && row[3]) ? row[3].toString().trim() : "control"
        def single_end = !fastq_R2
        def reads = single_end ? [fastq_R1] : [fastq_R1, fastq_R2]
        def meta = [
            id         : sample_meta.id,
            sample     : sample_meta,
            group      : group_label,
            single_end : single_end
        ]
        return tuple(meta, reads)
    }

    samples.view()

    // --- FASTQC ---
    fastqc = FASTQC(samples)

    // --- Adapter trimming ---
    if (params.trim_tool == "bbduk") {
        trimmed = BBMAP_BBDUK(samples, params.adapters)
        trim_logs = trimmed.stats
    } else if (params.trim_tool == "fastp") {
        trimmed = FASTP(samples, params.adapters, false, false, false)
        trim_logs = trimmed.json
    } else if (params.trim_tool == "trimmomatic") {
        trimmed = TRIMMOMATIC(samples, params.adapters)
        trim_logs = trimmed.out_log
    } else {
        error("Unsupported trim tool. Options are 'bbduk', 'fastp' or 'trimmomatic'.")
    }

    // --- Host removal ---
    if (params.host_removal_tool == 'bbmap') {
        def index
        if (!params.index_host) {
            host_map = BBMAP_ALIGN(trimmed.reads, params.host_bbmap_index)
        } else {
            index = BBMAP_INDEX(params.host_fasta).index
            host_map = BBMAP_ALIGN(trimmed.reads, index)
        }
        unmapped = host_map.unmapped
        mapping_logs = host_map.stats
    } else if (params.host_removal_tool == 'minimap2') {
        host_map = MINIMAP_HOST(trimmed.reads,
                                params.host_fasta,
                                true, false, false, true)
        unmapped = host_map.unmapped
        mapping_logs = host_map.log
    } else {
        error("Unsupported aligner. Options are 'bbmap' and 'minimap2'.")
    }

    // --- Taxonomic classification with Kraken2 ---
    kraken2_db_name = params.kraken2_db.tokenize('/').last()
    kraken2_db_name_ch = Channel.value(kraken2_db_name)

    kraken = KRAKEN2_KRAKEN2(unmapped, params.kraken2_db, params.kraken2_confidence, true, true, kraken2_db_name_ch)
    kraken_logs = kraken.report

    // Krona plots from krakentools
    kreport = KRAKEN2_KREPORT2KRONA(kraken.report)
    PLOTS_KRONA(kreport.txt)

    // Combine reports with krakentools combine_kreports.py
    kreports_ch = kraken.report.map { it -> it[1] }
    KRAKEN2_COMBINEKREPORTS(kreports_ch.collect())

    // Combine reports custom: with group variable.
    // Extract reports as tuples [id, group, report_path]
    kraken_reports_ch = kraken.report.map { tuple ->
        def meta = tuple[0]
        def report_path = tuple[1]
        return [meta.id, meta.group, report_path]
    }

    // Collect all reports as list of triplets [id, group, report_path]
    kraken_reports_ch
        .collect()
        .map { flat_list -> flat_list.collate(3) }
        .set { grouped_kraken_reports_ch }
	grouped_kraken_reports_ch.view()
    kraken_reports_combined = KRAKEN2_COMBINE_REPORTS(grouped_kraken_reports_ch)
	PLOTS_KRAKEN2(kraken_reports_combined.combine_long)

	// kraken_parse = KRAKEN2_PARSE(kraken.classified_reads_assignment, params.kraken2_parse_threshold, kraken2_db_name_ch)
	
	// kraken_taxo_pre = KRAKEN2_TAXONOMY_PRE(kraken_parse.stat_pre, params.kraken2_parse_threshold, kraken2_db_name_ch, file(params.taxonomy_db))
	// kraken_taxo_post = KRAKEN2_TAXONOMY_POST(kraken_parse.stat_post, params.kraken2_parse_threshold, kraken2_db_name_ch, file(params.taxonomy_db))

	// taxonomy_files_pre_ch = kraken_taxo_pre.taxonomy.map { it -> it[1] }
	// counts_pre = COUNTS_MERGE_PRE(taxonomy_files_pre_ch.collect())
	// taxonomy_files_post_ch = kraken_taxo_post.taxonomy.map { it -> it[1] }
	// counts_post = COUNTS_MERGE_POST(taxonomy_files_post_ch.collect())

	// PLOTS_HEATMAP_PRE(counts_pre.counts)
	// PLOTS_HEATMAP_POST(counts_post.counts)


    // --- Mash screen ---
    mash_db_name = params.mash_screen_db.tokenize('/').last()
    mash_db_name_ch = Channel.value(mash_db_name)

    MASH_SCREEN(unmapped, params.mash_screen_db, mash_db_name_ch)

    // --- Minimap2 candidate genomes mapping ---
    if (params.minimap2_candidates) {
        ncbi = NCBIGENOMEDOWNLOAD(params.candidates, params.genomes_filename)

        map_candidates = MINIMAP_CANDIDATES(unmapped,
                                           ncbi.cat_fna,
                                           true,
                                           false,
                                           false,
                                           false)
        SORT_INDEX_BAM(map_candidates.bam)
    }

    // --- Collect all reports for MultiQC ---
    collect_reports_input = fastqc.html
        .map { it[1] }
        .merge(
            fastqc.zip.map { it[1] },
            trim_logs.map { it[1] },
            mapping_logs.map { it[1] },
            kraken_logs.map { it[1] }
        )
        .collect()

    multiqc_input = COLLECT_REPORTS(collect_reports_input)

    MULTIQC(multiqc_input)
}
