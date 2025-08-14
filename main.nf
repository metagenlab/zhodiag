#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { MULTIQC_COLLECT_REPORTS } from './modules/local/multiqc/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main'
include { BBMAP_BBDUK } from './modules/nf-core/bbmap/bbduk/main'
include { FASTP } from './modules/nf-core/fastp/main'
include { BBMAP_INDEX } from './modules/nf-core/bbmap/index/main'
// include { BOWTIE2_BUILD } from './modules/nf-core/bowtie2/build/main'
include { BBMAP_ALIGN } from './modules/nf-core/bbmap/align/main'
// include { BOWTIE2_ALIGN } from './modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2HOST } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2FUNGI } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2BACTERIA } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2VIRUS } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2EZVIR } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2PROTOZOA } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2FPV } from './modules/nf-core/minimap2/align/main'

include { SAMTOOLS_DEPTH as MINIMAP2BACTERIA_DEPTH} from './modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH as MINIMAP2FPV_DEPTH} from './modules/nf-core/samtools/depth/main'

include { PAF_PREPARE as MINIMAP2PROTOZOA_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { PAF_PREPARE as MINIMAP2FUNGI_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { PAF_PREPARE as MINIMAP2BACTERIA_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { PAF_PREPARE as MINIMAP2VIRUS_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { PAF_PREPARE as MINIMAP2FPV_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'

include { CONCAT_PAFS as MINIMAP2PROTOZOA_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { CONCAT_PAFS as MINIMAP2FUNGI_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { CONCAT_PAFS as MINIMAP2BACTERIA_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { CONCAT_PAFS as MINIMAP2VIRUS_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { CONCAT_PAFS as MINIMAP2FPV_CONCAT_PAFS } from './modules/local/concat_pafs/main'

include { PLOTS_MINIMAP2 as PLOTS_PROTOZOA } from './modules/local/plots_minimap2/main'
include { PLOTS_MINIMAP2 as PLOTS_FUNGI } from './modules/local/plots_minimap2/main'
include { PLOTS_MINIMAP2 as PLOTS_BACTERIA } from './modules/local/plots_minimap2/main'
include { PLOTS_MINIMAP2 as PLOTS_VIRUS } from './modules/local/plots_minimap2/main'
include { PLOTS_MINIMAP2 as PLOTS_FPV } from './modules/local/plots_minimap2/main'

include { SORT_INDEX_BAM } from './modules/local/bam_sort_index/main'
include { KRAKEN2_BUILD } from './modules/nf-core/kraken2/build/main'

include { KRAKEN2_KRAKEN2 as KRAKEN2BACTERIA } from './modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2FPV } from './modules/nf-core/kraken2/kraken2/main'

include { KRAKEN2_COMBINEKREPORTS as KRAKEN2BACTERIA_COMBINEKREPORTS } from './modules/nf-core/krakentools/combinekreports/main'
include { KRAKEN2_COMBINEKREPORTS as KRAKEN2FPV_COMBINEKREPORTS } from './modules/nf-core/krakentools/combinekreports/main'

include { KRAKEN2_COMBINE_REPORTS as KRAKEN2BACTERIA_COMBINE_REPORTS} from './modules/local/kraken2_combineReports/main'
include { KRAKEN2_COMBINE_REPORTS as KRAKEN2FPV_COMBINE_REPORTS} from './modules/local/kraken2_combineReports/main'

include { KRONA_KREPORT2KRONA} from './modules/nf-core/krakentools/kreport2krona/main'

include { KRONA_PLOTS } from './modules/nf-core/krona/ktimporttext/main'

include { PLOTS_KRAKEN2 as PLOTS_KRAKEN2_BACTERIA } from './modules/local/plots_kraken2/main'
include { PLOTS_KRAKEN2 as PLOTS_KRAKEN2_FPV } from './modules/local/plots_kraken2/main'

// include { KRAKEN2_PARSE } from './modules/local/kraken2_parse/main'
// include { KRAKEN2_TAXONOMY as KRAKEN2_TAXONOMY_PRE } from './modules/local/kraken2_taxonomy/main'
// include { KRAKEN2_TAXONOMY as KRAKEN2_TAXONOMY_POST } from './modules/local/kraken2_taxonomy/main'
// include { COUNTS_MERGE as COUNTS_MERGE_PRE } from './modules/local/counts_merge/main'
// include { COUNTS_MERGE as COUNTS_MERGE_POST } from './modules/local/counts_merge/main'
// include { PLOTS_HEATMAP as PLOTS_HEATMAP_PRE } from './modules/local/plots_heatmap/main'
// include { PLOTS_HEATMAP as PLOTS_HEATMAP_POST } from './modules/local/plots_heatmap/main'
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

//    samples.view()

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
        host_map = MINIMAP2HOST(trimmed.reads,
                                params.host_fasta,
                                true, false, false, true, false)
        unmapped = host_map.unmapped
        mapping_logs = host_map.flagstat
    } else {
        error("Unsupported aligner. Options are 'bbmap' and 'minimap2'.")
    }

    // * --- Collect kingdoms of interest ---
    def selected_kingdoms = params.kingdoms.split(',').collect { it.trim().toLowerCase() }
    // dbs = Channel.fromList(selected_kingdoms)
    // dbs = dbs.map({
    //     it -> params.db_mapþ[it]
    //     })
    // dbs.view()
    def flagstat_channels = []

    // * --- Map to selected kingdoms with minimap2 ---
    // BACTERIA
    if (selected_kingdoms.contains('bacteria')) {
        bacteria_map = MINIMAP2BACTERIA(unmapped,
                                        params.bacteria_fasta,
                                        true, false, false, false, true)
        // log for multiqc
        bacteria_map_log = bacteria_map.flagstat
        flagstat_channels << bacteria_map_log
        // annotate paf table and concatenate
        bacteria_paf_ch = bacteria_map.paf.map { tuple ->
            def meta = tuple[0]
            def paf_path = tuple[1]
            return [meta.id, meta.group, paf_path]
        }
        bacteria_paf_ch.set { bacteria_grouped_paf_ch }
        bacteria_annotated_paf = MINIMAP2BACTERIA_ANNOTATE_PAF(bacteria_grouped_paf_ch)
        bacteria_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { bacteria_collected_annotated_paths }
        bacteria_concatenated_pafs = MINIMAP2BACTERIA_CONCAT_PAFS(bacteria_collected_annotated_paths)
        // plots
        if (params.bacteria_annotation) {
            PLOTS_BACTERIA(bacteria_concatenated_pafs.cat, params.mapq_cutoff, params.coverage_cutoff, "bacteria", params.bacteria_annotation)
        }
    }
    // FUNGI + PROTOZOA + VIRUS
    if (selected_kingdoms.contains('fpv')) {
        fpv_map = MINIMAP2FPV(unmapped, 
                                params.fpv_fasta,
                                true, false, false, false, true)
        // log for multiqc
        fpv_map_log = fpv_map.flagstat
        flagstat_channels << fpv_map_log
        // depth
        fpv_depth = MINIMAP2FPV_DEPTH(fpv_map.bam)
        // annotate paf table and concatenate
        fpv_paf_ch = fpv_map.paf.map { tuple ->
            def meta = tuple[0]
            def paf_path = tuple[1]
            return [meta.id, meta.group, paf_path]
        }
        fpv_paf_ch.set { fpv_grouped_paf_ch }
        fpv_annotated_paf = MINIMAP2FPV_ANNOTATE_PAF(fpv_grouped_paf_ch)
        fpv_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { fpv_collected_annotated_paths }
        fpv_concatenated_pafs = MINIMAP2FPV_CONCAT_PAFS(fpv_collected_annotated_paths)
        // plots
        if (params.fpv_annotation) {
            PLOTS_FPV(fpv_concatenated_pafs.cat, params.mapq_cutoff, params.coverage_cutoff, "fpv", params.fpv_annotation)
        }
    }

    // FUNGI
    if (selected_kingdoms.contains('fungi')) {
        fungi_map = MINIMAP2FUNGI(unmapped, 
                                params.fungi_fasta,
                                true, false, false, false, true)
        // log for multiqc
        fungi_map_log = fungi_map.flagstat
        flagstat_channels << fungi_map_log   
        // annotate paf table and concatenate
        fungi_paf_ch = fungi_map.paf.map { tuple ->
            def meta = tuple[0]
            def paf_path = tuple[1]
            return [meta.id, meta.group, paf_path]
        }
        fungi_paf_ch.set { fungi_grouped_paf_ch }
        fungi_annotated_paf = MINIMAP2FUNGI_ANNOTATE_PAF(fungi_grouped_paf_ch)
        fungi_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { fungi_collected_annotated_paths }
        fungi_concatenated_pafs = MINIMAP2FUNGI_CONCAT_PAFS(fungi_collected_annotated_paths)
        // plots
        PLOTS_FUNGI(fungi_concatenated_pafs.cat, params.mapq_cutoff, params.coverage_cutoff, "fungi", params.fungi_annotation)
    }
    // VIRUS
    if (selected_kingdoms.contains('virus')) {
        // REFSEQ
        virus_refseq_map = MINIMAP2VIRUS(unmapped,
                                        params.virus_fasta,
                                        true, false, false, false, true)
        // log for multiqc
        virus_refseq_map_log = virus_refseq_map.flagstat
        flagstat_channels << virus_refseq_map_log
        // annotate paf table and concatenate
        virus_refseq_paf_ch = virus_refseq_map.paf.map { tuple ->
            def meta = tuple[0]
            def paf_path = tuple[1]
            return [meta.id, meta.group, paf_path]
        }
        virus_refseq_paf_ch.set { virus_refseq_grouped_paf_ch }
        virus_refseq_annotated_paf = MINIMAP2VIRUS_ANNOTATE_PAF(virus_refseq_grouped_paf_ch)
        virus_refseq_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { virus_refseq_collected_annotated_paths }
        virus_refseq_concatenated_pafs = MINIMAP2VIRUS_CONCAT_PAFS(virus_refseq_collected_annotated_paths)
        // plots
        if (params.virus_refseq_annotation){
            PLOTS_VIRUS(virus_refseq_concatenated_pafs.cat, params.mapq_cutoff, params.coverage_cutoff, "virus_refseq", params.virus_refseq_annotation)
        }



        // // EZVIR
        // virus_ezvir_map = MINIMAPEZVIR(unmapped,
        //                                 params.ezvir_fasta,
        //                                 true, false, false, false, true)
        // // log for multiqc
        // virus_ezvir_map_log = virus_ezvir_map.flagstat
        // flagstat_channels << virus_ezvir_map_log
    }
    // PROTOZOA
    if (selected_kingdoms.contains('protozoa')) {
        protozoa_map = MINIMAP2PROTOZOA(unmapped,
                                        params.protozoa_fasta,
                                        true, false, false, false, true)
        protozoa_map_log = protozoa_map.flagstat
        // log for multiqc
        flagstat_channels << protozoa_map_log
        // annotate paf table and concatenate
        protozoa_paf_ch = protozoa_map.paf.map { tuple ->
            def meta = tuple[0]
            def paf_path = tuple[1]
            return [meta.id, meta.group, paf_path]
        }
        protozoa_paf_ch.set { protozoa_grouped_paf_ch }
        protozoa_annotated_paf = MINIMAP2PROTOZOA_ANNOTATE_PAF(protozoa_grouped_paf_ch)
        protozoa_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { protozoa_collected_annotated_paths }
        protozoa_concatenated_pafs = MINIMAP2PROTOZOA_CONCAT_PAFS(protozoa_collected_annotated_paths)
        // plots
        PLOTS_PROTOZOA(protozoa_concatenated_pafs.cat, params.mapq_cutoff, params.coverage_cutoff, "protozoa", params.protozoa_annotation)
    }
    all_flagstats = Channel.empty().mix(*flagstat_channels)



    // --- Taxonomic classification with Kraken2 for selected kingdoms ---
    def selected_kingdoms_k2 = params.kraken2_kingdoms.split(',').collect { it.trim().toLowerCase() }
    def kraken_channels = []

    if (selected_kingdoms_k2.contains('fpv')) {
        // get db name
        kraken2_db_name = params.fpv_kraken2db.tokenize('/').last().replaceFirst(/_.*/, '').replaceAll(/[0-9]/, '')
        kraken2_db_name_ch = Channel.value(kraken2_db_name)

        // run kraken2
        fpv_kraken = KRAKEN2FPV(unmapped, 
                                    params.fpv_kraken2db, 
                                    params.kraken2_confidence, 
                                    true, 
                                    true, 
                                    kraken2_db_name_ch)
        fpv_kraken_logs = fpv_kraken.report
        kraken_channels << fpv_kraken_logs

        // Combine reports with krakentools combine_kreports.py
        fpv_kreports_ch = fpv_kraken.report.map { it -> it[1] }
        KRAKEN2FPV_COMBINEKREPORTS(fpv_kreports_ch.collect())

        // Combine reports custom: with group variable.
        // Extract reports as tuples [id, group, report_path]
        fpv_kraken_reports_ch = fpv_kraken.report.map { tuple ->
            def meta = tuple[0]
            def report_path = tuple[1]
            return [meta.id, meta.group, report_path]
        }
        // Collect all reports as list of triplets [id, group, report_path]
        fpv_kraken_reports_ch
            .collect()
            .map { flat_list -> flat_list.collate(3) }
            .set { fpv_grouped_kraken_reports_ch }
        fpv_kraken_reports_combined = KRAKEN2FPV_COMBINE_REPORTS(fpv_grouped_kraken_reports_ch)
        PLOTS_KRAKEN2_FPV(fpv_kraken_reports_combined.combine_long, params.contaminant_taxids)
    }


    // BACTERIA KRAKEN2 PANDB
    if (selected_kingdoms_k2.contains('bacteria')) {
        // get db name
        kraken2_db_name = params.bacteria_kraken2db.tokenize('/').last().replaceFirst(/_.*/, '').replaceAll(/[0-9]/, '')
        kraken2_db_name_ch = Channel.value(kraken2_db_name)

        // run kraken2
        bacteria_kraken = KRAKEN2BACTERIA(unmapped, 
                                    params.bacteria_kraken2db, 
                                    params.kraken2_confidence, 
                                    true, 
                                    true, 
                                    kraken2_db_name_ch)
        bacteria_kraken_logs = bacteria_kraken.report
        kraken_channels << bacteria_kraken_logs

        // Combine reports with krakentools combine_kreports.py
        bacteria_kreports_ch = bacteria_kraken.report.map { tuple ->
            def meta = tuple[0]
            def file = tuple[1]
            return tuple(meta.id, file)
        }
        KRAKEN2BACTERIA_COMBINEKREPORTS(bacteria_kreports_ch.collect())


        // Combine reports custom: with group variable.
        // Extract reports as tuples [id, group, report_path]
        bacteria_kraken_reports_ch = bacteria_kraken.report.map { tuple ->
            def meta = tuple[0]
            def report_path = tuple[1]
            return [meta.id, meta.group, report_path]
        }
        // Collect all reports as list of triplets [id, group, report_path]
        bacteria_kraken_reports_ch
            .collect()
            .map { flat_list -> flat_list.collate(3) }
            .set { bacteria_grouped_kraken_reports_ch }
        bacteria_kraken_reports_combined = KRAKEN2BACTERIA_COMBINE_REPORTS(bacteria_grouped_kraken_reports_ch)
        PLOTS_KRAKEN2_BACTERIA(bacteria_kraken_reports_combined.combine_long, params.contaminant_taxids)
    }

    // --- Taxonomic classification with Kraken2 ---
    // kraken2_db_name = params.kraken2_db.tokenize('/').last().replaceAll(/[^a-zA-Z0-9]/, '')
    // kraken2_db_name_ch = Channel.value(kraken2_db_name)

    // kraken = KRAKEN2_KRAKEN2(unmapped, params.kraken2_db, params.kraken2_confidence, true, true, kraken2_db_name_ch)
    // kraken_logs = kraken.report

    // // Krona plots from krakentools
    // kreport = KRONA_KREPORT2KRONA(kraken.report)
    // KRONA_PLOTS(kreport.txt)

    // // Combine reports with krakentools combine_kreports.py
    // kreports_ch = kraken.report.map { it -> it[1] }
    // KRAKEN2_COMBINEKREPORTS(kreports_ch.collect())

    // // Combine reports custom: with group variable.
    // // Extract reports as tuples [id, group, report_path]
    // kraken_reports_ch = kraken.report.map { tuple ->
    //     def meta = tuple[0]
    //     def report_path = tuple[1]
    //     return [meta.id, meta.group, report_path]
    // }
    // // Collect all reports as list of triplets [id, group, report_path]
    // kraken_reports_ch
    //     .collect()
    //     .map { flat_list -> flat_list.collate(3) }
    //     .set { grouped_kraken_reports_ch }
    // kraken_reports_combined = KRAKEN2_COMBINE_REPORTS(grouped_kraken_reports_ch)
	// PLOTS_KRAKEN2(kraken_reports_combined.combine_long, params.contaminant_taxids)


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


    // Merge all into one channel
    kraken_logs_ch = Channel.empty()
    kraken_channels.each { ch ->
        kraken_logs_ch = kraken_logs_ch.mix(ch)
    }
    kraken_logs_ch = kraken_logs_ch.map { it[1] }

    // --- Collect all reports for MultiQC ---
    collect_reports_input = fastqc.html
        .map { it[1] }
        .merge(
            fastqc.zip.map { it[1] },
            trim_logs.map { it[1] },
            mapping_logs.map { it[1] },
            kraken_logs_ch,
            all_flagstats.map { it[1] }
        )
        .collect()

    multiqc_input = MULTIQC_COLLECT_REPORTS(collect_reports_input)

    MULTIQC(multiqc_input)
}
