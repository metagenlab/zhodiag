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

include { MINIMAP2_ALIGN as MINIMAP2_ALL } from './modules/nf-core/minimap2/align/main'

include { MINIMAP2_ALIGN as MINIMAP2FUNGI } from './modules/nf-core/minimap2/align/main'
// include { MINIMAP2_ALIGN as MINIMAP2BACTERIA } from './modules/nf-core/minimap2/align/main'
// include { MINIMAP2_ALIGN as MINIMAP2FPV } from './modules/nf-core/minimap2/align/main'

include { SAMTOOLS_SORT as SAMTOOLS_SORT} from './modules/nf-core/samtools/sort/main'                                                                                                                       
// include { SAMTOOLS_SORT as MINIMAP2BACTERIA_SORT} from './modules/nf-core/samtools/sort/main'                                                                                                                       
// include { SAMTOOLS_SORT as MINIMAP2FPV_SORT} from './modules/nf-core/samtools/sort/main'                                                                                                                       

include { SAMTOOLS_DEPTH as SAMTOOLS_DEPTH} from './modules/nf-core/samtools/depth/main'
// include { SAMTOOLS_DEPTH as MINIMAP2BACTERIA_DEPTH} from './modules/nf-core/samtools/depth/main'
// include { SAMTOOLS_DEPTH as MINIMAP2FPV_DEPTH} from './modules/nf-core/samtools/depth/main'

include { PAF_PREPARE as MINIMAP2_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'

// include { PAF_PREPARE as MINIMAP2BACTERIA_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
// include { PAF_PREPARE as MINIMAP2FPV_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'

include { CONCAT_PAFS as MINIMAP2_CONCAT_PAFS } from './modules/local/concat_pafs/main'
// include { CONCAT_PAFS as MINIMAP2BACTERIA_CONCAT_PAFS } from './modules/local/concat_pafs/main'
// include { CONCAT_PAFS as MINIMAP2FPV_CONCAT_PAFS } from './modules/local/concat_pafs/main'

include { PLOTS_MINIMAP2 as PLOTS_MINIMAP2 } from './modules/local/plots_minimap2/main'

// include { PLOTS_MINIMAP2 as PLOTS_BACTERIA } from './modules/local/plots_minimap2/main'
// include { PLOTS_MINIMAP2 as PLOTS_FPV } from './modules/local/plots_minimap2/main'

include { SORT_INDEX_BAM } from './modules/local/bam_sort_index/main'
include { KRAKEN2_BUILD } from './modules/nf-core/kraken2/build/main'

include { KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'
// include { KRAKEN2_KRAKEN2 as KRAKEN2BACTERIA } from './modules/nf-core/kraken2/kraken2/main'
// include { KRAKEN2_KRAKEN2 as KRAKEN2FPV } from './modules/nf-core/kraken2/kraken2/main'

include { KRAKEN2_COMBINE_REPORTS as KRAKEN2_COMBINE_REPORTS} from './modules/local/kraken2_combineReports/main'
// include { KRAKEN2_COMBINE_REPORTS as KRAKEN2BACTERIA_COMBINE_REPORTS} from './modules/local/kraken2_combineReports/main'
// include { KRAKEN2_COMBINE_REPORTS as KRAKEN2FPV_COMBINE_REPORTS} from './modules/local/kraken2_combineReports/main'

// include { KRONA_KREPORT2KRONA} from './modules/nf-core/krakentools/kreport2krona/main'

// include { KRONA_PLOTS } from './modules/nf-core/krona/ktimporttext/main'

include { PLOTS_KRAKEN2 as PLOTS_KRAKEN2 } from './modules/local/plots_kraken2/main'
// include { PLOTS_KRAKEN2 as PLOTS_KRAKEN2_BACTERIA } from './modules/local/plots_kraken2/main'
// include { PLOTS_KRAKEN2 as PLOTS_KRAKEN2_FPV } from './modules/local/plots_kraken2/main'

include { MASH_SCREEN } from './modules/nf-core/mash/screen/main'
include { MASH_SKETCH } from './modules/nf-core/mash/sketch/main'

include { KRAKENTOOLS_EXTRACTKRAKENREADS } from './modules/nf-core/krakentools/extractkrakenreads/main'
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

    // * --- MINIMAP2 on all kingdoms ---
    if (params.map_all) {
        map = MINIMAP2_ALL(unmapped,
                                        params.reference_fasta,
                                        true, false, false, false, true)
        // log for multiqc
        minimap_log = map.flagstat
        // annotate paf table and concatenate
        paf_ch = map.paf.map { tuple ->
            def meta = tuple[0]
            def paf_path = tuple[1]
            return [meta.id, meta.group, paf_path]
        }
        paf_ch.set { grouped_paf_ch }
        annotated_paf = MINIMAP2_ANNOTATE_PAF(grouped_paf_ch)
        annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { collected_annotated_paths }
        concatenated_pafs = MINIMAP2_CONCAT_PAFS(collected_annotated_paths)
        // plots
        if (params.reference_annotation) {
            PLOTS_MINIMAP2(concatenated_pafs.cat, 
                            params.mapq_cutoff, 
                            params.coverage_cutoff, 
                            "everything", 
                            params.reference_annotation)
        }
    }
    
    // // * --- MINIMAP2 on kingdoms of interest ---
    // def selected_kingdoms = params.kingdoms.split(',').collect { it.trim().toLowerCase() }
    // def flagstat_channels = []
    // // * --- Map to selected kingdoms with minimap2 ---
    // // BACTERIA
    // if (selected_kingdoms.contains('bacteria')) {
    //     bacteria_map = MINIMAP2BACTERIA(unmapped,
    //                                     params.bacteria_fasta,
    //                                     true, false, false, false, true)
    //     // log for multiqc
    //     bacteria_map_log = bacteria_map.flagstat
    //     flagstat_channels << bacteria_map_log
    //     // annotate paf table and concatenate
    //     bacteria_paf_ch = bacteria_map.paf.map { tuple ->
    //         def meta = tuple[0]
    //         def paf_path = tuple[1]
    //         return [meta.id, meta.group, paf_path]
    //     }
    //     bacteria_paf_ch.set { bacteria_grouped_paf_ch }
    //     bacteria_annotated_paf = MINIMAP2BACTERIA_ANNOTATE_PAF(bacteria_grouped_paf_ch)
    //     bacteria_annotated_paf.paf
    //         .map { it[1] }
    //         .collect()
    //         .set { bacteria_collected_annotated_paths }
    //     bacteria_concatenated_pafs = MINIMAP2BACTERIA_CONCAT_PAFS(bacteria_collected_annotated_paths)
    //     // plots
    //     if (params.bacteria_annotation) {
    //         PLOTS_BACTERIA(bacteria_concatenated_pafs.cat, params.mapq_cutoff, params.coverage_cutoff, "bacteria", params.bacteria_annotation)
    //     }
    // }
    // // FUNGI + PROTOZOA + VIRUS
    // if (selected_kingdoms.contains('fpv')) {
    //     fpv_map = MINIMAP2FPV(unmapped, 
    //                             params.fpv_fasta,
    //                             true, false, false, false, true)
    //     // log for multiqc
    //     fpv_map_log = fpv_map.flagstat
    //     flagstat_channels << fpv_map_log
    //     // sort bam
    //     fpv_sorted_bam = MINIMAP2FPV_SORT(fpv_map.bam)
    //     // depth
    //     fpv_depth = MINIMAP2FPV_DEPTH(fpv_sorted_bam.sorted_bam)
    //     // annotate paf table and concatenate
    //     fpv_paf_ch = fpv_map.paf.map { tuple ->
    //         def meta = tuple[0]
    //         def paf_path = tuple[1]
    //         return [meta.id, meta.group, paf_path]
    //     }
    //     fpv_paf_ch.set { fpv_grouped_paf_ch }
    //     fpv_annotated_paf = MINIMAP2FPV_ANNOTATE_PAF(fpv_grouped_paf_ch)
    //     fpv_annotated_paf.paf
    //         .map { it[1] }
    //         .collect()
    //         .set { fpv_collected_annotated_paths }
    //     fpv_concatenated_pafs = MINIMAP2FPV_CONCAT_PAFS(fpv_collected_annotated_paths)
    //     // plots
    //     if (params.fpv_annotation) {
    //         PLOTS_FPV(fpv_concatenated_pafs.cat, params.mapq_cutoff, params.coverage_cutoff, "fpv", params.fpv_annotation)
    //     }
    // }

    // all_flagstats = Channel.empty().mix(*flagstat_channels)

    ////////////
    // KRAKEN //
    ////////////
    kraken2_db_name = params.kraken2_db.tokenize('/').last() //.replaceFirst(/_.*/, '').replaceAll(/[0-9]/, '')
    kraken2_db_name_ch = Channel.value(kraken2_db_name)

    // run kraken2
    kraken = KRAKEN2_KRAKEN2(unmapped, 
                                params.kraken2_db, 
                                params.kraken2_confidence, 
                                true, 
                                true, 
                                kraken2_db_name_ch)
    // collect log for multiqc
    kraken_logs = kraken.report
    // combine sample reports into one
    kreports_ch = kraken.report.map { it -> it[1] }
    def metadata_file = file(params.input, checkExists: true)
    Channel.fromPath(metadata_file).set { metadata_ch }
    kraken_reports_combined = KRAKEN2_COMBINE_REPORTS(kreports_ch.collect(), metadata_ch)
    PLOTS_KRAKEN2(kraken_reports_combined.combine_long, params.contaminant_taxids)


    // // --- Taxonomic classification with Kraken2 for selected kingdoms ---
    // def selected_kingdoms_k2 = params.kraken2_kingdoms.split(',').collect { it.trim().toLowerCase() }
    // def kraken_channels = []

    // if (selected_kingdoms_k2.contains('fpv')) {
    //     // get db name
    //     kraken2_db_name = params.fpv_kraken2db.tokenize('/').last().replaceFirst(/_.*/, '').replaceAll(/[0-9]/, '')
    //     kraken2_db_name_ch = Channel.value(kraken2_db_name)

    //     // run kraken2
    //     kraken = KRAKEN2FPV(unmapped, 
    //                                 params.fpv_kraken2db, 
    //                                 params.kraken2_confidence, 
    //                                 true, 
    //                                 true, 
    //                                 kraken2_db_name_ch)
    //     kraken_logs = kraken.report
    //     kraken_channels << kraken_logs

    //     // Combine reports
    //     kreports_ch = kraken.report.map { it -> it[1] }
    //     def metadata_file = file(params.input, checkExists: true)
    //     Channel.fromPath(metadata_file).set { metadata_ch }
    //     kraken_reports_combined = KRAKEN2FPV_COMBINE_REPORTS(kreports_ch.collect(), metadata_ch)
    //     PLOTS_KRAKEN2_FPV(kraken_reports_combined.combine_long, params.contaminant_taxids)
    // }


    // // BACTERIA KRAKEN2 PANDB
    // if (selected_kingdoms_k2.contains('bacteria')) {
    //     // get db name
    //     kraken2_db_name = params.bacteria_kraken2db.tokenize('/').last().replaceFirst(/_.*/, '').replaceAll(/[0-9]/, '')
    //     kraken2_db_name_ch = Channel.value(kraken2_db_name)

    //     // run kraken2
    //     kraken = KRAKEN2BACTERIA(unmapped, 
    //                                 params.bacteria_kraken2db, 
    //                                 params.kraken2_confidence, 
    //                                 true, 
    //                                 true, 
    //                                 kraken2_db_name_ch)
    //     kraken_logs = kraken.report
    //     kraken_channels << kraken_logs

    //     // Combine reports
    //     kreports_ch = kraken.report.map { it -> it[1] }
    //     def metadata_file = file(params.input, checkExists: true)
    //     Channel.fromPath(metadata_file).set { metadata_ch }
    //     kraken_reports_combined = KRAKEN2BACTERIA_COMBINE_REPORTS(kreports_ch.collect(), metadata_ch)
    //     PLOTS_KRAKEN2_BACTERIA(kraken_reports_combined.combine_long, params.contaminant_taxids)
    // }


    // --- Mash screen ---
    mash_db_name = params.mash_screen_db.tokenize('/').last()
    mash_db_name_ch = Channel.value(mash_db_name)

    MASH_SCREEN(unmapped, params.mash_screen_db, mash_db_name_ch)

    // --- Minimap2 selected candidates ---
    if (params.minimap2_candidates) {
    // extract reads from kraken
        k2_extracted_reads = KRAKENTOOLS_EXTRACTKRAKENREADS(params.candidates,
                                    kraken.classified_reads_assignment,
                                    kraken.classified_reads_fastq,
                                    kraken.report)
        // minimap2 candidates
        map_candidates = MINIMAP_CANDIDATES(k2_extracted_reads.extracted_kraken2_reads,
                                           params.reference_fasta,
                                           true,
                                           false,
                                           false,
                                           false,
                                           true)
        SAMTOOLS_SORT(map_candidates.bam)
        SAMTOOLS_DEPTH(map_candidates.bam)
    }


    // // Merge all into one channel
    // kraken_logs_ch = Channel.empty()
    // kraken_channels.each { ch ->
    //     kraken_logs_ch = kraken_logs_ch.mix(ch)
    // }
    // kraken_logs_ch = kraken_logs_ch.map { it[1] }

    // --- Collect all reports for MultiQC ---
    collect_reports_input = fastqc.html
        .map { it[1] }
        .merge(
            fastqc.zip.map { it[1] },
            trim_logs.map { it[1] },
            mapping_logs.map { it[1] },
            // kraken_logs_ch,
            kraken_logs.map { it[1] },
            // minimap_log.map { it[1] } // TO BE ADDED WHEN READY
            // all_flagstats.map { it[1] }
        )
        .collect()

    multiqc_input = MULTIQC_COLLECT_REPORTS(collect_reports_input)

    MULTIQC(multiqc_input)
}
