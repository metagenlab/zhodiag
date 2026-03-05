#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'
// qc
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { MULTIQC_COLLECT_REPORTS } from './modules/local/multiqc/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
// trimming
include { FASTP } from './modules/nf-core/fastp/main'

// host removal
include { BOWTIE2_ALIGN as BOWTIE2HOST } from './modules/nf-core/bowtie2/align/main'

// mapping
include { BOWTIE2_ALIGN as BOWTIE2DB } from './modules/nf-core/bowtie2/align/main'
include { FILTER_SAM2BAM } from './modules/local/filter_sam2bam/main'
include { SAMTOOLS_DEPTH as MAP_DEPTH } from './modules/nf-core/samtools/depth/main'
include { MAP_SUMMARY } from './modules/local/map_summary/main'
include { MAP_SUMMARY_STATS } from './modules/local/map_summary_stats/main'
include { PLOTS_MAPPING } from './modules/local/plots_mapping/main'

//  kraken2
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_COMBINE_REPORTS } from './modules/local/kraken2_combineReports/main'
include { PLOTS_KRAKEN2 } from './modules/local/plots_kraken2/main'

// krakenuniq
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ } from './modules/nf-core/krakenuniq/preloadedkrakenuniq/main'                                                                                     
include { KRAKENUNIQ_COMBINE_REPORTS } from './modules/local/krakenuniq_combineReports/main'
include { PLOTS_KRAKENUNIQ } from './modules/local/plots_krakenuniq/main'

// compare classifiers
include { COMPARE_CLASSIFIERS } from './modules/local/compare_classifiers/main'

// run report
include { FINAL_STAT_REPORT } from './modules/local/final_stat_report/main'

    // --------------------------------------------- //
    // --------------------------------------------- //
    // -------------- WORKFLOW --------------------- //
    // --------------------------------------------- //
    // --------------------------------------------- //

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


    // --------------------------------------------- //
    // --- FASTQC ---
    // --------------------------------------------- //
    fastqc = FASTQC(samples)

    // --------------------------------------------- //
    // --- Adapter trimming ---
    // --------------------------------------------- //
    trimmed = FASTP(samples, params.adapters, false, false, false)
    trim_logs = trimmed.json

    // --------------------------------------------- //
    // --- Host removal ---
    // --------------------------------------------- //
    if (params.host_removal) {
        host_map = BOWTIE2HOST(trimmed.reads,
                                 params.host_bowtie2_index,
                                 params.host_fasta,
                                 true, true)
        unmapped = host_map.fastq
        mapping_logs = host_map.log
    }

    // --------------------------------------------- //
    // --- Input for Taxonomic classification ---
    // --------------------------------------------- //
    if (params.host_removal) {
        input_reads = unmapped
    } else {
        input_reads = trimmed.reads
    }
    nclassifiers = 0

    // --------------------------------------------- //
    // --- Taxonomic classification by mapping ---
    // --------------------------------------------- //
    if (params.run_mapping) {
        nclassifiers += 1
        map = BOWTIE2DB(input_reads,
                        params.reference_bowtie2_index,
                        params.reference_fasta,
                        true, false)
        map_log = map.log
        // Filter out human hits; convert to slim bam
        map_filter = FILTER_SAM2BAM(map.sam,
                                    params.mapq_cutoff,
                                    params.coverage_cutoff)
        map_noHuman_log = map_filter.flagstat
        // Calculate depth and map summary
        map_depth = MAP_DEPTH(map_filter.filtered)
        map_summary = MAP_SUMMARY(map_filter.filtered)
        // join depth stats and map stats before analysis
        def depth_ch = map_depth.depth.map { meta, file ->
            [ meta.id, [meta, file] ]
        }
        def map_ch = map_summary.mapSummary.map { meta, file ->
            [ meta.id, [meta, file] ]
        }
        depth_ch.join(map_ch)
            .map { id, a, b ->
                def (meta1, depth_file) = a
                def (meta2, map_file)  = b
                tuple(meta1, depth_file, map_file)
            }
            .set { summary_input }
        analysis_samples = MAP_SUMMARY_STATS(summary_input)
        // join all samples
        stats_accession_ch = analysis_samples.accession_table.map { it -> it[1] }
        stats_taxid_ch = analysis_samples.taxid_table.map { it -> it[1] }
        def metadata_file = file(params.input, checkExists: true)
        Channel.fromPath(metadata_file).set { metadata_ch }
        results_mapping = PLOTS_MAPPING(metadata_ch,
                                            params.min_reads,
                                            params.contaminants,
                                            stats_accession_ch.collect(),
                                            stats_taxid_ch.collect()
                                            )
    }

    // --------------------------------------------- //
    // --- Taxonomic classification with Kraken2 ---
    // --------------------------------------------- //
    // run kraken2
    if (params.run_kraken2) {
        nclassifiers += 1
        kraken = KRAKEN2_KRAKEN2(input_reads, 
                                    params.kraken2_db, 
                                    params.kraken2_confidence, 
                                    true, 
                                    true)

        // collect log for multiqc
        kraken_logs = kraken.report

        // combine sample reports into one
        kreports_ch = kraken.report.map { it -> it[1] }
        def metadata_file = file(params.input, checkExists: true)
        Channel.fromPath(metadata_file).set { metadata_ch }
        kraken_reports_combined = KRAKEN2_COMBINE_REPORTS(kreports_ch.collect(), metadata_ch)

        // plot results
        results_kraken = PLOTS_KRAKEN2(kraken_reports_combined.combine_long,
                        params.min_reads,
                        params.contaminants)
    }

    // --------------------------------------------- //
    // --- Taxonomic classification with KrakenUniq ---
    // --------------------------------------------- //
    if (params.run_krakenuniq) {
        nclassifiers += 1
        unmapped_for_krakenuniq = input_reads.map { meta, reads ->
            tuple(
                meta,
                reads,
                [ meta.id ]   // add prefixes for krakenuniq module...
            )
        }
        krakenuniq = KRAKENUNIQ_PRELOADEDKRAKENUNIQ(unmapped_for_krakenuniq,
                                                    'fastq',
                                                    params.krakenuniq_db,
                                                    true, true, true)
        // collect log for multiqc
        // multiqc cannot use krakenuniq log: would need a module to modify report to be readable by multiqc....
        krakenuniq_logs = krakenuniq.report

        // combine sample reports into one
        ku_reports_ch = krakenuniq.report.map { it -> it[1] }
        def ku_metadata_file = file(params.input, checkExists: true)
        Channel.fromPath(ku_metadata_file).set { ku_metadata_ch }
        krakenuniq_reports_combined = KRAKENUNIQ_COMBINE_REPORTS(ku_reports_ch.collect(), ku_metadata_ch)

        // plot
        results_krakenuniq = PLOTS_KRAKENUNIQ(krakenuniq_reports_combined.combine_long,
                                                params.min_reads,
                                                params.contaminants)
    }


    // --------------------------------------------- //
    // --- COMPARISON OF CLASSIFIERS ---
    // --------------------------------------------- //
    if (nclassifiers >= 2) {
        krakenuniq_ch = params.run_krakenuniq  ? results_krakenuniq.clean_reads : []
        kraken2_ch    = params.run_kraken2     ? results_kraken.clean_reads     : []
        mapping_ch    = params.run_mapping     ? results_mapping.clean_reads    : []

        COMPARE_CLASSIFIERS(params.run_krakenuniq,
                            krakenuniq_ch,
                            params.run_kraken2,
                            kraken2_ch,
                            params.run_mapping,
                            mapping_ch)
    }
    
    // --------------------------------------------- //
    // --- MultiQC ---
    // --------------------------------------------- //
    collect_reports_input = fastqc.html
        .map { it[1] }
        .merge(
            fastqc.zip.map { it[1] },
            trim_logs.map { it[1] }
        )

    if (params.host_removal) {
        collect_reports_input = collect_reports_input
            .merge(mapping_logs.map { it[1] } )
    }

    if (params.run_mapping) {
        collect_reports_input = collect_reports_input
            .merge(
                map_log.map { it[1] },
                map_noHuman_log.map { it [1] }
                )
    }

    if (params.run_kraken2) {
        collect_reports_input = collect_reports_input
            .merge(kraken_logs.map { it[1] } )
    }

    if (params.run_krakenuniq) {
        collect_reports_input = collect_reports_input
            .merge(krakenuniq_logs.map { it[1] } )
    }

    multiqc_input = MULTIQC_COLLECT_REPORTS(collect_reports_input.collect())

    multiqc = MULTIQC(multiqc_input)


    // --------------------------------------------- //
    // --- Custom Stat Report ---
    // --------------------------------------------- //
    krakenuniq_report_ch = params.run_krakenuniq ? krakenuniq_reports_combined.combine_long : []
    krakenuniq_kingdoms_ch = params.run_krakenuniq ? results_krakenuniq.kingdoms : []
    krakenuniq_removal_ch = params.run_krakenuniq ? results_krakenuniq.removedReadsFromPlots : []

    kraken2_report_ch = params.run_kraken2 ? kraken_reports_combined.combine_long : []
    kraken2_kingdoms_ch = params.run_kraken2 ? results_kraken.kingdoms : []
    kraken2_removal_ch = params.run_kraken2 ? results_kraken.removedReadsFromPlots : []

    host_name = params.host_removal ? params.host_bowtie2_index : []
    db_name = params.run_mapping ? params.reference_bowtie2_index : []

    bowtie2_kingdoms_ch = params.run_mapping ? results_mapping.kingdoms : []
    bowtie2_removal_ch = params.run_mapping ? results_mapping.removedReadsFromPlots : []

    FINAL_STAT_REPORT(multiqc.data,
                        params.host_removal,
                        host_name,
                        params.run_krakenuniq,
                        krakenuniq_report_ch,
                        krakenuniq_kingdoms_ch,
                        krakenuniq_removal_ch,
                        params.run_kraken2,
                        kraken2_report_ch,
                        kraken2_kingdoms_ch,
                        kraken2_removal_ch,
                        params.run_mapping,
                        db_name,
                        bowtie2_kingdoms_ch,
                        bowtie2_removal_ch
                        )
}
  
    // --------------------------------------------- //
    // --- Mapping classified reads or selected candidates ---
    // --------------------------------------------- //
    // if (params.map_classified) {
    //     // if (params.candidate_mode == 'manual') {
    //     //     // extract reads from taxid of interest from kraken !!!!! NEEDS UPDATING CODE !!!!
    //     //     k2_extracted_reads = MANUAL_CANDIDATES(
    //     //                                 kraken.classified_reads_fastq)
    //     //     // minimap2 candidates
    //     //     map_candidates = MINIMAP_CANDIDATES_MANUAL(k2_extracted_reads.extracted_kraken2_reads,
    //     //                                     params.reference_fasta,
    //     //                                     false,
    //     //                                     false,
    //     //                                     false,
    //     //                                     false,
    //     //                                     true)
    //     //     manual_candidate_mapping_logs = map_candidates.flagstat
    //     //     candidate_sorted_bam = CANDIDATES_SAMTOOLS_SORT_MANUAL(map_candidates.bam)
    //     //     CANDIDATES_SAMTOOLS_DEPTH_MANUAL(candidate_sorted_bam.sorted_bam)
    //     // } else if (params.candidate_mode == 'automatic') {


    //     // extract all non-human reads from kraken2/uniq
    //     if (params.which_classified == 'kraken2') {
    //         k2_extracted_reads = EXTRACT_CLASSIFIED(
    //                                     kraken.classified_reads_fastq)
    //     } else if (params.which_classified == 'krakenuniq') {
    //         k2_extracted_reads = EXTRACT_CLASSIFIED(
    //                                     krakenuniq.classified_reads)
    //     }
    //     // map
    //     if (params.mapper == 'minimap2') {
    //         map_candidates = MINIMAP_CANDIDATES_AUTO(k2_extracted_reads.extracted_kraken2_reads,
    //                                         params.reference_fasta,
    //                                         false,
    //                                         false,
    //                                         false,
    //                                         false,
    //                                         false)
    //         auto_candidate_mapping_logs = map_candidates.flagstat
    //     } else if (params.mapper == 'bowtie2') {
    //         map_candidates = BOWTIE_CANDIDATES_AUTO(k2_extracted_reads.extracted_kraken2_reads,
    //                                                 params.bowtie2_index,
    //                                                 params.reference_fasta,
    //                                                 false, true)
    //         auto_candidate_mapping_logs = map_candidates.log
    //     }

    //     // Filter (remove human), analysis
    //     map_candidates_filter = MAP_FILTER_AUTO(map_candidates.sam,
    //                                                     params.mapq_cutoff,
    //                                                     params.coverage_cutoff)
    //     auto_candidate_mapping_noHuman_logs = map_candidates_filter.flagstat
    //     // candidate_sorted_bam = CANDIDATES_SAMTOOLS_SORT_AUTO(map_candidates_filter.filtered)
    //     candidates_depth = CANDIDATES_SAMTOOLS_DEPTH_AUTO(map_candidates_filter.filtered)
    //     // SLIM_SAM2BAM_AUTO(candidate_sorted_bam.sorted)
    //     map_summary = SUMMARY_MAP_CANDIDATES_AUTO(map_candidates_filter.filtered)
    //     // join depth stats and map stats before analysis
    //     def depth_ch = candidates_depth.depth.map { meta, file ->
    //         [ meta.id, [meta, file] ]
    //     }
    //     def map_ch = map_summary.mapSummary.map { meta, file ->
    //         [ meta.id, [meta, file] ]
    //     }
    //     depth_ch.join(map_ch)
    //         .map { id, a, b ->
    //             def (meta1, depth_file) = a
    //             def (meta2, map_file)  = b
    //             tuple(meta1, depth_file, map_file)
    //         }
    //         .set { summary_input }
    //     analysis_samples = SUMMARY_STATS_CANDIDATES_AUTO(summary_input)
    //     // join all samples
    //     stats_accession_ch = analysis_samples.accession_table.map { it -> it[1] }
    //     stats_taxid_ch = analysis_samples.taxid_table.map { it -> it[1] }
    //     def metadata_file = file(params.input, checkExists: true)
    //     Channel.fromPath(metadata_file).set { metadata_ch }
    //     analysis_combined = ANALYSIS_COMBINED_AUTO(stats_accession_ch.collect(),
    //                                                 stats_taxid_ch.collect(),
    //                                                 metadata_ch)
        // }
    // }

