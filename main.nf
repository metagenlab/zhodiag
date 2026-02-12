#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'
// qc
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { MULTIQC_COLLECT_REPORTS } from './modules/local/multiqc/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
// trimming
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main'
//include { BBMAP_BBDUK } from './modules/nf-core/bbmap/bbduk/main'
include { FASTP } from './modules/nf-core/fastp/main'

// host removal
//include { BBMAP_INDEX } from './modules/nf-core/bbmap/index/main'
// include { BOWTIE2_BUILD } from './modules/nf-core/bowtie2/build/main'
//include { BBMAP_ALIGN } from './modules/nf-core/bbmap/align/main'
include { BOWTIE2_ALIGN  as BOWTIE2HOST } from './modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2HOST } from './modules/nf-core/minimap2/align/main'

// mapping
include { BOWTIE2_ALIGN  as BOWTIE2DB } from './modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALL } from './modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT as MINIMAP2_SORT} from './modules/nf-core/samtools/sort/main'                                                                                                                       
include { SAMTOOLS_DEPTH as MINIMAP2_DEPTH} from './modules/nf-core/samtools/depth/main'
include { MINIMAP2_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { MINIMAP2_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { MINIMAP2_TAXONOMY } from './modules/local/minimap2_taxonomy/main'
include { PLOTS_MINIMAP2 } from './modules/local/plots_minimap2/main'

//  kraken2
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_COMBINE_REPORTS } from './modules/local/kraken2_combineReports/main'
include { PLOTS_KRAKEN2 } from './modules/local/plots_kraken2/main'
include { KRAKEN2_TAXONOMY } from './modules/local/kraken2_taxonomy/main'

// krakenuniq
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ } from './modules/nf-core/krakenuniq/preloadedkrakenuniq/main'                                                                                     
include { KRAKENUNIQ_COMBINE_REPORTS } from './modules/local/krakenuniq_combineReports/main'
include { PLOTS_KRAKENUNIQ } from './modules/local/plots_krakenuniq/main'

// manual mapping of candidates
// include { EXTRACT_NONHUMAN_READS as MANUAL_CANDIDATES } from './modules/local/extract_nonHuman_reads/main'
// include { MINIMAP2_ALIGN as MINIMAP_CANDIDATES_MANUAL } from './modules/nf-core/minimap2/align/main'
// include { SAMTOOLS_SORT as CANDIDATES_SAMTOOLS_SORT_MANUAL } from './modules/nf-core/samtools/sort/main'                                                                                                                       
// include { SAMTOOLS_DEPTH as CANDIDATES_SAMTOOLS_DEPTH_MANUAL } from './modules/nf-core/samtools/depth/main'
// automatic mapping of candidates
include { EXTRACT_CLASSIFIED } from './modules/local/extract_nonHuman_reads/main'
include { MINIMAP2_ALIGN as MINIMAP_CANDIDATES_AUTO } from './modules/nf-core/minimap2/align/main'
include { SAMTOOLS_DEPTH as CANDIDATES_SAMTOOLS_DEPTH_AUTO } from './modules/nf-core/samtools/depth/main'
include { SAMTOOLS_VIEW as MAP_FILTER_AUTO } from './modules/local/samtools_view/main'
include { SUMMARY_STATS_CANDIDATES as SUMMARY_STATS_CANDIDATES_AUTO } from './modules/local/summary_stats/main'
include { SUMMARY_MAP_CANDIDATES as SUMMARY_MAP_CANDIDATES_AUTO } from './modules/local/summary_map/main'
include { ANALYSIS_COMBINED as ANALYSIS_COMBINED_AUTO } from './modules/local/analysis_combined/main'

include { BOWTIE2_ALIGN  as BOWTIE_CANDIDATES_AUTO } from './modules/nf-core/bowtie2/align/main'


include { CUSTOM_STAT_REPORT } from './modules/local/custom_stat_report/main'

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

//    samples.view()

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
    if (params.mapper == 'minimap2') {
        host_map = MINIMAP2HOST(trimmed.reads,
                                params.host_minimap2_index,
                                true, false, false, true, false)
        unmapped = host_map.unmapped
        mapping_logs = host_map.flagstat
    } else if (params.mapper == 'bowtie2') {
        host_map = BOWTIE2HOST(trimmed.reads,
                                 params.host_bowtie2_index,
                                 params.host_fasta,
                                 true, false)
        unmapped = host_map.fastq
        mapping_logs = host_map.log
    } else {
        error("Unsupported aligner. Options are 'bowtie2' and 'minimap2'.")
    }

    // --------------------------------------------- //
    // --- Taxonomic classification by mapping ---
    // --------------------------------------------- //
    if (params.run_mapping) {
        map = BOWTIE2DB(unmapped,
                        params.bowtie2_index,
                        params.reference_fasta,
                        false, false)
        map_log = map.log
        // map = MINIMAP2_ALL(unmapped,
        //                     params.reference_fasta,
        //                     true, false, false, false, true)
        // // log for multiqc
        // map_log = map.flagstat
        // // samtools sort, index, depth
        // map_sorted = MINIMAP2_SORT(map.bam)
        // MINIMAP2_DEPTH(map_sorted.sorted_bam)
        // // annotate paf table and concatenate
        // paf_ch = map.paf.map { tuple ->
        //     def meta = tuple[0]
        //     def paf_path = tuple[1]
        //     return [meta.id, meta.group, paf_path]
        // }
        // paf_ch.set { grouped_paf_ch }
        // annotated_paf = MINIMAP2_ANNOTATE_PAF(grouped_paf_ch)
        // annotated_paf.paf
        //     .map { it[1] }
        //     .collect()
        //     .set { collected_annotated_paths }
        // concatenated_pafs = MINIMAP2_CONCAT_PAFS(collected_annotated_paths)
        // // add full taxonomy 
        // minimap2_report_taxonomy = MINIMAP2_TAXONOMY(concatenated_pafs.cat)
        // // plots
        // contaminants_ch = params.contaminants ? Channel.fromPath(params.contaminants) :  Channel.value("")
        // PLOTS_MINIMAP2(minimap2_report_taxonomy.taxonomy, 
        //                 params.mapq_cutoff, 
        //                 params.coverage_cutoff, 
        //                 contaminants_ch,
        //                 params.taxonomy_level)
    }

    // --------------------------------------------- //
    // --- Taxonomic classification with Kraken2 ---
    // --------------------------------------------- //
    // run kraken2
    if (params.run_kraken2) {
        kraken = KRAKEN2_KRAKEN2(unmapped, 
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
        // kraken_report_taxonomy = KRAKEN2_TAXONOMY(kraken_reports_combined.combine_long)
        // contaminants_ch = params.contaminants ? Channel.fromPath(params.contaminants) :  Channel.value("")
        results_kraken = PLOTS_KRAKEN2(kraken_reports_combined.combine_long,
                        params.min_reads)
        // KRAKEN2_TAXONOMY(results_kraken.counts, 
        //                     results_kraken.minimizers)
    }

    // --------------------------------------------- //
    // --- Taxonomic classification with KrakenUniq ---
    // --------------------------------------------- //
    if (params.run_krakenuniq) {
        unmapped_for_krakenuniq = unmapped.map { meta, reads ->
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
                                                params.min_reads)
    }


    
    // --------------------------------------------- //
    // --- Mapping classified reads or selected candidates ---
    // --------------------------------------------- //
    if (params.map_classified) {
        // if (params.candidate_mode == 'manual') {
        //     // extract reads from taxid of interest from kraken !!!!! NEEDS UPDATING CODE !!!!
        //     k2_extracted_reads = MANUAL_CANDIDATES(
        //                                 kraken.classified_reads_fastq)
        //     // minimap2 candidates
        //     map_candidates = MINIMAP_CANDIDATES_MANUAL(k2_extracted_reads.extracted_kraken2_reads,
        //                                     params.reference_fasta,
        //                                     false,
        //                                     false,
        //                                     false,
        //                                     false,
        //                                     true)
        //     manual_candidate_mapping_logs = map_candidates.flagstat
        //     candidate_sorted_bam = CANDIDATES_SAMTOOLS_SORT_MANUAL(map_candidates.bam)
        //     CANDIDATES_SAMTOOLS_DEPTH_MANUAL(candidate_sorted_bam.sorted_bam)
        // } else if (params.candidate_mode == 'automatic') {


        // extract all non-human reads from kraken2/uniq
        if (params.which_classified == 'kraken2') {
            k2_extracted_reads = EXTRACT_CLASSIFIED(
                                        kraken.classified_reads_fastq)
        } else if (params.which_classified == 'krakenuniq') {
            k2_extracted_reads = EXTRACT_CLASSIFIED(
                                        krakenuniq.classified_reads)
        }
        // map
        if (params.mapper == 'minimap2') {
            map_candidates = MINIMAP_CANDIDATES_AUTO(k2_extracted_reads.extracted_kraken2_reads,
                                            params.reference_fasta,
                                            false,
                                            false,
                                            false,
                                            false,
                                            false)
            auto_candidate_mapping_logs = map_candidates.flagstat
        } else if (params.mapper == 'bowtie2') {
            map_candidates = BOWTIE_CANDIDATES_AUTO(k2_extracted_reads.extracted_kraken2_reads,
                                                    params.bowtie2_index,
                                                    params.reference_fasta,
                                                    false, true)
            auto_candidate_mapping_logs = map_candidates.log
        }

        // Filter (remove human), analysis
        map_candidates_filter = MAP_FILTER_AUTO(map_candidates.sam,
                                                        params.mapq_cutoff,
                                                        params.coverage_cutoff)
        auto_candidate_mapping_noHuman_logs = map_candidates_filter.flagstat
        // candidate_sorted_bam = CANDIDATES_SAMTOOLS_SORT_AUTO(map_candidates_filter.filtered)
        candidates_depth = CANDIDATES_SAMTOOLS_DEPTH_AUTO(map_candidates_filter.filtered)
        // SLIM_SAM2BAM_AUTO(candidate_sorted_bam.sorted)
        map_summary = SUMMARY_MAP_CANDIDATES_AUTO(map_candidates_filter.filtered)
        // join depth stats and map stats before analysis
        def depth_ch = candidates_depth.depth.map { meta, file ->
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
        analysis_samples = SUMMARY_STATS_CANDIDATES_AUTO(summary_input)
        // join all samples
        stats_accession_ch = analysis_samples.accession_table.map { it -> it[1] }
        stats_taxid_ch = analysis_samples.taxid_table.map { it -> it[1] }
        def metadata_file = file(params.input, checkExists: true)
        Channel.fromPath(metadata_file).set { metadata_ch }
        analysis_combined = ANALYSIS_COMBINED_AUTO(stats_accession_ch.collect(),
                                                    stats_taxid_ch.collect(),
                                                    metadata_ch)
        // }
    }

    // --------------------------------------------- //
    // --- MultiQC ---
    // --------------------------------------------- //
    collect_reports_input = fastqc.html
        .map { it[1] }
        .merge(
            fastqc.zip.map { it[1] },
            trim_logs.map { it[1] },
            mapping_logs.map { it[1] }
        )

    if (params.run_mapping) {
        collect_reports_input = collect_reports_input
            .merge(map_log.map { it[1] } )
    }

    if (params.run_kraken2) {
        collect_reports_input = collect_reports_input
            .merge(kraken_logs.map { it[1] } )
    }

    if (params.run_krakenuniq) {
        collect_reports_input = collect_reports_input
            .merge(krakenuniq_logs.map { it[1] } )
    }

    if (params.map_classified) {
        // if (params.candidate_mode == 'manual') {
        //     collect_reports_input = collect_reports_input
        //         .merge(manual_candidate_mapping_logs.map { it[1] } )
        // } else if (params.candidate_mode == 'automatic') {
        collect_reports_input = collect_reports_input
            .merge(auto_candidate_mapping_logs.map { it[1] } )
        collect_reports_input = collect_reports_input
            .merge(auto_candidate_mapping_noHuman_logs.map { it[1] } )              
        // }
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

    if (params.mapper == "bowtie2") {
        host_name = params.host_bowtie2_index
    } else if (params.mapper == 'minimap2') {
        host_name = params.host_minimap2_index
    }

    CUSTOM_STAT_REPORT(multiqc.data,
                        params.mapper,
                        host_name,
                        params.run_krakenuniq,
                        krakenuniq_report_ch,
                        krakenuniq_kingdoms_ch,
                        krakenuniq_removal_ch,
                        params.run_kraken2,
                        kraken2_report_ch,
                        kraken2_kingdoms_ch,
                        kraken2_removal_ch,
                        params.map_classified
                        )
}
