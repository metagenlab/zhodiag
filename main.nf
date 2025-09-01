#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'
// qc
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { MULTIQC_COLLECT_REPORTS } from './modules/local/multiqc/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
// trimming
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main'
include { BBMAP_BBDUK } from './modules/nf-core/bbmap/bbduk/main'
include { FASTP } from './modules/nf-core/fastp/main'

// host removal
include { BBMAP_INDEX } from './modules/nf-core/bbmap/index/main'
// include { BOWTIE2_BUILD } from './modules/nf-core/bowtie2/build/main'
include { BBMAP_ALIGN } from './modules/nf-core/bbmap/align/main'
// include { BOWTIE2_ALIGN } from './modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2HOST } from './modules/nf-core/minimap2/align/main'

// minimap2
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

// include { KRONA_KREPORT2KRONA} from './modules/nf-core/krakentools/kreport2krona/main'
// include { KRONA_PLOTS } from './modules/nf-core/krona/ktimporttext/main'

// mash
include { MASH_SCREEN } from './modules/nf-core/mash/screen/main'
include { MASH_SKETCH } from './modules/nf-core/mash/sketch/main'
// manual mapping of candidates
include { KRAKENTOOLS_EXTRACTKRAKENREADS as MANUAL_CANDIDATES } from './modules/nf-core/krakentools/extractkrakenreads/main'
include { MINIMAP2_ALIGN as MINIMAP_CANDIDATES_MANUAL } from './modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT as CANDIDATES_SAMTOOLS_SORT_MANUAL } from './modules/nf-core/samtools/sort/main'                                                                                                                       
include { SAMTOOLS_DEPTH as CANDIDATES_SAMTOOLS_DEPTH_MANUAL } from './modules/nf-core/samtools/depth/main'
// automatic mapping of candidates
include { KRAKENTOOLS_EXTRACTKRAKENREADS as AUTOMATIC_CANDIDATES } from './modules/nf-core/krakentools/extractkrakenreads/main'
include { MINIMAP2_ALIGN as MINIMAP_CANDIDATES_AUTO } from './modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT as CANDIDATES_SAMTOOLS_SORT_AUTO } from './modules/nf-core/samtools/sort/main'                                                                                                                       
include { SAMTOOLS_DEPTH as CANDIDATES_SAMTOOLS_DEPTH_AUTO } from './modules/nf-core/samtools/depth/main'
include { SAMTOOLS_VIEW as MINIMAP_FILTER_AUTO } from './modules/nf-core/samtools/view/main'
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

    // --------------------------------------------- //
    // --- Host removal ---
    // --------------------------------------------- //
    if (params.host_removal_tool == 'bbmap') {
        host_map = BBMAP_ALIGN(trimmed.reads, params.host_bbmap_index)
        unmapped = host_map.unmapped
        mapping_logs = host_map.stats
    } else if (params.host_removal_tool == 'minimap2') {
        host_map = MINIMAP2HOST(trimmed.reads,
                                params.host_minimap2_index,
                                true, false, false, true, false)
        unmapped = host_map.unmapped
        mapping_logs = host_map.flagstat
    } else {
        error("Unsupported aligner. Options are 'bbmap' and 'minimap2'.")
    }

    // --------------------------------------------- //
    // ---------- MINIMAP2 on all kingdoms ---------
    // --------------------------------------------- //
    if (params.run_minimap2) {
        map = MINIMAP2_ALL(unmapped,
                            params.reference_fasta,
                            true, false, false, false, true)
        // log for multiqc
        minimap_log = map.flagstat
        // samtools sort, index, depth
        map_sorted = MINIMAP2_SORT(map.bam)
        MINIMAP2_DEPTH(map_sorted.sorted_bam)
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
        // add full taxonomy 
        minimap2_report_taxonomy = MINIMAP2_TAXONOMY(concatenated_pafs.cat)
        // plots
        contaminants_ch = params.contaminants ? Channel.fromPath(params.contaminants) :  Channel.value("")
        PLOTS_MINIMAP2(minimap2_report_taxonomy.taxonomy, 
                        params.mapq_cutoff, 
                        params.coverage_cutoff, 
                        contaminants_ch,
                        params.taxonomy_level)
    }

    // --------------------------------------------- //
    // --- Taxonomic classification with Kraken2 ---
    // --------------------------------------------- //
    // kraken2_db_name = params.kraken2_db.tokenize('/').last() //.replaceFirst(/_.*/, '').replaceAll(/[0-9]/, '')
    // kraken2_db_name_ch = Channel.value(kraken2_db_name)

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
        contaminants_ch = params.contaminants ? Channel.fromPath(params.contaminants) :  Channel.value("")
        PLOTS_KRAKEN2(kraken_reports_combined.combine_long,
                        contaminants_ch, 
                        params.taxonomy_level)
    }

    // --------------------------------------------- //
    // --- Mash screen ---
    // --------------------------------------------- //
    if (params.run_mash) {
        mash_db_name = params.mash_screen_db.tokenize('/').last()
        mash_db_name_ch = Channel.value(mash_db_name)

        MASH_SCREEN(unmapped, params.mash_screen_db, mash_db_name_ch)
    }

    // --------------------------------------------- //
    // --- Minimap2 selected candidates ---
    // --------------------------------------------- //
    if (params.minimap2_candidates) {
        if (params.candidate_mode == 'manual') {
            // extract reads from taxid of interest from kraken
            k2_extracted_reads = MANUAL_CANDIDATES(params.candidates,
                                        kraken.classified_reads_assignment,
                                        kraken.classified_reads_fastq,
                                        kraken.report)
            // minimap2 candidates
            map_candidates = MINIMAP_CANDIDATES_MANUAL(k2_extracted_reads.extracted_kraken2_reads,
                                            params.reference_fasta,
                                            true,
                                            false,
                                            false,
                                            false,
                                            true)
            manual_candidate_mapping_logs = map_candidates.flagstat
            candidate_sorted_bam = CANDIDATES_SAMTOOLS_SORT_MANUAL(map_candidates.bam)
            CANDIDATES_SAMTOOLS_DEPTH_MANUAL(candidate_sorted_bam.sorted_bam)
        } else if (params.candidate_mode == 'automatic') {
            // extract all reads except human from kraken
            k2_extracted_reads = AUTOMATIC_CANDIDATES("9606",
                                        kraken.classified_reads_assignment,
                                        kraken.classified_reads_fastq,
                                        kraken.report)
            // minimap2 candidates
            map_candidates = MINIMAP_CANDIDATES_AUTO(k2_extracted_reads.extracted_kraken2_reads,
                                            params.reference_fasta,
                                            true,
                                            false,
                                            false,
                                            false,
                                            false)
            auto_candidate_mapping_logs = map_candidates.flagstat
            // map_candidates_filter = MINIMAP_FILTER_AUTO(map_candidates.bam,,,)
            // candidate_sorted_bam = CANDIDATES_SAMTOOLS_SORT_AUTO(map_candidates_filter.bam)
            // candidates_depth = CANDIDATES_SAMTOOLS_DEPTH_AUTO(candidate_sorted_bam.sorted_bam)
        }
    }

    // --------------------------------------------- //
    // --- MultiQC ---
    // --------------------------------------------- //
    // collect_reports_input = fastqc.html
    //     .map { it[1] }
    //     .merge(
    //         fastqc.zip.map { it[1] },
    //         trim_logs.map { it[1] },
    //         mapping_logs.map { it[1] },
    //         kraken_logs.map { it[1] },
    //         minimap_log.map { it[1] }
    //     )
    //     .collect()


    collect_reports_input = fastqc.html
        .map { it[1] }
        .merge(
            fastqc.zip.map { it[1] },
            trim_logs.map { it[1] },
            mapping_logs.map { it[1] }
        )

    if (params.run_minimap2) {
        collect_reports_input = collect_reports_input
            .merge(minimap_log.map { it[1] } )
    }

    if (params.run_kraken2) {
        collect_reports_input = collect_reports_input
            .merge(kraken_logs.map { it[1] } )
    }

    if (params.minimap2_candidates) {
        if (params.candidate_mode == 'manual') {
            collect_reports_input = collect_reports_input
                .merge(manual_candidate_mapping_logs.map { it[1] } )
        } else if (params.candidate_mode == 'automatic') {
            collect_reports_input = collect_reports_input
                .merge(auto_candidate_mapping_logs.map { it[1] } )
        }
    }


    multiqc_input = MULTIQC_COLLECT_REPORTS(collect_reports_input.collect())

        //     kraken_logs.map { it[1] },
        //     minimap_log.map { it[1] }
        // )
        // .collect()

    // def report_channels = [
    //     fastqc.html.map { it[1] },
    //     fastqc.zip.map { it[1] },
    //     trim_logs.map { it[1] },
    //     mapping_logs.map { it[1] }
    //     ]

    // collect_reports_input = Channel
    //     .merge(*report_channels)
    //     .collect()


    MULTIQC(multiqc_input)
}
