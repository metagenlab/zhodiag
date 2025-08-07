#!/usr/bin/env nextflow
nextflow.enable.dsl=2
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
include { MINIMAP2_ALIGN as MINIMAPHOST } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAPFUNGI } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAPBACTERIA } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAPVIRUS } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAPEZVIR } from './modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAPPROTOZOA } from './modules/nf-core/minimap2/align/main'
include { PAF_PREPARE as MINIMAPPROTOZOA_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { PAF_PREPARE as MINIMAPFUNGI_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { PAF_PREPARE as MINIMAPBACTERIA_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'
include { PAF_PREPARE as MINIMAPVIRUS_ANNOTATE_PAF } from './modules/local/minimap2_pafPrepare/main'

include { CONCAT_PAFS as MINIMAPPROTOZOA_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { CONCAT_PAFS as MINIMAPFUNGI_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { CONCAT_PAFS as MINIMAPBACTERIA_CONCAT_PAFS } from './modules/local/concat_pafs/main'
include { CONCAT_PAFS as MINIMAPVIRUS_CONCAT_PAFS } from './modules/local/concat_pafs/main'

include { PLOTS_MINIMAP2 as PLOTS_PROTOZOA } from './modules/local/plots_minimap2/main'
include { PLOTS_MINIMAP2 as PLOTS_FUNGI } from './modules/local/plots_minimap2/main'
include { PLOTS_MINIMAP2 as PLOTS_BACTERIA } from './modules/local/plots_minimap2/main'
include { PLOTS_MINIMAP2 as PLOTS_VIRUS } from './modules/local/plots_minimap2/main'

include { SORT_INDEX_BAM } from './modules/local/bam_sort_index/main'
include { KRAKEN2_BUILD } from './modules/nf-core/kraken2/build/main'
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_COMBINEKREPORTS } from './modules/nf-core/krakentools/combinekreports/main'
include { KRAKEN2_COMBINE_REPORTS } from './modules/local/kraken2_combineReports/main'
include { KRONA_KREPORT2KRONA} from './modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_PLOTS } from './modules/nf-core/krona/ktimporttext/main'
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
        host_map = MINIMAPHOST(trimmed.reads,
                                params.host_fasta,
                                true, false, false, true, false)
        unmapped = host_map.unmapped
        mapping_logs = host_map.flagstat
    } else {
        error("Unsupported aligner. Options are 'bbmap' and 'minimap2'.")
    }

    unmapped_broadcast = unmapped.broadcast()

    // * --- Map to kingdoms ---
    def selected_kingdoms = params.kingdoms.split(',').collect { it.trim().toLowerCase() }
    def flagstat_channels = []

    // FUNGI
    if (selected_kingdoms.contains('fungi')) {
        fungi_map = MINIMAPFUNGI(unmapped_broadcast, 
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
        fungi_annotated_paf = MINIMAPFUNGI_ANNOTATE_PAF(fungi_grouped_paf_ch)
        fungi_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { fungi_collected_annotated_paths }
        fungi_concatenated_pafs = MINIMAPFUNGI_CONCAT_PAFS(fungi_collected_annotated_paths)
        // plots
        PLOTS_FUNGI(fungi_concatenated_pafs.cat, params.mapq_cutoff, "fungi", params.fungi_annotation)
    }
    // BACTERIA
    if (selected_kingdoms.contains('bacteria')) {
        bacteria_map = MINIMAPBACTERIA(unmapped_broadcast,
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
        bacteria_annotated_paf = MINIMAPBACTERIA_ANNOTATE_PAF(bacteria_grouped_paf_ch)
        bacteria_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { bacteria_collected_annotated_paths }
        bacteria_concatenated_pafs = MINIMAPBACTERIA_CONCAT_PAFS(bacteria_collected_annotated_paths)
        // plots
        if (params.bacteria_annotation) {
            PLOTS_BACTERIA(bacteria_concatenated_pafs.cat, params.mapq_cutoff, "bacteria", params.bacteria_annotation)
        }
    }
    // VIRUS
    if (selected_kingdoms.contains('virus')) {
        // REFSEQ
        virus_refseq_map = MINIMAPVIRUS(unmapped_broadcast,
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
        virus_refseq_annotated_paf = MINIMAPVIRUS_ANNOTATE_PAF(virus_refseq_grouped_paf_ch)
        virus_refseq_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { virus_refseq_collected_annotated_paths }
        virus_refseq_concatenated_pafs = MINIMAPVIRUS_CONCAT_PAFS(virus_refseq_collected_annotated_paths)
        // plots
        if (params.virus_refseq_annotation){
            PLOTS_VIRUS(virus_refseq_concatenated_pafs.cat, params.mapq_cutoff, "virus_refseq", params.virus_refseq_annotation)
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
        protozoa_map = MINIMAPPROTOZOA(unmapped_broadcast,
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
        protozoa_annotated_paf = MINIMAPPROTOZOA_ANNOTATE_PAF(protozoa_grouped_paf_ch)
        protozoa_annotated_paf.paf
            .map { it[1] }
            .collect()
            .set { protozoa_collected_annotated_paths }
        protozoa_concatenated_pafs = MINIMAPPROTOZOA_CONCAT_PAFS(protozoa_collected_annotated_paths)
        // plots
        PLOTS_PROTOZOA(protozoa_concatenated_pafs.cat, params.mapq_cutoff, "protozoa", params.protozoa_annotation)
    }
    all_flagstats = Channel.empty().mix(*flagstat_channels)

    // --- Taxonomic classification with Kraken2 ---
    kraken2_db_name = params.kraken2_db.tokenize('/').last().replaceAll(/[^a-zA-Z0-9]/, '')
    kraken2_db_name_ch = Channel.value(kraken2_db_name)

    kraken = KRAKEN2_KRAKEN2(unmapped, params.kraken2_db, params.kraken2_confidence, true, true, kraken2_db_name_ch)
    kraken_logs = kraken.report

    // Krona plots from krakentools
    kreport = KRONA_KREPORT2KRONA(kraken.report)
    KRONA_PLOTS(kreport.txt)

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
    kraken_reports_combined = KRAKEN2_COMBINE_REPORTS(grouped_kraken_reports_ch)
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

    // --- Collect all reports for MultiQC ---
    collect_reports_input = fastqc.html
        .map { it[1] }
        .merge(
            fastqc.zip.map { it[1] },
            trim_logs.map { it[1] },
            mapping_logs.map { it[1] },
            kraken_logs.map { it[1] },
            all_flagstats.map { it[1] }
        )
        .collect()

    multiqc_input = MULTIQC_COLLECT_REPORTS(collect_reports_input)

    MULTIQC(multiqc_input)
}
