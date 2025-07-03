#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include {COLLECT_REPORTS} from './modules/local/multiqc/main'
include { FASTQC } from './modules/nf-core/fastqc/main'                                                                                                                                     
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main'                                                                                                                           
include { BBMAP_BBDUK } from './modules/nf-core/bbmap/bbduk/main'                                                                                                                                                
include { FASTP } from './modules/nf-core/fastp/main'                                                                                                                                       
include { BBMAP_INDEX } from './modules/nf-core/bbmap/index/main'                                                                                                                                                
// include { BOWTIE2_BUILD } from './modules/nf-core/bowtie2/build/main'         
include { BBMAP_ALIGN } from './modules/nf-core/bbmap/align/main'                                                                                                                                                                                                            
// include { BOWTIE2_ALIGN } from './modules/nf-core/bowtie2/align/main'                                                                                                                       
include { MINIMAP2_ALIGN as MINIMAP_HOST } from './modules/nf-core/minimap2/align/main'                                                                                                                     
include {SORT_INDEX_BAM} from './modules/local/bam_sort_index/main'
include { KRAKEN2_BUILD } from './modules/nf-core/kraken2/build/main'                                                                                                                       
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'                                                                                                                   
include { KRAKEN2_PARSE } from './modules/local/kraken2_parse/main'                                                                                                                   
include { KRAKEN2_TAXONOMY as KRAKEN2_TAXONOMY_PRE} from './modules/local/kraken2_taxonomy/main'                                                                                                                   
include { KRAKEN2_TAXONOMY as KRAKEN2_TAXONOMY_POST} from './modules/local/kraken2_taxonomy/main'                                                                                                                   
include { COUNTS_MERGE as  COUNTS_MERGE_PRE} from './modules/local/counts_merge/main'                                                                                                                   
include { COUNTS_MERGE as  COUNTS_MERGE_POST} from './modules/local/counts_merge/main'                                                                                                                   
include { PLOTS_HEATMAP as PLOTS_HEATMAP_PRE} from './modules/local/plots_heatmap/main'                                                                                                                   
include { PLOTS_HEATMAP as PLOTS_HEATMAP_POST} from './modules/local/plots_heatmap/main'                                                                                                                   
include { MASH_SCREEN } from './modules/nf-core/mash/screen/main'                                                                                                                           
include { MASH_SKETCH } from './modules/nf-core/mash/sketch/main'                                                                                                                           
include { NCBIGENOMEDOWNLOAD } from './modules/nf-core/ncbigenomedownload/main'
include { MINIMAP2_ALIGN as MINIMAP_CANDIDATES } from './modules/nf-core/minimap2/align/main'                                                                                                                     


// params.input = "data/samples.csv"

workflow {
    // Create a new channel of metadata from the sample sheet passed to the pipeline through the --input parameter
    samples = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    // Update the metadata with the single_end parameter and put reads files in a list
    samples = samples.map( {
        row ->
            if (!row[2]) {
                return new Tuple (row[0] + [ single_end:true ], [ row[1] ])
            } else {
                return new Tuple (row[0] + [ single_end:false ], [ row[1], row[2] ])
            }
        })
	samples.view()

	
	// * fastqc * //
    fastqc = FASTQC(samples)

	// * adapter trimming * //
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


	// * host removal * //
	if (params.host_removal_tool == 'bbmap') {
		if (!params.index_host) {
			host_map = BBMAP_ALIGN(trimmed.reads, params.host_bbmap_index)
			unmapped = host_map.unmapped
			mapping_logs = host_map.stats
		} else {
			index = BBMAP_INDEX(params.host_fasta).index
			host_map = BBMAP_ALIGN(trimmed.reads, index)
			unmapped = host_map.unmapped
			mapping_logs = host_map.stats
		}
	// } else if (params.host_removal_tool == 'bowtie2') {
	// 	if (!params.index_host) {
	// 		unmapped = BOWTIE2_ALIGN(trimmed, 
	// 					params.host_bowtie2_index,
	// 					params.host_fasta,
	// 					true, false).fastq
	// 	} else {
	// 		index = BOWTIE2_BUILD(params.host_fasta).index
	// 		unmapped = BOWTIE2_ALIGN(trimmed, 
	// 					index,
	// 					params.host_fasta,
	// 					true, false).fastq
	// 	}		
	} else if (params.host_removal_tool == 'minimap2') {
		host_map = MINIMAP_HOST(trimmed.reads,
						params.host_fasta,
						true, false, false, true)
		unmapped = host_map.unmapped
		mapping_logs = host_map.log
	} else {
		error("Unsupported aligner. Options are 'bbmap' and 'minimap2'.")
	}
	
	// * ncbi db preparation * //
	// ...

	// * taxonomic classification * //
	if (params.taxonomy_tool == 'all') {
		kraken2_db_name = params.kraken2_db.tokenize('/').last()
		kraken2_db_name_ch = Channel.value(kraken2_db_name)

		kraken = KRAKEN2_KRAKEN2(unmapped, params.kraken2_db, true, true, kraken2_db_name_ch)
		kraken_logs = kraken.report
		
		kraken_parse = KRAKEN2_PARSE(kraken.classified_reads_assignment, params.kraken2_parse_threshold, kraken2_db_name_ch)
		
		kraken_taxo_pre = KRAKEN2_TAXONOMY_PRE(kraken_parse.stat_pre, params.kraken2_parse_threshold, kraken2_db_name_ch, file(params.taxonomy_db))
		kraken_taxo_post = KRAKEN2_TAXONOMY_POST(kraken_parse.stat_post, params.kraken2_parse_threshold, kraken2_db_name_ch, file(params.taxonomy_db))

		taxonomy_files_pre_ch = kraken_taxo_pre.taxonomy.map { it -> it[1] }
		counts_pre = COUNTS_MERGE_PRE(taxonomy_files_pre_ch.collect())
		taxonomy_files_post_ch = kraken_taxo_post.taxonomy.map { it -> it[1] }
		counts_post = COUNTS_MERGE_POST(taxonomy_files_pre_ch.collect())

		PLOTS_HEATMAP_PRE(counts_pre.counts)
		PLOTS_HEATMAP_POST(counts_post.counts)

		mash_db_name = params.mash_screen_db.tokenize('/').last()
		mash_db_name_ch = Channel.value(mash_db_name)
		
		MASH_SCREEN(unmapped, params.mash_screen_db, mash_db_name_ch)
	} else if (params.taxonomy_tool == 'kraken2') {
		kraken2_db_name = params.kraken2_db.tokenize('/').last()
		kraken2_db_name_ch = Channel.value(kraken2_db_name)
		if (!params.ncbi_db_prepare) {
			kraken = KRAKEN2_KRAKEN2(unmapped, params.kraken2_db, true, true, kraken2_db_name_ch)
			kraken_logs = kraken.report
			kraken_parse = KRAKEN2_PARSE(kraken.classified_reads_assignment, params.kraken2_parse_threshold, kraken2_db_name_ch)
			kraken_taxo = KRAKEN2_TAXONOMY(kraken_parse.stat_post, params.kraken2_parse_threshold, kraken2_db_name_ch, file(params.taxonomy_db))
		} else {
			// *** untested *** //
			db = KRAKEN2_BUILD(params.ncbi_db, true).db
			kraken = KRAKEN2_KRAKEN2(unmapped, db, true, true, kraken2_db_name_ch)
			kraken_parse = KRAKEN2_PARSE(kraken.classified_reads_assignment, params.kraken2_parse_threshold, kraken2_db_name_ch)
			kraken_taxo = KRAKEN2_TAXONOMY(kraken_parse.stat_post, params.kraken2_parse_threshold, kraken2_db_name_ch, file(params.taxonomy_db))
			kraken_logs = kraken.report
		}
	} else if (params.taxonomy_tool == 'mash') {
		mash_db_name = params.mash_screen_db.tokenize('/').last()
		mash_db_name_ch = Channel.value(mash_db_name)
		if (!params.ncbi_db_prepare) {
			MASH_SCREEN(unmapped, params.mash_screen_db, mash_db_name_ch)
		} else {
			db = MASH_SKETCH(params.ncbi_db).mash
			MASH_SCREEN(unmapped, db, mash_db_name_ch)
		}
	} else {
		error("Unsupported taxonomy tool. Options are 'kraken2', 'mash', and 'all'.")
	}

	// // Channel with all taxonomy files
	// taxonomy_files_ch = Channel.fromPath('./kraken2/*_taxonomy.tsv')
	// COUNTS_MERGE(taxonomy_files_ch)


	if (params.minimap2_candidates) {
		//* Download genomes of interest * //
		// NCBIGENOMEDOWNLOAD(candidates)
		ncbi = NCBIGENOMEDOWNLOAD(params.candidates, params.genomes_filename)

		//* Mapping against candidates * //
		map_candidates = MINIMAP_CANDIDATES(unmapped,
						ncbi.cat_fna,
						true,
						false,
						false,
						false)
		SORT_INDEX_BAM(map_candidates.bam)
	}

    // -- COLLECT REPORTS FOR MULTIQC --
    if (params.taxonomy_tool == 'mash') {
			collect_reports_input = fastqc.html
		.map { it[1] }
		.merge(
			fastqc.zip.map { it[1] },
			trim_logs.map { it[1] },
			mapping_logs.map { it[1] },
		)
		.collect()
	} else {
			collect_reports_input = fastqc.html
		.map { it[1] }
		.merge(
			fastqc.zip.map { it[1] },
			trim_logs.map { it[1] },
			mapping_logs.map { it[1] },
			kraken_logs.map { it[1] }
		)
		.collect()
	}

    multiqc_input = COLLECT_REPORTS(collect_reports_input)

    MULTIQC(multiqc_input)
    
}

