#!/usr/bin/env nextflow

include { MULTIQC } from './modules/nf-core/multiqc/main'
include { FASTQC } from './modules/nf-core/fastqc/main'                                                                                                                                     
include { TRIMMOMATIC } from './modules/nf-core/trimmomatic/main'                                                                                                                           
include { BBMAP_BBDUK } from './modules/nf-core/bbmap/bbduk/main'                                                                                                                                                
include { FASTP } from './modules/nf-core/fastp/main'                                                                                                                                       
include { BBMAP_INDEX } from './modules/nf-core/bbmap/index/main'                                                                                                                                                
include { BOWTIE2_BUILD } from './modules/nf-core/bowtie2/build/main'         
include { BBMAP_ALIGN } from './modules/nf-core/bbmap/align/main'                                                                                                                                                                                                            
include { BOWTIE2_ALIGN } from './modules/nf-core/bowtie2/align/main'                                                                                                                       
include { MINIMAP2_ALIGN } from './modules/nf-core/minimap2/align/main'                                                                                                                     
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'                                                                                                                   
include { MASH_SCREEN } from './modules/nf-core/mash/screen/main'                                                                                                                           
include { MASH_SKETCH } from './modules/nf-core/mash/sketch/main'                                                                                                                           
include { samplesheetToList } from 'plugin/nf-schema'


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
    FASTQC(samples)

	// * adapter trimming * //
	if (params.trim_tool == "bbduk") {
		trimmed = BBMAP_BBDUK(samples, params.adapters).reads
	} else if (params.trim_tool == "fastp") {
	    trimmed = FASTP(samples, params.adapters, false, false, false).reads
	} else if (params.trim_tool == "trimmomatic") {
		trimmed = TRIMMOMATIC(samples).trimmed_reads
	} else {
		error("Unsupported trim tool. Options are 'bbduk', 'fastp' or 'trimmomatic'.")
	}


	// * host removal * //
	if (params.host_removal_tool == 'bbmap') {
		if (!params.index_host) {
			unmapped = BBMAP_ALIGN(trimmed, params.host_bbmap_index).unmapped
			unmapped.view()
		} else {
			index = BBMAP_INDEX(params.host_fasta).index
			index.view()
			unmapped = BBMAP_ALIGN(trimmed, index).unmapped
		}
	} else if (params.host_removal_tool == 'bowtie2') {
		if (!params.index_host) {
			unmapped = BOWTIE2_ALIGN(trimmed, 
						params.host_bowtie2_index,
						params.host_fasta,
						true, true).fastq
		} else {
			index = BOWTIE2_BUILD(params.host_fasta).index
			unmapped = BOWTIE2_ALIGN(trimmed.reads, 
						index,
						params.host_fasta,
						true, true).fastq
		}		
	} else if (params.host_removal_tool == 'minimap2') {
		MINIMAP2_ALIGN(trimmed.reads,
						params.host_fasta,
						true, "bai", false, false)
	} else {
		error("Unsupported aligner. Options are 'bbmap', 'bowtie2' and 'minimap2'.")
	}

	// * ncbi db preparation * //
	// * taxonomic classification * //
	if (params.taxonomy_tool == 'all') {
		KRAKEN2_KRAKEN2(unmapped, params.kraken2_db, true, true)
		MASH_SCREEN(unmapped, params.mash_screen_db)
	} else if (params.taxonomy_tool == 'kraken2'){
		KRAKEN2_KRAKEN2(unmapped, params.kraken2_db, true, true)
	} else if (params.taxonomy_tool == 'mash'){
		if (!params.ncbi_db_prepare){
			MASH_SCREEN(unmapped, params.mash_screen_db)
		} else {
			db = MASH_SKETCH(params.ncbi_db).mash
			MASH_SCREEN(unmapped, db)
		}
	} else {
		error("Unsupported taxonomy tool. Options are 'kraken2', 'mash', and 'all'.")
	}
}
