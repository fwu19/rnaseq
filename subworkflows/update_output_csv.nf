nextflow.enable.dsl=2

/*
* Writes out csv files containing output paths
*/

include { WRITE_CSV as WRITE_CSV_ALIGN_FASTQ } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_FC_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_STAR_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SALMON } from '../modules/write_csv.nf'

workflow UPDATE_OUTPUT_CSV {
    take:
    samplesheet

    main: 

    samplesheet
        .splitCsv( header: true )
        .map { row -> [ id: row.id, sample_group: row.sample_group ] }
        .set { ch_metadata }

    /* regular or exome workflow, mapped.csv */
    if (params.workflow != 'pdx' ){
        def subdir = params.aligner == 'star' ? "${params.aligner.toUpperCase()}/${params.genome}" :
            params.aligner == 'bwa-mem' ? "${params.aligner}/${params.genome}" :
            null

        WRITE_CSV_ALIGN_FASTQ(
            ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [bam: "${subdir}/${it.id}.bam"] + [bai: "${subdir}/${it.id}.bam.bai"]
                }
                .collect(),
            "mapped.${params.aligner}.csv"        
        )
        
    }

    /* pdx workflow, mapped.csv */
    if(params.workflow == 'pdx'){
        def subdir = params.aligner == 'star' ? "${params.aligner.toUpperCase()}" :
            params.aligner == 'bwa-mem' ? "${params.aligner}" :
            null

        WRITE_CSV_ALIGN_FASTQ(
            ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [graft_bam: "${subdir}/${params.genome}/_unfiltered/${it.id}.bam" ] + [graft_bai: "${subdir}/${params.genome}/_unfiltered/${it.id}.bam.bai" ] + [host_bam: "${subdir}/${params.host_genome}/_unfiltered/${it.id}.bam" ] + [host_bai: "${subdir}/${params.host_genome}/_unfiltered/${it.id}.bam.bai" ] + [bam: "${subdir}/${params.genome}/${it.id}.bam" ] + [bai: "${subdir}/${params.genome}/${it.id}.bam.bai" ]
                    }
                .collect(),
            "mapped.${params.aligner}.csv"        
        )
    }     

    /* gene_counts.csv */
    if (params.run_featurecounts ){
        def subdir = "featureCounts/${params.genome}" 
        WRITE_CSV_FC_COUNTS(
            ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [gene_count: "${subdir}/${it.id}.exonicFragmentCounts.txt" ] 
                }
                .collect(),
            "gene_counts.featurecounts.csv"        
        )

    }
    
    if (params.workflow != 'pdx' && params.aligner == 'star'){
        def subdir = "${params.aligner.toUpperCase()}/${params.genome}" 
        WRITE_CSV_STAR_COUNTS(
            ch_metadata
                .map { 
                        it -> [id: it.id] + [sample_group: it.sample_group] + [gene_count: "${subdir}/_work/${it.id}/${it.id}.ReadsPerGene.out.tab" ] 
                }
                .collect(),
            "gene_counts.${params.aligner}.csv"        
        )

    }

    /* salmon.csv */
    if (params.transcript_level && params.run_salmon ){
        WRITE_CSV_SALMON(
            ch_metadata
                .map { 
                    it -> [id: it.id] + [sample_group: it.sample_group] + [salmon: "salmon/${params.genome}/${it.id}/"]
                }
                .collect(),
        "salmon.csv"        
        )
    }

}
