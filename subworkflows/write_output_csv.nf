nextflow.enable.dsl=2

/*
* Writes out csv files containing output paths
*/

include { WRITE_CSV as WRITE_CSV_ALIGN_FASTQ } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_FC_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_STAR_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SALMON } from '../modules/write_csv.nf'

workflow WRITE_OUTPUT_CSV {
    take:
    ch_bam_bai
    ch_bam_bai_host
    ch_bam_bai_xeno
    ch_star_counts
    ch_fc_counts
    ch_salmon

    main: 

    /* regular or exome workflow, mapped.csv */
    if (params.gene_level && params.run_alignment && params.workflow != 'pdx' ){
        def subdir = params.aligner == 'star' ? "${params.aligner.toUpperCase()}/${params.genome}" :
            params.aligner == 'bwa-mem' ? "${params.aligner}/${params.genome}" :
            null

        WRITE_CSV_ALIGN_FASTQ(
            ch_bam_bai
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${subdir}/${it[2].name}"] + [bai: "${subdir}/${it[3].name}"]
                }
                .collect(),
            "mapped.${params.aligner}.csv"        
        )
        
    }

    /* pdx workflow, mapped.csv */
    if(params.gene_level && params.run_alignment && params.workflow == 'pdx'){
        def subdir = params.aligner == 'star' ? "${params.aligner.toUpperCase()}" :
            params.aligner == 'bwa-mem' ? "${params.aligner}" :
            null

        WRITE_CSV_ALIGN_FASTQ(
            ch_bam_bai
                .join(ch_bam_bai_host)
                .join(ch_bam_bai_xeno)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [graft_bam: "${subdir}/${params.genome}/_unfiltered/${it[2].name}" ] + [graft_bai: "${subdir}/${params.genome}/_unfiltered/${it[3].name}" ] + [host_bam: "${subdir}/${params.host_genome}/_unfiltered/${it[4].name}" ] + [host_bai: "${subdir}/${params.host_genome}/_unfiltered/${it[5].name}" ] + [bam: "${subdir}/${params.genome}/${it[6].name}" ] + [bai: "${subdir}/${params.genome}/${it[7].name}" ]
                    }
                .collect(),
            "mapped.${params.aligner}.csv"        
        )
    }     

    /* gene_counts.csv */
    if (params.gene_level && params.run_gene_count && params.run_featurecounts ){
        def subdir = "featureCounts/${params.genome}" 
        WRITE_CSV_FC_COUNTS(
            ch_fc_counts
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${subdir}/${it[2].name}" ] 
                }
                .collect(),
            "gene_counts.featurecounts.csv"        
        )

    }
    
    if (params.gene_level && params.run_alignment && params.workflow != 'pdx' && params.aligner == 'star'){
        def subdir = "${params.aligner.toUpperCase()}/${params.genome}" 
        WRITE_CSV_STAR_COUNTS(
            ch_star_counts
                .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${subdir}/_work/${it[0].id}/${it[2].name}" ] 
                }
                .collect(),
            "gene_counts.${params.aligner}.csv"        
        )

    }

    /* salmon.csv */
    if (params.transcript_level && params.run_salmon ){
        WRITE_CSV_SALMON(
                ch_salmon
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [salmon: "salmon/${params.genome}/${it[2].name}/"]
                }
                .collect(),
        "salmon.csv"        
        )
    }

}
