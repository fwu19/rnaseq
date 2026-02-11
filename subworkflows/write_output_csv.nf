nextflow.enable.dsl=2

/*
* Writes out csv files containing output paths
*/

include { WRITE_CSV as WRITE_CSV_ALIGN_FASTQ } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_FC_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_STAR_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MAP_TX } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SALMON } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_ARRIBA } from '../modules/write_csv.nf'

workflow WRITE_OUTPUT_CSV {
    take:
    ch_bam_bai
    ch_bam_bai_host
    ch_bam_bai_xeno
    ch_tx_bam
    ch_star_counts
    ch_fc_counts
    ch_salmon
    ch_arriba

    main: 

    /* regular or exome workflow, map2genome.csv and map2transcriptome.csv */
    if (params.run_alignment && params.workflow != 'pdx' ){
        if (params.aligner == 'star'){
            WRITE_CSV_ALIGN_FASTQ(
                ch_bam_bai
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.star_dir}/${it[2].name}"] + [bai: "${params.star_dir}/${it[3].name}"]
                }
                .collect(),
                "map2genome.${params.aligner}.csv"        
            )
        
            WRITE_CSV_MAP_TX(
                ch_tx_bam
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [tx_bam: "${params.star_dir}/_work/${it[0].id}/${it[2].name}"]
                }
                .collect(),
                "map2transcriptome.${params.aligner}.csv"        
            )
        }

        if (params.aligner == 'bwa-mem'){
            WRITE_CSV_ALIGN_FASTQ(
                ch_bam_bai
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.bwa_dir}/${it[2].name}"] + [bai: "${params.bwa_dir}/${it[3].name}"]
                }
                .collect(),
                "map2genome.${params.aligner}.csv"        
            )
        }

    }

    /* pdx workflow, map2genome.csv */
    if(params.run_alignment && params.workflow == 'pdx'){
        if (params.aligner == 'star'){
            WRITE_CSV_ALIGN_FASTQ(
                ch_bam_bai
                .join(ch_bam_bai_host)
                .join(ch_bam_bai_xeno)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [graft_bam: "${params.star_dir}/${it[2].name}" ] + [graft_bai: "${params.star_dir}/${it[3].name}" ] + [host_bam: "${params.star_host_dir}/${it[5].name}" ] + [host_bai: "${params.star_host_dir}/${it[6].name}" ] + [bam: "${params.star_xeno_dir}/${it[8].name}" ] + [bai: "${params.star_xeno_dir}/${it[9].name}" ]
                    }
                .collect(),
                "map2genome.${params.aligner}.csv"        
            )
        }

        if (params.aligner == 'bwa-mem'){
            WRITE_CSV_ALIGN_FASTQ(
                ch_bam_bai
                .join(ch_bam_bai_host)
                .join(ch_bam_bai_xeno)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [graft_bam: "${params.star_dir}/${it[2].name}" ] + [graft_bai: "${params.star_dir}/${it[3].name}" ] + [host_bam: "${params.star_host_dir}/${it[5].name}" ] + [host_bai: "${params.star_host_dir}/${it[6].name}" ] + [bam: "${params.star_xeno_dir}/${it[8].name}" ] + [bai: "${params.star_xeno_dir}/${it[9].name}" ]
                    }
                .collect(),
                "map2genome.${params.aligner}.csv"        
            )
        }
    }     

    /* gene_counts.csv */
    if (params.run_quant_genes && params.run_gene_count && params.run_featurecounts ){
        WRITE_CSV_FC_COUNTS(
            ch_fc_counts
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.featurecounts_dir}/${it[2].name}" ] 
                }
                .collect(),
            "gene_counts.featurecounts.csv"        
        )
        
    }
    
    if (params.run_alignment && params.workflow != 'pdx' && params.aligner == 'star'){
        WRITE_CSV_STAR_COUNTS(
            ch_star_counts
                .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.star_dir}/_work/${it[0].id}/${it[2].name}" ] 
                }
                .collect(),
            "gene_counts.${params.aligner}.csv"        
        )

    }

    /* salmon.csv */
    if (params.run_quant_transcripts && params.run_salmon && params.workflow != 'pdx'){
        WRITE_CSV_SALMON(
            ch_salmon
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [salmon: "${params.salmon_dir}/${it[2].name}/"]
                }
                .collect(),
            "salmon.csv"        
        )
    }

    /* arriba.csv */
    if (params.run_fusion && params.run_arriba ){
        WRITE_CSV_ARRIBA(
            ch_arriba
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [arriba: "${params.arriba_dir}/${it[2].name}/"]
                }
                .collect(),
            "arriba.csv"        
        )
    }
}
