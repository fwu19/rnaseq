nextflow.enable.dsl=2

/*
* Writes out csv files containing output paths
*/

include { WRITE_CSV as WRITE_CSV_BAM } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_BAM_HOST } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_BAM_XENO } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_FC_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_FC_COUNTS_HOST } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_STAR_COUNTS } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_STAR_COUNTS_HOST } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_TX_BAM } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_TX_BAM_HOST } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SALMON } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SALMON_HOST } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_ARRIBA } from '../modules/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_ARRIBA_HOST } from '../modules/write_csv.nf'

workflow WRITE_OUTPUT_CSV {
    take:
    ch_bam_bai
    ch_bam_bai_host
    ch_bam_bai_xeno
    ch_tx_bam
    ch_tx_bam_host
    ch_star_counts
    ch_star_counts_host
    ch_fc_counts
    ch_fc_counts_host
    ch_salmon
    ch_salmon_host
    ch_arriba

    main:     

    /* for STAR, map2genome.csv, map2transcriptome.csv, gene_counts.csv  */
    if (params.run_alignment && params.aligner == 'star'){
        WRITE_CSV_BAM(
            ch_bam_bai
                .map {
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.star_dir}/${it[2].name}"] + [bai: "${params.star_dir}/${it[3].name}"]  
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.csv"        
        )
            
        WRITE_CSV_TX_BAM(
            ch_tx_bam
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [tx_bam: "${params.star_dir}/_work/${it[0].id}/${it[2].name}"]
                }
                .collect(),
            "map2transcriptome.${params.aligner}.${params.genome}.csv"        
        )

        WRITE_CSV_STAR_COUNTS(
            ch_star_counts
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.star_dir}/_work/${it[0].id}/${it[2].name}" ] 
                    }
                    .collect(),
            "gene_counts.${params.aligner}.${params.genome}.csv"        
        )
            
        WRITE_CSV_BAM_HOST(
            ch_bam_bai_host
                .map {
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.star_host_dir}/${it[2].name}"] + [bai: "${params.star_host_dir}/${it[3].name}"] 
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome_host}.csv"        
        )

        WRITE_CSV_TX_BAM_HOST(
            ch_tx_bam_host
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [tx_bam: "${params.star_host_dir}/_work/${it[0].id}/${it[2].name}"]
                }
                .collect(),
            "map2transcriptome.${params.aligner}.${params.genome_host}.csv"        
        )

        WRITE_CSV_STAR_COUNTS_HOST(
            ch_star_counts_host
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.star_host_dir}/_work/${it[0].id}/${it[2].name}" ] 
                    }
                    .collect(),
            "gene_counts.${params.aligner}.${params.genome_host}.csv"        
        )

        WRITE_CSV_BAM_XENO(
            ch_bam_bai_xeno
                .map {
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.star_xeno_dir}/${it[2].name}"] + [bai: "${params.star_xeno_dir}/${it[3].name}"] 
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.xenofilteR_filtered.csv"        
        )

    }

    /* for BWA, map2genome.csv */
    if (params.run_alignment && params.aligner == 'bwa-mem'){
        WRITE_CSV_BAM(
            ch_bam_bai
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.bwa_dir}/${it[2].name}"]
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.csv"        
        )
    
        WRITE_CSV_BAM_HOST(
            ch_bam_bai_host
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.bwa_host_dir}/${it[2].name}"]
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome_host}.csv"        
        )

        WRITE_CSV_BAM_XENO(
            ch_bam_bai_xeno
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.bwa_xeno_dir}/${it[2].name}"]
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.xenofilteR_filtered.csv"        
        )

    }

    /* gene_counts.csv */
    if (params.run_quant_genes && params.run_gene_count && params.run_featurecounts ){
        WRITE_CSV_FC_COUNTS(
            ch_fc_counts
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.featurecounts_dir}/${it[2].name}" ] 
                }
                .collect(),
            "gene_counts.featurecounts.${params.genome}.csv"        
            )
    }

    if (params.run_quant_genes && params.run_gene_count_host && params.run_featurecounts_host ){
        WRITE_CSV_FC_COUNTS_HOST(
            ch_fc_counts_host
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.featurecounts_host_dir}/${it[2].name}" ] 
                }
                .collect(),
            "gene_counts.featurecounts.${params.genome_host}.csv"        
            )

    }


    /* salmon.csv */
    if (params.run_quant_transcripts && params.run_tx_count && params.run_salmon){
        WRITE_CSV_SALMON(
            ch_salmon
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [salmon: "${params.salmon_dir}/${it[2].name}/"]
                    }
                .collect(),
            "salmon.${params.genome}.csv"        
        )
    }

    if (params.run_quant_transcripts && params.run_tx_count_host && params.run_salmon_host){    
        WRITE_CSV_SALMON_HOST(
            ch_salmon_host
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [salmon: "${params.salmon_host_dir}/${it[2].name}/"]
                    }
                .collect(),
            "salmon.${params.genome_host}.csv"        
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
            "arriba.${params.genome}.csv"        
        )
    }


}
