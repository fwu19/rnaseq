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
    samplesheet
    srcdir
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

    samplesheet
        .splitCsv( header: true )
        .map { row -> [ id: row.id, sample_group: row.sample_group ] }
        .set { ch_metadata }

    /* for STAR, map2genome.csv, map2transcriptome.csv, gene_counts.csv  */
    if (params.aligner == 'star'){
        // target/graft alignment
        star_dir = file("${srcdir}/${params.star_dir}/")
        bam_bai  = star_dir.exists() ? ch_metadata
                .map { 
                    meta -> [meta, meta.id, file("${meta.id}.bam"), file("${meta.id}.bam.bai") ]                    
                } : Channel.empty()
        WRITE_CSV_BAM(
            (params.run_alignment ? ch_bam_bai : bam_bai)
                .map {
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.star_dir}/${it[2].name}"] + [bai: "${params.star_dir}/${it[3].name}"]  
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.csv"        
        )

        tx_bam = star_dir.exists() ? ch_metadata
                .map { 
                    it -> [ it, it.id, file("${it.id}.Aligned.toTranscriptome.out.bam") ]
                } : Channel.empty()
        WRITE_CSV_TX_BAM(
            (params.run_alignment ? ch_tx_bam : tx_bam)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [tx_bam: "${params.star_dir}/_work/${it[0].id}/${it[2].name}"]
                }
                .collect(),
            "map2transcriptome.${params.aligner}.${params.genome}.csv"        
        )
        
        star_counts = star_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.ReadsPerGene.out.tab") ]
                } : Channel.empty()
        WRITE_CSV_STAR_COUNTS(
            (params.run_alignment ? ch_star_counts : star_counts)
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.star_dir}/_work/${it[0].id}/${it[2].name}" ] 
                    }
                    .collect(),
            "gene_counts.${params.aligner}.${params.genome}.csv"        
        )

        // host alignment
        star_host_dir = file("${srcdir}/${params.star_host_dir}/")
        bam_bai_host = star_host_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.bam"), file("${it.id}.bam.bai") ]
                } : Channel.empty()
        WRITE_CSV_BAM_HOST(
            (params.run_alignment ? ch_bam_bai_host : bam_bai_host)
                .map {
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.star_host_dir}/${it[2].name}"] + [bai: "${params.star_host_dir}/${it[3].name}"] 
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome_host}.csv"        
        )

        tx_bam_host = star_host_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.Aligned.toTranscriptome.out.bam") ]
                } : Channel.empty()
        WRITE_CSV_TX_BAM_HOST(
            (params.run_alignment ? ch_tx_bam_host : tx_bam_host)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [tx_bam: "${params.star_host_dir}/_work/${it[0].id}/${it[2].name}"]
                }
                .collect(),
            "map2transcriptome.${params.aligner}.${params.genome_host}.csv"        
        )

        star_counts_host = star_host_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.ReadsPerGene.out.tab") ]
                } : Channel.empty()
        WRITE_CSV_STAR_COUNTS_HOST(
            (params.run_alignment ? ch_star_counts_host : star_counts_host)
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.star_host_dir}/_work/${it[0].id}/${it[2].name}" ] 
                    }
                    .collect(),
            "gene_counts.${params.aligner}.${params.genome_host}.csv"        
        )

        // xenofilteR filtered alignment
        star_xeno_dir = file("${srcdir}/${params.star_xeno_dir}/")
        bam_bai_xeno = star_xeno_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.bam"), file("${it.id}.bam.bai") ]
                } : Channel.empty()
        WRITE_CSV_BAM_XENO(
            (params.run_alignment ? ch_bam_bai_xeno : bam_bai_xeno)
                .map {
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.star_xeno_dir}/${it[2].name}"] + [bai: "${params.star_xeno_dir}/${it[3].name}"] 
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.xenofilteR_filtered.csv"        
        )

    }

    /* for BWA, map2genome.csv */
    if (params.aligner == 'bwa-mem'){
        // target/graft alignment
        bwa_dir = file("${srcdir}/${params.bwa_dir}/")
        bam_bai  = bwa_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.bam"), file("${it.id}.bam.bai") ]
                } : Channel.empty()

        WRITE_CSV_BAM(
            (params.run_alignment ? ch_bam_bai : bam_bai)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.bwa_dir}/${it[2].name}"]
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.csv"        
        )

        // host alignment
        bwa_host_dir = file("${srcdir}/${params.bwa_host_dir}/")
        bam_bai_host  = bwa_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.bam"), file("${it.id}.bam.bai") ]
                } : Channel.empty()

        WRITE_CSV_BAM_HOST(
            (params.run_alignment ? ch_bam_bai_host : bam_bai_host)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.bwa_host_dir}/${it[2].name}"]
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome_host}.csv"        
        )

        // xenofilteR filtered alignment
        bwa_xeno_dir = file("${srcdir}/${params.bwa_xeno_dir}/")
        bam_bai_xeno  = bwa_dir.exists() ? ch_metadata
                .map { 
                    it -> [it, it.id, file("${it.id}.bam"), file("${it.id}.bam.bai") ]
                } : Channel.empty()

        WRITE_CSV_BAM_XENO(
            (params.run_alignment ? ch_bam_bai_xeno : bam_bai_xeno)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [bam: "${params.bwa_xeno_dir}/${it[2].name}"]
                }
                .collect(),
            "map2genome.${params.aligner}.${params.genome}.xenofilteR_filtered.csv"        
        )

    }

    /* gene_counts.csv */
    fc_counts_dir = file("${srcdir}/${params.featurecounts_dir}/")
    fc_counts = fc_counts_dir.exists() ? ch_metadata
        .map { 
            it -> [it, it.id, file("${it.id}.exonicFragmentCounts.txt") ]
        } : Channel.empty()

    WRITE_CSV_FC_COUNTS(
        (params.run_quant_genes && params.run_gene_count && params.run_featurecounts ? ch_fc_counts : fc_counts)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.featurecounts_dir}/${it[2].name}" ] 
                }
                .collect(),
        "gene_counts.featurecounts.${params.genome}.csv"        
    )

    fc_counts_host_dir = file("${srcdir}/${params.featurecounts_host_dir}/")
    fc_counts_host = fc_counts_host_dir.exists() ? ch_metadata
        .map { 
            it -> [it, it.id, file("${it.id}.exonicFragmentCounts.txt")]
        } : Channel.empty()

    WRITE_CSV_FC_COUNTS_HOST(
        (params.run_quant_genes && params.run_gene_count_host && params.run_featurecounts_host ? ch_fc_counts_host : fc_counts_host)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [gene_count: "${params.featurecounts_host_dir}/${it[2].name}" ] 
                }
                .collect(),
            "gene_counts.featurecounts.${params.genome_host}.csv"        
    )


    /* salmon.csv */
    salmon_dir = file("${srcdir}/${params.salmon_dir}/")
    salmon = salmon_dir.exists() ? ch_metadata
        .map { 
            it -> [it, it.id, file("${salmon_dir}/${it.id}/") ]
        } : Channel.empty()

    WRITE_CSV_SALMON(
        (params.run_quant_transcripts && params.run_tx_count && params.run_salmon ? ch_salmon : salmon)
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [salmon: "${params.salmon_dir}/${it[2].name}/"]
                    }
                .collect(),
        "salmon.${params.genome}.csv"        
    )
    
    salmon_host_dir = file("${srcdir}/${params.salmon_host_dir}/")
    salmon_host = salmon_host_dir.exists() ? ch_metadata
        .map { 
            it -> [it, it.id, file("${salmon_host_dir}/${it.id}/") ]
        } : Channel.empty()

    WRITE_CSV_SALMON_HOST(
        (params.run_quant_transcripts && params.run_tx_count_host && params.run_salmon_host ? ch_salmon_host : salmon_host)
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [salmon: "${params.salmon_host_dir}/${it[2].name}/"]
                    }
                .collect(),
        "salmon.${params.genome_host}.csv"        
    )


    /* arriba.csv */
    arriba_dir = file("${srcdir}/${params.arriba_dir}/")
    arriba = arriba_dir.exists() ? ch_metadata
        .map { 
            it -> [it, it.id, file("${arriba_dir}/${it.id}/") ]
        } : Channel.empty()

    WRITE_CSV_ARRIBA(
        (params.run_fusion && params.run_arriba ? ch_arriba : arriba)
                .map { 
                    it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [arriba: "${params.arriba_dir}/${it[2].name}/"]
                }
                .collect(),
        "arriba.${params.genome}.csv"        
    )


}
