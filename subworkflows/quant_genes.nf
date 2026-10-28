/*
* Analyze gene-level expression 
*/

include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { GENERATE_GENE_COUNT_MATRIX } from '../modules/generate_gene_count_matrix.nf'
include { DIFFERENTIAL_GENES } from '../modules/differential_genes.nf'
include { FEATURECOUNTS as FEATURECOUNTS_HOST } from '../modules/featureCounts.nf'
include { GENERATE_GENE_COUNT_MATRIX as GENERATE_GENE_COUNT_MATRIX_HOST } from '../modules/generate_gene_count_matrix.nf'
include { DIFFERENTIAL_GENES as DIFFERENTIAL_GENES_HOST } from '../modules/differential_genes.nf'

workflow QUANT_GENES {
    take:
    samplesheet
    ch_bam_bai
    ch_star_counts
    gene_txt
    ch_bam_bai_host
    ch_star_counts_host
    gene_txt_host
    srcdir
    infer_experiment

    main: 
    ch_fc_counts = Channel.empty()
    ch_fc_summary = Channel.empty()
    ch_expr = Channel.empty()
    ch_de = Channel.empty()
    de_csv = Channel.empty()
    ch_fc_counts_host = Channel.empty()
    ch_fc_summary_host = Channel.empty()
    ch_expr_host = Channel.empty()
    ch_de_host = Channel.empty()
    de_csv_host = Channel.empty()
    ch_versions  = Channel.empty()

    /*
    * Generate read count matrix
    */
    if (params.run_gene_count){
        if (params.run_featurecounts){
            if (!params.run_alignment){
                ch_bam_bai = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.${params.genome}.csv")
                    .splitCsv(header: true)
                    .map { it -> [ [ it ], it.id , "${srcdir}/${it.bam}", "${srcdir}/${it.bai}" ] }
            }

            FEATURECOUNTS(
                ch_bam_bai, 
                params.gtf,
                params.run_infer_experiment ? infer_experiment : file("${srcdir}/csv/infer_experiment.csv")
            )     
            ch_fc_counts = FEATURECOUNTS.out.counts // [ [meta], val(out_prefix), path("count.txt") ]
            ch_fc_summary = FEATURECOUNTS.out.summary
            ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)            
        }else{      
            if (!params.run_alignment){      
                ch_star_counts = Channel.fromPath("${srcdir}/csv/gene_counts.${params.aligner}.${params.genome}.csv")
                    .splitCsv(header: true)
                    .map { it -> [ [ it ], it.id , "${srcdir}/${it.gene_count}" ] }
            }
        }

        GENERATE_GENE_COUNT_MATRIX(
            samplesheet, 
            (params.run_featurecounts ? ch_fc_counts : ch_star_counts).map{it[2]}.collect(),
            gene_txt,
            params.length_col, // default: gene_length
            params.run_infer_experiment ? infer_experiment : file("${srcdir}/csv/infer_experiment.csv")
        )
        ch_expr = GENERATE_GENE_COUNT_MATRIX.out.rds
        ch_versions = ch_versions.mix(GENERATE_GENE_COUNT_MATRIX.out.versions)
    }

    /*
    * differential genes
    */
    if (params.run_de){
        DIFFERENTIAL_GENES(
            samplesheet, 
            params.comparison ? file(params.comparison) : file("$projectDir/assets/dummy_file.csv"),
            params.run_gene_count ? ch_expr : Channel.fromPath("${srcdir}/${params.gene_expr_dir}/all_samples.gene_raw_counts.txt"), 
            "gene_length",
            gene_txt,
            params.de_method
        )
        de_csv = DIFFERENTIAL_GENES.out.csv
        ch_de = DIFFERENTIAL_GENES.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_GENES.out.versions)
    }

    /*
    * Generate read count matrix
    */
    if (params.workflow == 'pdx' && params.run_gene_count_host){
        if (params.run_featurecounts_host){
            if (!params.run_alignment){
                ch_bam_bai_host = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.${params.genome_host}.csv")
                    .splitCsv(header: true)
                    .map { it -> [ [ it ], it.id , "${srcdir}/${it.bam}", "${srcdir}/${it.bai}" ] }
            }

            FEATURECOUNTS_HOST(
                ch_bam_bai_host, 
                params.gtf_host,
                params.run_infer_experiment ? infer_experiment : file("${srcdir}/csv/infer_experiment.csv")
            )     
            ch_fc_counts_host = FEATURECOUNTS_HOST.out.counts // [ [meta], val(out_prefix), path("count.txt") ]
            ch_fc_summary_host = FEATURECOUNTS_HOST.out.summary
            ch_versions = ch_versions.mix(FEATURECOUNTS_HOST.out.versions)            
        }else{      
            if (!params.run_alignment){      
                ch_star_counts_host = Channel.fromPath("${srcdir}/csv/gene_counts.${params.aligner}.${params.genome_host}.csv")
                    .splitCsv(header: true)
                    .map { it -> [ [ it ], it.id , "${srcdir}/${it.gene_count}" ] }
            }
        }

        GENERATE_GENE_COUNT_MATRIX_HOST(
            samplesheet, 
            (params.run_featurecounts_host ? ch_fc_counts_host : ch_star_counts_host).map{it[2]}.collect(),
            gene_txt_host,
            params.length_col, // default: gene_length
            params.run_infer_experiment ? infer_experiment : file("${srcdir}/csv/infer_experiment.csv")
        )
        ch_expr_host = GENERATE_GENE_COUNT_MATRIX_HOST.out.rds
        ch_versions = ch_versions.mix(GENERATE_GENE_COUNT_MATRIX_HOST.out.versions)
    }

    /*
    * differential genes
    */
    if (params.workflow == 'pdx' && params.run_de_host){
        DIFFERENTIAL_GENES_HOST(
            samplesheet, 
            params.comparison ? file(params.comparison) : file("$projectDir/assets/dummy_file.csv"),
            params.run_gene_count_host ? ch_expr_host : Channel.fromPath("${srcdir}/${params.gene_expr_host_dir}/all_samples.gene_raw_counts.txt"), 
            "gene_length",
            gene_txt_host,
            params.de_method
        )
        de_csv_host = DIFFERENTIAL_GENES_HOST.out.csv
        ch_de_host = DIFFERENTIAL_GENES_HOST.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_GENES_HOST.out.versions)
    }

    
    emit:
    expr = ch_expr
    fc_counts = ch_fc_counts
    fc_summary = ch_fc_summary
    de = ch_de
    de_csv
    expr_host = ch_expr_host
    fc_counts_host = ch_fc_counts_host
    fc_summary_host = ch_fc_summary_host
    de_host = ch_de_host
    de_csv_host
    versions = ch_versions

}
