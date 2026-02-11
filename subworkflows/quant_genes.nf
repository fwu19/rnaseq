/*
* Analyze gene-level expression 
*/

include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { GENERATE_GENE_COUNT_MATRIX } from '../modules/generate_gene_count_matrix.nf'
include { DIFFERENTIAL_GENES } from '../modules/differential_genes.nf'

workflow QUANT_GENES {
    take:
    samplesheet
    ch_bam_bai
    ch_star_counts
    gene_txt
    infer_experiment
    srcdir

    main: 
    ch_fc_counts = Channel.empty()
    ch_fc_summary = Channel.empty()
    ch_expr = Channel.empty()
    ch_de = Channel.empty()
    de_csv = Channel.empty()
    ch_versions  = Channel.empty()

    /*
    * Generate read count matrix
    */
    if (params.run_gene_count){
        if (!params.run_alignment){
            def my_dir = new File("${srcdir}")
            def srcdir = my_dir.absolutePath
            ch_bam_bai = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.csv")
                .splitCsv(header: true)
                .map { it -> [ [ it ], it.id , "${srcdir}/${it.bam}", "${srcdir}/${it.bai}" ] }
        }

        if (params.run_featurecounts){
            FEATURECOUNTS(
                    ch_bam_bai, 
                    params.gtf,
                    params.run_infer_experiment ? infer_experiment : file("${srcdir}/csv/infer_experiment.csv")
            )     
            ch_fc_counts = FEATURECOUNTS.out.counts // [ [meta], val(out_prefix), path("count.txt") ]
            ch_fc_summary = FEATURECOUNTS.out.summary
            ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)            
        }else{            
            ch_star_counts = Channel.fromPath("${srcdir}/csv/gene_counts.${params.aligner}.csv")
                    .splitCsv(header: true)
                    .map { it -> [ [ it ], it.id , "${srcdir}/${it.gene_count}" ] }
        }

        ch_counts = params.run_featurecounts ? ch_fc_counts : ch_star_counts
        GENERATE_GENE_COUNT_MATRIX(
            samplesheet, 
            ch_counts.map{it[2]}.collect(),
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
            gene_txt
        )
        de_csv = DIFFERENTIAL_GENES.out.csv
        ch_de = DIFFERENTIAL_GENES.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_GENES.out.versions)
    }

    
    emit:
    expr = ch_expr
    fc_counts = ch_fc_counts
    fc_summary = ch_fc_summary
    de = ch_de
    de_csv
    versions = ch_versions

}
