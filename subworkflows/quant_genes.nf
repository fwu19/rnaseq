/*
* Analyze gene-level expression 
*/

include { WRITE_CSV as WRITE_CSV_FEATURECOUNTS} from '../modules/write_csv.nf'
include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { GENERATE_GENE_COUNT_MATRIX } from '../modules/generate_gene_count_matrix.nf'
include { DIFFERENTIAL_GENES } from '../modules/differential_genes.nf'

workflow QUANT_GENES {
    take:
    samplesheet
    ch_bam
    ch_bai
    ch_counts
    gene_txt
    tx_bed
    experiment


    main: 
    ch_expr = Channel.empty()
    ch_de = Channel.empty()
    ch_versions  = Channel.empty()



    /*
    * Generate read count matrix
    */
    if (params.run_gene_count){
        
        if (params.run_featurecounts){
                FEATURECOUNTS(
                ch_bam, 
                params.gtf,
                experiment
                )     
                ch_counts = FEATURECOUNTS.out.counts
                // [ [meta], val(out_prefix), path("count.txt") ]
                ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)

                WRITE_CSV_FEATURECOUNTS(
                    ch_counts
                    .map { 
                        it -> [id: it[0].id] + [sample_group: it[0].sample_group] + [ gene_count: "featureCounts/${params.genome}/${it[2].name}" ] 
                    }
                    .collect(),
                "gene_counts.featureCounts.csv"        
                )
            
        }

        GENERATE_GENE_COUNT_MATRIX(
            samplesheet, 
            ch_counts.map{it[2]}.collect(),
            gene_txt,
            params.length_col, // default: gene_length
            experiment
        )
        ch_expr = GENERATE_GENE_COUNT_MATRIX.out.rds
        ch_versions = ch_versions.mix(GENERATE_GENE_COUNT_MATRIX.out.versions)
    }else{
        if (params.run_de){
            ch_expr = Channel.fromPath("${params.outdir}/expression_quantification/all_samples.gene_raw_counts.txt")
        }
        
    }

    /*
    * differential genes
    */
    if (params.run_de){
       DIFFERENTIAL_GENES(
            samplesheet, 
            file(params.comparison),
            ch_expr, 
            "gene_length",
            gene_txt
        )
        ch_de = DIFFERENTIAL_GENES.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_GENES.out.versions)
    }

    
    emit:
    expr = ch_expr
    counts = ch_counts
    de = ch_de
    versions = ch_versions

}
