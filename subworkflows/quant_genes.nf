/*
* Analyze gene-level expression 
*/

include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { GENERATE_GENE_COUNT_MATRIX } from '../modules/generate_gene_count_matrix.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'

workflow QUANT_GENES {
    take:
    samplesheet
    ch_bam
    ch_counts
    gene_txt

    main: 
    ch_gene_rds = Channel.empty()
    ch_de = Channel.empty()

    /*
    * Parse saved alignments
    */
    if (params.step in [ "expression_quanification", "differential_expression"] ){     
            align_csv = file("${params.outdir}/csv/align_fastq.csv", checkIfExists: true)       
            if (params.workflow == 'pdx'){
                align_csv
                    .splitCsv(header: true)
                    .map {
                        row -> [ row, row.id, row.filtered_graft_bam ]
                    }
                    .set { ch_bam}
            }else{
                align_csv
                    .splitCsv(header: true)
                    .map {
                        row -> [ row, row.id, row.bam ]
                    }
                    .set { ch_bam}

            }

            if (params.workflow != 'pdx' && params.aligner == 'star'){
                align_csv
                    .splitCsv(header: true)
                    .map {
                        row -> [ row, row.id, row.gene_count ]
                    }
                    .set { ch_counts}
            }
            
    }

    /*
    * run featureCounts if needed
    */
    if (params.run_featurecounts){
            FEATURECOUNTS(
                ch_bam, 
                params.gtf,
                params.read_type,
                params.strand
            )
        
            ch_counts = FEATURECOUNTS.out.counts
            // [ [meta], val(out_prefix), path("count.txt") ]
    }


    /*
    * Generate read count matrix
    */
    GENERATE_GENE_COUNT_MATRIX(
            samplesheet, 
            ch_counts.map{it[2]}.flatten().collect(),
            gene_txt,
            params.length_col?:"gene_length",
            params.strand,
            params.workflow
    )
    ch_gene_rds = GENERATE_GENE_COUNT_MATRIX.out.rds
    
    
    /*
    * differential genes
    */
    if (params.run_de && params.comparison){
        DIFFERENTIAL_EXPRESSION(
            samplesheet, 
            file(params.comparison, checkIfExists: true),
            ch_gene_rds, 
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2,
            params.gene_txt ?: gene_txt,
            params.de_gene_type ?: 'all'
        )
        ch_de = DIFFERENTIAL_EXPRESSION.out.rds
    }

    
    emit:
    gene_rds = ch_gene_rds
    de = ch_de
    counts = ch_counts
    gene_txt = gene_txt

}
