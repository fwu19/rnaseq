/*
* Analyze gene-level expression 
*/

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
    strand
    read_type

    main: 
    ch_gene_rds = Channel.empty()
    ch_de = Channel.empty()
    ch_versions  = Channel.empty()

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
                read_type,
                strand
            )
        
            ch_counts = FEATURECOUNTS.out.counts
            // [ [meta], val(out_prefix), path("count.txt") ]
            ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions)
    }


    /*
    * Generate read count matrix
    */
    GENERATE_GENE_COUNT_MATRIX(
            samplesheet, 
            ch_counts.map{it[2]}.flatten().collect(),
            gene_txt,
            params.length_col, // default: gene_length
            strand,
            params.workflow
    )
    ch_gene_rds = GENERATE_GENE_COUNT_MATRIX.out.rds
    ch_versions = ch_versions.mix(GENERATE_GENE_COUNT_MATRIX.out.versions)
    
    /*
    * differential genes
    */
    if (params.run_de && params.comparison){
        DIFFERENTIAL_GENES(
            samplesheet, 
            file(params.comparison, checkIfExists: true),
            ch_gene_rds, 
            "gene_length",
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2,
            gene_txt,
            params.de_gene_type
        )
        ch_de = DIFFERENTIAL_GENES.out.rds
        ch_versions = ch_versions.mix(DIFFERENTIAL_GENES.out.versions)
    }

    
    emit:
    gene_rds = ch_gene_rds
    de = ch_de
    counts = ch_counts
    gene_txt = gene_txt
    versions = ch_versions

}
