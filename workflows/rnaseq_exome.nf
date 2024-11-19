#!/usr/bin/env nextflow 

include { CHECK_INPUT  } from '../modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc'
include { CAT_FASTQ  } from '../modules/cat_fastq'
include { STAR  } from '../modules/star'
include { RNASEQC  } from '../modules/rnaseqc'
include { RSEQC  } from '../modules/rseqc'
include { HS_METRICS  } from '../modules/hs_metrics'
include { MULTIQC  } from '../modules/multiqc'
include { FEATURECOUNTS } from '../modules/featureCounts'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression'
include { OUTPUT_PARAMS  } from '../modules/output_params'
//include { TEST  } from './modules/test.nf'

ch_dummy_csv = Channel.fromPath("$projectDir/assets/dummy_file.csv", checkIfExists: true)

if(!params.input){ exit 1, "Provide params.input!" }else{ ch_input = Channel.fromPath( params.input, checkIfExists:true ) }
ch_metadata = params.metadata ? Channel.fromPath( params.metadata, checkIfExists: true ) : Channel.empty()



workflow RNASEQ_EXOME {
    /*
    * Run input check
    */
    samplesheet = ch_input
    if (params.run_input_check){
        CHECK_INPUT(
            ch_input,
            ch_metadata
        )
        samplesheet = CHECK_INPUT.out.csv

        CHECK_INPUT.out.fq
        .splitCsv(header: true)
        .map {
            row -> [ row.id, row, row.fastq_1, row.fastq_2]
        }
        .groupTuple (by: [0])
        .set {ch_fastq}
        //ch_fastq.view() 
        // [ sample_id, [paths/to/R1.fastq.gz], [paths/to/R2.fastq.gz] ]
    }

    /*
    * cat fastq files if required
    */
    if (params.cat_fastq){
        CAT_FASTQ(
            ch_fastq
        )
        ch_reads = CAT_FASTQ.out.reads

    }else{
        ch_reads = ch_fastq
    }

    /*
    * run STAR alignment
    */
    ch_bam = Channel.empty()
    ch_log = Channel.empty()
    if (params.run_alignment){
        STAR(
            ch_reads, 
            params.genome, 
            params.star, 
            params.gtf
        )
        ch_bam = STAR.out.bam 
        // [ [meta], path(bam) ]
        ch_log = STAR.out.log
        // [ [meta], path(log) ]
    }

    /*
    * run featureCounts
    */
    ch_counts = Channel.empty()
    if (params.run_featurecounts){
        FEATURECOUNTS(
            ch_bam, 
            params.gtf
        )
        ch_counts = FEATURECOUNTS.out.counts
        // [ [meta], path("count.txt") ]
    }

    /*
    * differential expression
    */
    ch_comparison = params.comparison ? Channel.fromPath(params.comparison, checkIfExists: true) : Channel.empty()
    if (params.run_de){
        DIFFERENTIAL_EXPRESSION(
            samplesheet, 
            ch_comparison, 
            ch_counts.map{it[1]}.collect(), 
            params.gene_txt
        )

    }

    /*
    * run QC
    */    
    if (params.run_qc){
        /*
        * FastQC
        */
        ch_fastqc = Channel.empty()
        if (params.run_fastqc){
            FASTQC(
                CHECK_INPUT.out.fq
                    .splitCsv(header: true)
                    .map {
                        row -> [ row, row.fastq_1, row.fastq_2]
                    }
            )
            ch_fastqc = FASTQC.out.qc
        }

        /*
        * RNASeQC
        */
        ch_rnaseqc = Channel.empty()
        if (params.run_rnaseqc){
            RNASEQC(
                ch_bam,
                params.rnaseqc_gtf,
                params.strand,
                params.read_type
            )
            ch_rnaseqc = RNASEQC.out.qc
            // [ [meta], path("*") ]
        }

        /*
        * RSeQC
        */
        ch_rseqc = Channel.empty()
        if (params.run_rseqc){
            RSEQC(
                ch_bam,
                params.rseqc_bed,
                params.tx_bed,
                params.gene_bed
            )
            ch_rseqc = RSEQC.out.qc
            // [ [meta], path("*") ]
        }

        /*
        * GATK collect_hs_metrics
        */
        ch_hs_metrics = Channel.empty()
        if (params.run_hs_metrics){
            HS_METRICS(
                ch_bam,
                params.genome_fa,
                params.target_region
            )
            ch_hs_metrics = HS_METRICS.out.qc
            // [ [meta], path("*") ]
        }

    }


    /*
    * MultiQC
    */
    if (params.run_multiqc){
        MULTIQC(
            ch_log.collect(), 
            ch_fastqc.map{it[1]}.flatten().collect(), 
            ch_rseqc.map{it[1]}.flatten().collect(), 
            ch_rnaseqc.map{it[1]}.flatten().collect(),
            ch_hs_metrics.map{it[1]}.flatten().collect()
        )
    
    }
    
}
