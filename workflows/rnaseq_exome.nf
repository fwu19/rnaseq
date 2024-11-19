#!/usr/bin/env nextflow 

//include { CHECK_INPUT  } from './modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc'
include { CAT_FASTQ  } from '../modules/cat_fastq'
include { STAR  } from '../modules/star'
//include { CAT_STAR_COUNTS  } from '../modules/cat_star_counts.nf'
include { RNASEQC  } from '../modules/rnaseqc'
include { RSEQC  } from '../modules/rseqc'
include { HS_METRICS  } from '../modules/hs_metrics'
include { MULTIQC  } from '../modules/multiqc'
include { FEATURECOUNTS } from '../modules/featureCounts'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression'
include { OUTPUT_PARAMS  } from '../modules/output_params'
//include { TEST  } from './modules/test.nf'

if(!params.input){ exit 1, "Provide params.input!" }
Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .set {ch_samples}


workflow RNASEQ_EXOME {
    //CHECK_INPUT(params.input)

    ch_samples
        .map {
            row -> [row.sample_id, row.fastq_1, row.fastq_2]
        }
        .groupTuple (by: [0])
        .set {ch_fastq}
    //ch_fastq.view() 
    // [ sample_id, [paths/to/R1.fastq.gz], [paths/to/R2.fastq.gz] ]

    /*
    * cat fastq files if required
    */
    if (params.cat_fastq){
        CAT_FASTQ(ch_fastq)
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
        STAR(ch_reads, params.genome, params.star, params.gtf)
        ch_bam = STAR.out.bam
        ch_log = STAR.out.log
    }

    /*
    * run featureCounts
    */
    ch_counts = Channel.empty()
    if (params.run_alignment & params.run_featurecounts){
        FEATURECOUNTS(ch_bam, params.gtf)
        ch_counts = FEATURECOUNTS.out.counts.collect()
        // [path/to/all/counts.txt]
    }

    /*
    * differential expression
    */
    if (params.run_de){
        DIFFERENTIAL_EXPRESSION(params.input, params.comparison, ch_counts, params.gene_txt)

    }

    /*
    * run QC
    */    
    if (params.run_qc){
        /*
        * FastQC
        */
        fastqc = Channel.empty()
        if (params.run_fastqc){
            FASTQC(
                    ch_samples
                        .map { it -> [it.sample_id, it.fastq_1, it.fastq_2]}

            )
            ch_fastqc = FASTQC.out.qc
        }

        /*
        * RNASeQC
        */
        rnaseqc = Channel.empty()
        if (params.run_rnaseqc){
            RNASEQC(ch_bam)
            ch_rnaseqc = RNASEQC.out.qc
        }

        /*
        * RSeQC
        */
        rseqc = Channel.empty()
        if (params.run_rseqc){
            RSEQC(ch_bam)
            ch_rseqc = RSEQC.out.qc
        }

        /*
        * GATK collect_hs_metrics
        */
        hs_metrics = Channel.empty()
        if (params.run_hs_metrics){
            HS_METRICS(
                ch_bam,
                params.genome_fa,
                params.target_region
            )
            ch_hs_metrics = HS_METRICS.out.qc
        }

    }


    /*
    * MultiQC
    */
    if (params.run_multiqc){
        MULTIQC(
            ch_log.collect(), 
            ch_fastqc.flatten().collect(), 
            ch_rseqc.flatten().collect(), 
            ch_rnaseqc.flatten().collect(),
            ch_hs_metrics.flatten().collect()
        )
    
    }
    
}
