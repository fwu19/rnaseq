#!/usr/bin/env nextflow 

//include { CHECK_INPUT  } from './modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc.nf'
include { CAT_FASTQ  } from '../modules/cat_fastq.nf'
include { STAR  } from '../modules/star.nf'
//include { CAT_STAR_COUNTS  } from '../modules/cat_star_counts.nf'
include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { MULTIQC  } from '../modules/multiqc.nf'
include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'
include { OUTPUT_PARAMS  } from '../modules/output_params.nf'
//include { TEST  } from './modules/test.nf'

/*
how to parse output to input of the next process
how to implement dependency?
*/

Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        // row is a list object
        //.view { row -> "${row.sample_id}, ${row.fastq_1}, ${row.fastq_2}" }
        //.take (3)
        .set {ch_samples}


workflow RNASEQ_REGULAR {
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
    if (params.run_alignment){
        STAR(ch_reads, params.genome, params.star, params.gtf)
        ch_bam = STAR.out.bam
        ch_star = STAR.out.star
    }

    /*
    * run featureCounts
    */
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
        ch_samples
        .map {
            row -> [row.sample_id, row.fastq_1, row.fastq_2]
        }
        .set {ch_fastqc}
        //ch_fastqc.view()

        fastqc = Channel.empty()
        if (params.run_fastqc){
            FASTQC(ch_fastqc)
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
    }


    /*
    * MultiQC
    */
    if (params.run_multiqc){
        MULTIQC(
            ch_star.flatten().collect(), 
            ch_fastqc.flatten().collect(), 
            ch_rseqc.flatten().collect(), 
            ch_rnaseqc.flatten().collect()
        )
    
    }
    
    /* 
    * Report params used for the workflow 
    */
    OUTPUT_PARAMS( params.outdir )

}
