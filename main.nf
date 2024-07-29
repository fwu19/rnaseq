#!/usr/bin/env nextflow 

//include { CHECK_INPUT  } from './modules/check_input.nf'
include { FASTQC  } from './modules/fastqc.nf'
include { STAR  } from './modules/star.nf'
//include { CAT_STAR_COUNTS  } from './modules/cat_star_counts.nf'
include { RNASEQC  } from './modules/rnaseqc.nf'
include { RSEQC  } from './modules/rseqc.nf'
include { MULTIQC  } from './modules/multiqc.nf'
include { OUTPUT_PARAMS  } from './modules/output_params.nf'
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
        .set {samples_ch}


workflow {
    //CHECK_INPUT(params.input)
    
    fastqc = FASTQC(samples_ch, params.outdir)
    
    star = STAR(samples_ch, params.star, params.gtf, params.outdir)
    
    //CAT_STAR_COUNTS()
    
    if (params.gtfQC){
        rnaseqc = RNASEQC(star, params.gtfQC, params.strand, params.readType, params.outdir)
    }else{
        rnaseqc = Channel.empty()
    }
    
    if (params.rseqcBed){
        rseqc = RSEQC(star, params.rseqcBed, params.txBed, params.geneBed, params.outdir)
    }else{
        rseqc = Channel.empty()
    }
    
    
    MULTIQC(star, fastqc, rseqc, rnaseqc, params.outdir)
    //TEST(fastqc_path, params.outdir)

    /* Report params used for the workflow */
    OUTPUT_PARAMS( params.outdir )
}
