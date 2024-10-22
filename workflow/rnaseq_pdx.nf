#!/usr/bin/env nextflow 

//include { CHECK_INPUT  } from './modules/check_input.nf'
include { FASTQC  } from './modules/fastqc.nf'
include { CAT_FASTQ  } from './modules/cat_fastq.nf'
include { STAR } from './modules/star.nf'
include { STAR  as  STAR_HOST } from './modules/star.nf'
include { XENOFILTER } from './modules/xenofilter.nf'
include { RNASEQC  } from './modules/rnaseqc.nf'
include { RSEQC  } from './modules/rseqc.nf'
include { MULTIQC  } from './modules/multiqc.nf'
include { FEATURECOUNTS } from './modules/featureCounts.nf'
include { OUTPUT_PARAMS  } from './modules/output_params.nf'
//include { TEST  } from './modules/test.nf'


Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .set {ch_samples} // all fields in sample sheet


workflow {
    //CHECK_INPUT(params.input)

    ch_samples
        .map {
            row -> [row.sample_id, [row.fastq_1, row.fastq_2]]
        }
        .groupTuple (by: [0])
        .map {
            sample_id,fastq -> 
                [sample_id, fastq.flatten()]
        }
        .set {ch_fastq}
    //ch_fastq.view() 
    // [ sample_id, path/to/*.fastq.gz ]

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
    * run STAR alignment to graft and host genomes
    * run XenofilteR to remove reads with host origin
    */
    if (params.run_alignment){
        STAR(ch_reads, params.genome)
        ch_bam_graft = STAR.out.bam
        ch_star_graft = STAR.out.star

        STAR_HOST(ch_reads, params.genome_host)
        ch_bam_host = STAR_HOST.out.bam
        ch_star_host = STAR_HOST.out.star

        ch_bam_graft
            .cross(ch_bam_host)
            .set(ch_bam_paired)
        XENOFILTER(ch_bam_paired)
        ch_bam = XENOFILTER.out.bam // [ sample_id, path/to/filtered.bam ]
    }


    /*
    * run featureCounts
    */
    if (params.run_alignment & params.run_featurecounts){
        FEATURECOUNTS(ch_bam, params.gtf)
        ch_counts = FEATURECOUNTS.out.counts // [ sample_id, path/to/counts.txt ]
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
            row -> [row.sample_id, [row.fastq_1, row.fastq_2]]
        }
        .map {
            sample_id,fastq -> 
                [sample_id, fastq.flatten()]
        }
        .set {ch_fastqc}

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
            ch_star.collect(), 
            ch_fastqc.collect(), 
            ch_rseqc.collect(), 
            ch_rnaseqc.collect()
        )
    
    }
    
    /*
    * MultiQC
    */
    if (params.run_de){
        DIFFERENTIAL_EXPRESSION()

    }
    /* 
    * Report params used for the workflow 
    */
    OUTPUT_PARAMS( params.outdir )

}
