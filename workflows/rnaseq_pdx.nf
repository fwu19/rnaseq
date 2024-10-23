#!/usr/bin/env nextflow 

//include { CHECK_INPUT  } from './modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc.nf'
include { CAT_FASTQ  } from '../modules/cat_fastq.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { MULTIQC  } from '../modules/multiqc.nf'
include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { OUTPUT_PARAMS  } from '../modules/output_params.nf'
//include { TEST  } from './modules/test.nf'


Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .set {ch_samples} // all fields in sample sheet


workflow RNASEQ_PDX {
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
        STAR(ch_reads, params.genome, params.star, params.gtf)
        ch_bam_graft = STAR.out.bam
        ch_star_graft = STAR.out.star

        STAR_HOST(ch_reads, params.genome_host, params.star_host, params.gtf_host)
        ch_bam_host = STAR_HOST.out.bam
        ch_star_host = STAR_HOST.out.star

        ch_bam_graft
            .join (ch_bam_host)
            .set {ch_bam_paired}
        // ch_bam_paired.view()
        // [sample_id, [path/to/graft.{bam,bai}], [path/to/host.{bam,bai}]]
        XENOFILTER(ch_bam_paired, params.genome, params.mm_threshold)
        ch_bam = XENOFILTER.out.bam // [ sample_id, path/to/filtered.bam ]
    }


    /*
    * run featureCounts
    */
    if (params.run_featurecounts){
        FEATURECOUNTS(ch_bam, params.gtf)
        ch_counts = FEATURECOUNTS.out.counts // [ sample_id, path/to/counts.txt ]
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
    * Report params used for the workflow 
    */
    OUTPUT_PARAMS( params.outdir )

}
