#!/usr/bin/env nextflow 

/*
* SUBWORKFLOWS
*/
include { SPLIT_READS } from '../subworkflows/split_reads'

/*
* MODULES
*/
include { GET_FASTQ_PATHS } from '../modules/get_fastq_paths'
include { CHECK_INPUT  } from '../modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED  } from '../modules/fastqc.nf'
include { CUTADAPT  } from '../modules/cutadapt.nf'
include { CAT_FASTQ  } from '../modules/cat_fastq.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { SAMTOOLS_VIEW} from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HOST } from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_XENO } from '../modules/samtools_view.nf'
include { HS_METRICS  } from '../modules/hs_metrics'
include { MULTIQC  } from '../modules/multiqc.nf'
include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'
include { GENERATE_REPORT } from '../modules/generate_report.nf'
include { OUTPUT_PARAMS  } from '../modules/output_params.nf'


ch_dummy_csv = Channel.fromPath("$projectDir/assets/dummy_file.csv", checkIfExists: true)

ch_metadata = params.metadata ? Channel.fromPath( params.metadata, checkIfExists: true ) : ch_dummy_csv

//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml")

workflow RNASEQ_PDX {
    /*
    * Run input check
    */
    
    if (params.run_input_check){
        if ( params.input_dir =~ 'dummy' ){
            if ( params.input =~ 'dummy' ){
                exit 1, 'Neither --input nor --input_dir is specified!'
            }else {
                ch_input = Channel.fromPath( params.input, checkIfExists: true )
            }
        }else {
            GET_FASTQ_PATHS (
                Channel.fromPath("${params.input_dir}/", type: 'dir', checkIfExists: true)
            )
            ch_input = GET_FASTQ_PATHS.out.csv
        }

        CHECK_INPUT(
            ch_input,
            ch_metadata
        )
        samplesheet = CHECK_INPUT.out.csv
        // final samplesheet, one row per id
        //CHECK_INPUT.out.fq.view()
        // samplesheet, one row per pair of fastq
    }

    /*
    * cat fastq files if needed
    */
    if (params.run_cat_fastq){
        CAT_FASTQ(
            CHECK_INPUT.out.fq
                .splitCsv(header: true)
                .map {
                    row -> [ row.id, row.fastq_1, row.fastq_2 ]
                }
                .groupTuple (by: [0])
        )
        samplesheet
            .splitCsv( header: true )
            .map {
                row -> [ row.id, row ]
            }
            .join ( CAT_FASTQ.out.reads )
            .map { it -> [ it[1], it[2], it[3] ]}
            .set { ch_reads }
        
    }else{
        CHECK_INPUT.out.fq
                .splitCsv(header: true)
                .map {
                    row -> [ row, row.fastq_1, row.fastq_2 ]
                }
                .set { ch_reads}
    }
    // ch_reads.view()
    // [ [meta], path("R1.fastq.gz"), path("R2.fastq.gz") ]
    
    /*
    * FastQC on raw fastq files
    */
    ch_fastqc = Channel.empty()
    if (params.run_qc && params.run_fastqc){
        FASTQC(
            ch_reads
        )
        ch_fastqc = FASTQC.out.qc
    }

    /*
    * run cutadapt if needed
    */
    ch_fastqc_trimmed = Channel.empty()
    if (params.run_cut_adapt){
        CUTADAPT(
            ch_reads
        )
        ch_reads = CUTADAPT.out.fq

        if (params.run_qc && params.run_fastqc){
            FASTQC_TRIMMED(
                ch_reads
            )
            ch_fastqc_trimmed = FASTQC.out.qc
        }

    }
    // ch_reads.view()
    // [ [meta], path("R1.fastq.gz"), path("R2.fastq.gz") ]
    
    /*
    * run STAR alignment 
    * for pdx workflow, align to both graft and host genomes 
    * for pdx workflow, run XenofilteR to remove reads with host origin
    */
    ch_bam = Channel.empty()
    if (params.run_alignment){
        STAR(
            ch_reads
                .map{ it -> [ it[0], it[0].id, it[1], it[2] ]}, 
            params.genome, 
            params.star, 
            params.gtf
        )
        ch_bam = STAR.out.bam
        ch_star_log = STAR.out.log
        // [ [meta], val(out_prefix), path(bam) ]

        STAR_HOST(
            ch_reads
                .map{ it -> [ it[0], it[0].id, it[1], it[2] ]}, 
            params.genome_host, 
            params.star_host, 
            params.gtf_host
        )
        ch_bam_host = STAR_HOST.out.bam
        ch_star_host_log = STAR_HOST.out.log
        // [ [meta], val(out_prefix), path(bam) ]

        ch_bam_xeno = Channel.empty()
        if (params.run_split_fastq){
            SPLIT_READS(
                ch_reads
            )
            ch_bam_xeno = SPLIT_READS.out.bam_xeno

        }else{
            ch_bam
                .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                .join (
                    ch_bam_host
                        .map{ it -> [ [ it[0], it[1] ], it[2] ]}
                )
                .map{ it -> [ it[0][0], it[0][1], it[1], it[2] ] }
                .set {ch_bam_paired}
            // ch_bam_paired.view()
            // [ [meta], val(out_prefix), [path/to/graft.{bam,bai}], [path/to/host.{bam,bai}]]

            XENOFILTER(
                ch_bam_paired, 
                params.genome, 
                params.mm_threshold
            )
            ch_bam_xeno = XENOFILTER.out.bam 
            // ch_bam_xeno.view()
            // [ [meta], val(out_prefix), path("*.{bam,bai}") ]      

        }

    }


    /*
    * run featureCounts
    */
    ch_counts = Channel.empty()
    if (params.run_featurecounts){
        FEATURECOUNTS(
            ch_bam_xeno, 
            params.gtf,
            params.read_type,
            params.strand
        )
        ch_counts = FEATURECOUNTS.out.counts
        // [ [meta], path("count.txt") ]
    }

    /*
    * differential expression
    */
    ch_dp = Channel.empty()
    if (params.run_de){
        DIFFERENTIAL_EXPRESSION(
            samplesheet, 
            Channel.fromPath(params.comparison, checkIfExists: true), 
            ch_counts.map{it[1]}.collect().ifEmpty([]), 
            params.gene_txt,
            params.length_col,
            params.strand,
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2
        )
        ch_dp = DIFFERENTIAL_EXPRESSION.out.rds
    }

    /*
    * run QC
    */    
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()
    ch_bam_stat = Channel.empty()
    ch_bam_host_stat = Channel.empty()
    ch_bam_xeno_stat = Channel.empty()
    if (params.run_qc){
        /*
        * RNASeQC
        */
        if (params.run_rnaseqc){
            RNASEQC(
                ch_bam_xeno,
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
        if (params.run_rseqc){
            RSEQC(
                ch_bam_xeno,
                params.rseqc_bed,
                params.tx_bed,
                params.gene_bed
            )
            ch_rseqc = RSEQC.out.qc
            // [ [meta], path("*") ]
        }

        /*
        * samtools
        */
        SAMTOOLS_VIEW(
            ch_bam
        )
        ch_bam_stat = SAMTOOLS_VIEW.out.data

        SAMTOOLS_VIEW_HOST(
            ch_bam_host
        )
        ch_bam_host_stat = SAMTOOLS_VIEW_HOST.out.data

        SAMTOOLS_VIEW_XENO(
            ch_bam_xeno
        )
        ch_bam_xeno_stat = SAMTOOLS_VIEW_XENO.out.data
    }


    /*
    * MultiQC
    */
    ch_multiqc = Channel.empty()
    if (params.run_multiqc){
        MULTIQC(
            ch_star_log.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_star_host_log.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_fastqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_fastqc_trimmed.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_rseqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rnaseqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_bam_stat.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_bam_host_stat.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_bam_xeno_stat.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_counts.map{it[1]}.collect().ifEmpty([]),
        )
        ch_multiqc = MULTIQC.out.data
        //ch_multiqc.view()
    }

    /*
    * Generate a report
    */
    if (params.run_report){
        ch_report_rmd = params.local_assets ? Channel.fromPath("${params.local_assets}/report/", type: 'dir', checkIfExists: true) : Channel.fromPath("$projectDir/assets/report/", type: 'dir', checkIfExists: true)
        GENERATE_REPORT(
            params.workflow,
            samplesheet,
            ch_multiqc.ifEmpty([]),
            ch_dp.ifEmpty([]),
            ch_report_rmd
        )
    }
}