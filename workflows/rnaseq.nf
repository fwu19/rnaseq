#!/usr/bin/env nextflow 

include { GET_FASTQ_PATHS } from '../modules/get_fastq_paths'
include { CHECK_INPUT  } from '../modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc.nf'
include { CAT_FASTQ  } from '../modules/cat_fastq.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { XENOFILTER } from '../modules/xenofilter.nf'
include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { HS_METRICS  } from '../modules/hs_metrics'
include { MULTIQC  } from '../modules/multiqc.nf'
include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'
include { GENERATE_REPORT } from '../modules/generate_report.nf'
include { OUTPUT_PARAMS  } from '../modules/output_params.nf'
//include { TEST  } from './modules/test.nf'


ch_dummy_csv = Channel.fromPath("$projectDir/assets/dummy_file.csv", checkIfExists: true)

ch_metadata = params.metadata ? Channel.fromPath( params.metadata, checkIfExists: true ) : ch_dummy_csv

//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml")

workflow RNASEQ {
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

    /*
    * run STAR alignment 
    * for pdx workflow, align to both graft and host genomes 
    * for pdx workflow, run XenofilteR to remove reads with host origin
    */
    ch_counts = Channel.empty()
    ch_bam = Channel.empty()
    ch_star_log = Channel.empty()
    if (params.run_alignment){
        STAR(
            ch_reads, 
            params.genome, 
            params.star, 
            params.gtf
        )
        ch_bam = STAR.out.bam
        // [ [meta], path(bam) ]
        ch_star_log = STAR.out.log
        // [ [meta], path(log) ]
        ch_counts = STAR.out.counts
        // [ [meta], path("ReadsPerGene.tab") ]

        if (params.workflow == 'pdx'){
            STAR_HOST(
                ch_reads, 
                params.genome_host, 
                params.star_host, 
                params.gtf_host
            )
            ch_bam_host = STAR_HOST.out.bam
            // [ [meta], path(bam) ]
            ch_star_log_host = STAR_HOST.out.log
            // [ [meta], path("*") ]

            ch_bam
                .join (ch_bam_host)
                .set {ch_bam_paired}
            // ch_bam_paired.view()
            // [ [meta], [path/to/graft.{bam,bai}], [path/to/host.{bam,bai}]]

            XENOFILTER(
                ch_bam_paired, 
                params.genome, 
                params.mm_threshold
            )
            ch_bam = XENOFILTER.out.bam 
            // [ [meta], path/to/filtered.bam ]
        }
    }


    /*
    * run featureCounts
    */
    if (params.run_featurecounts){
        FEATURECOUNTS(
            ch_bam, 
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
    ch_fastqc = Channel.empty()
    ch_rnaseqc = Channel.empty()
    ch_rseqc = Channel.empty()
    ch_hs_metrics = Channel.empty()

    if (params.run_qc){
        /*
        * FastQC
        */
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
        if (params.workflow == 'exome' && params.run_hs_metrics){
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
    ch_multiqc = Channel.empty()
    if (params.run_multiqc){
        MULTIQC(
            ch_star_log.map{it[1]}.flatten().collect().ifEmpty([]), 
            ch_fastqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rseqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rnaseqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([])
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
            //ch_hs_metrics.collect{it[1]}.ifEmpty([]),
            ch_dp.ifEmpty([]),
            ch_report_rmd
        )
    }
}
