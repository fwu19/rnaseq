#!/usr/bin/env nextflow 

include { GET_FASTQ_PATHS } from '../modules/get_fastq_paths'
include { CHECK_INPUT  } from '../modules/check_input.nf'
include { FASTQC  } from '../modules/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED  } from '../modules/fastqc.nf'
include { CUTADAPT  } from '../modules/cutadapt.nf'
include { CAT_FASTQ  } from '../modules/cat_fastq.nf'
include { STAR } from '../modules/star.nf'
include { STAR  as  STAR_HOST } from '../modules/star.nf'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq.nf'
include { SALMON  } from '../modules/salmon.nf'
include { ARRIBA  } from '../modules/arriba.nf'
include { RNASEQC  } from '../modules/rnaseqc.nf'
include { RSEQC  } from '../modules/rseqc.nf'
include { HS_METRICS  } from '../modules/hs_metrics'
include { SAMTOOLS_VIEW} from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HOST } from '../modules/samtools_view.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_XENO } from '../modules/samtools_view.nf'
include { MULTIQC  } from '../modules/multiqc.nf'
include { MULTIQC_PDX  } from '../modules/multiqc_pdx.nf'
include { FEATURECOUNTS } from '../modules/featureCounts.nf'
include { DIFFERENTIAL_EXPRESSION } from '../modules/differential_expression.nf'
include { DIFFERENTIAL_TRANSCRIPTS } from '../modules/differential_transcripts.nf'
include { GENERATE_GENE_COUNT_MATRIX } from '../modules/generate_gene_count_matrix.nf'
include { GENERATE_TRANSCRIPT_COUNT_MATRIX } from '../modules/generate_transcript_count_matrix.nf'
include { GENERATE_REPORT } from '../modules/generate_report.nf'
//include { OUTPUT_PARAMS  } from '../modules/output_params.nf'
//include { TEST  } from './modules/test.nf'

include { SINGLE_LIB } from '../subworkflows/single_lib.nf'
include { SPLIT_LIB } from '../subworkflows/split_lib.nf'

ch_dummy_csv = Channel.fromPath("$projectDir/assets/dummy_file.csv", checkIfExists: true)

ch_metadata = params.metadata ? Channel.fromPath( params.metadata, checkIfExists: true ) : ch_dummy_csv

ch_multiqc_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/multiqc_config.yml")

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
                Channel.fromPath("${params.input_dir}", checkIfExists: true)
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
    * run cutadapt if needed
    */
    ch_fastqc_trimmed = Channel.empty()
    ch_cutadapt_js = Channel.empty()
    if (params.run_cut_adapt){
        CUTADAPT(
            ch_reads
        )
        ch_reads = CUTADAPT.out.fq
        ch_cutadapt_js = CUTADAPT.out.js

        if (params.run_qc && params.run_fastqc){
            FASTQC_TRIMMED(
                ch_reads
            )
            ch_fastqc_trimmed = FASTQC_TRIMMED.out.qc
        }

    }
    // ch_reads.view()
    // [ [meta], path("R1.fastq.gz"), path("R2.fastq.gz") ]

    /*
    * run STAR alignment 
    * for pdx workflow, align to both graft and host genomes 
    * for pdx workflow, run XenofilteR to remove reads with host origin
    */
    ch_counts = Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_star_log = Channel.empty()
    if (params.run_alignment){
        STAR(
            ch_reads
                .map{ it -> [ it[0], it[0].id, it[1], it[2] ]}, 
            params.genome, 
            params.star, 
            params.gtf
        )
        ch_bam = STAR.out.bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_bai = STAR.out.bai
        // [ [meta], val(out_prefix), path(bai) ]
        ch_tx_bam = STAR.out.tx_bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_star_log = STAR.out.log
        // [ [meta], val(out_prefix), path(log) ]
        ch_counts = STAR.out.counts
        // [ [meta], val(out_prefix), path("ReadsPerGene.tab") ]
    }

    if (params.run_alignment & params.run_salmon){
        SALMON(
            ch_tx_bam,
            params.tx_fa
        )
        ch_salmon = SALMON.out.sf
        // [ [meta], val(out_prefix), path("${out_prefix}/") ]
    }

    if (params.run_alignment & params.run_arriba){
        ARRIBA(
            ch_reads
                .map{ it -> [ it[0], it[0].id, it[1], it[2] ]}, 
            params.genome, 
            params.star, 
            params.gtf,
            params.genome_fa,
            params.blacklist,
            params.known_fusions,
            params.protein_domains
        )
        
    }

    if (params.run_alignment & params.workflow == 'pdx'){
        /*
        * align to host genome
        */
        ch_bam_host = Channel.empty()
        ch_star_host_log = Channel.empty()

        STAR_HOST(
            ch_reads
                .map{ it -> [ it[0], it[0].id, it[1], it[2] ]}, 
            params.genome_host, 
            params.star_host, 
            params.gtf_host
        )
        ch_bam_host = STAR_HOST.out.bam
        // [ [meta], val(out_prefix), path(bam) ]
        ch_star_log_host = STAR_HOST.out.log
        // [ [meta], val(out_prefix), path(bam) ]

        /*
        * filter out reads of host origin
        */
        ch_bam_xeno = Channel.empty()
        if (params.run_split_fastq){
            SPLIT_LIB(
                ch_reads
            )
            ch_bam_xeno = SPLIT_LIB.out.bam_xeno
            // [ [meta], val(out_prefix), path/to/filtered.bam ]

        }else{
            SINGLE_LIB(
                ch_bam,
                ch_bam_host
            )
            ch_bam_xeno = SINGLE_LIB.out.bam_xeno
            // [ [meta], val(out_prefix), path/to/filtered.bam ]

        }
    }

    if(params.workflow == 'pdx' & params.save_fastq){
        BAM_TO_FASTQ(
            ch_bam_xeno
        )
            
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
                ch_bai,
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
                ch_bai,
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
                ch_bai,
                params.genome_fa,
                params.target_region
            )
            ch_hs_metrics = HS_METRICS.out.qc
            // [ [meta], path("*") ]
        }

        /*
        * pdx
        */
        if (params.workflow == 'pdx'){
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
            ch_bam_stat_host = SAMTOOLS_VIEW_HOST.out.data

            SAMTOOLS_VIEW_XENO(
                ch_bam_xeno
            )
            ch_bam_stat_xeno = SAMTOOLS_VIEW_XENO.out.data
        }


    }


    /*
    * MultiQC
    */
    ch_multiqc = Channel.empty()
    if (params.run_multiqc){
        if (params.workflow == 'pdx'){
            MULTIQC_PDX(
            ch_multiqc_config,
            ch_fastqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_cutadapt_js.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_fastqc_trimmed.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rseqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rnaseqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_star_log.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_star_log_host.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_bam_stat.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_bam_stat_host.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_bam_stat_xeno.map{it[1]}.flatten().collect().ifEmpty([])
            )
            ch_multiqc = MULTIQC_PDX.out.data
            //ch_multiqc.view()

        }else{
            MULTIQC(
            ch_multiqc_config,
            ch_fastqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_cutadapt_js.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_fastqc_trimmed.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rseqc.map{it[1]}.flatten().collect().ifEmpty([]),  
            ch_rnaseqc.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_hs_metrics.map{it[1]}.flatten().collect().ifEmpty([]),
            ch_star_log.map{it[2]}.flatten().collect().ifEmpty([]), 
            ch_counts.map{it[2]}.flatten().collect().ifEmpty([])
            )
            ch_multiqc = MULTIQC.out.data
            //ch_multiqc.view()
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
        // [ [meta], val(out_prefix), path("count.txt") ]

    }

    /*
    *   collect gene-level count matrix and run PCA
    */
    ch_gene_rds = Channel.empty()
    if (params.run_alignment & !params.only_fastq){
        GENERATE_GENE_COUNT_MATRIX(
            samplesheet, 
            ch_counts.map{it[2]}.collect().ifEmpty([]), 
            params.gene_txt,
            params.length_col,
            params.strand,
            params.workflow
        )
        ch_gene_rds = GENERATE_GENE_COUNT_MATRIX.out.rds
    }
    
    /*
    * differential genes
    */
    ch_de = Channel.empty()
    if (params.run_de){
        DIFFERENTIAL_EXPRESSION(
            samplesheet, 
            Channel.fromPath(params.comparison, checkIfExists: true), 
            ch_gene_rds, 
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2
        )
        ch_de = DIFFERENTIAL_EXPRESSION.out.rds
    }

    /*
    *   collect transcript-level count matrix and run PCA
    */
    ch_tx_rds = Channel.empty()
    if (params.run_alignment & params.run_salmon & !params.only_fastq){
        GENERATE_TRANSCRIPT_COUNT_MATRIX(
            samplesheet, 
            Channel.fromPath(params.comparison, checkIfExists: true), 
            ch_salmon.map{it[2]}.collect().ifEmpty([]), 
            params.tx_txt,
            "EffectiveLength"
        )
        ch_tx_rds = GENERATE_TRANSCRIPT_COUNT_MATRIX.out.rds
    }


    /*
    * differential transcripts
    */
    ch_dt = Channel.empty()
    if (params.run_dt & params.run_salmon){
        DIFFERENTIAL_TRANSCRIPTS(
            samplesheet, 
            Channel.fromPath(params.comparison, checkIfExists: true), 
            ch_tx_rds, 
            "EffectiveLength",
            params.fdr,
            params.fc,
            params.fdr2,
            params.fc2
        )
        ch_dt = DIFFERENTIAL_TRANSCRIPTS.out.rds
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
            ch_de.ifEmpty([]),
            ch_dt.ifEmpty([]),
            ch_report_rmd
        )
    }
}
